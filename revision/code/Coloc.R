## ================================
## 0. Set up environment & paths
## ================================
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(biomaRt)
library(coloc)
library(arrow)

base_dir  <- "/Users/adithyamadduri/Desktop/Projects/Patel_Lab/Revision/Coloc"
gwas_file <- file.path(base_dir, "GCST90027158_buildGRCh38.tsv")
gtex_dir  <- file.path(base_dir, "GTEx_Analysis_v10_eQTL_updated")
pred_dir  <- file.path(base_dir, "Final_Predictive_Genes_for_Coloc")

## Tissues we will use
# tissues <- c("Brain_Cortex", "Brain_Frontal_Cortex_BA9")
## Automatically detect all GTEx brain tissues present in the folder
all_files <- list.files(gtex_dir, full.names = TRUE)

tissues <- unique(
  sub("\\.v10.*", "", basename(all_files))        # remove .v10.* suffix
)

## Keep only brain tissues
tissues <- tissues[grepl("^Brain_", tissues)]

tissues

## Map tissue -> eQTL parquet file
eqtl_files <- setNames(
  vapply(tissues, function(t) {
    f <- list.files(
      gtex_dir,
      pattern = paste0("^", t, ".*eQTLs\\.signif_pairs\\.(parquet|txt\\.gz)$"),
      full.names = TRUE
    )
    if (length(f) == 0L) stop("No eQTL file found for tissue: ", t)
    f[1]
  }, character(1)),
  tissues
)

eqtl_files

# ----
## ================================
## 1. Load predictor gene lists
## ================================

## Assumes you saved:
##   predictive_genes_Ast.csv, predictive_genes_Mic.csv, ...
cell_types <- c("Ast", "Mic", "In", "Ex", "Oli", "Opc")

list_predictors <- lapply(cell_types, function(ct) {
  f <- file.path(pred_dir, paste0("predictive_genes_", ct, ".csv"))
  if (!file.exists(f)) {
    stop("Predictor file not found for ", ct, ": ", f)
  }
  df <- fread(f)
  
  ## <- only check for 'gene' and use that
  if (!"gene" %in% names(df)) {
    stop("Predictor file ", f, " must have a column named 'gene'")
  }
  
  df %>%
    transmute(
      cell_type   = ct,
      gene_symbol = toupper(trimws(gene))
    )
})

predictors <- bind_rows(list_predictors) %>% distinct()
print(head(predictors))

# -----
## ================================
## 2. Map symbols -> Ensembl GRCh38
## ================================
library(biomaRt)
library(dplyr)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  GRCh    = 38
)

gene_annot <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id",
                 "chromosome_name", "start_position", "end_position"),
  filters    = "hgnc_symbol",
  values     = unique(predictors$gene_symbol),
  mart       = mart
)

## Set clean column names explicitly
colnames(gene_annot) <- c(
  "gene_symbol",      # hgnc_symbol
  "ensembl_gene_id",  # ensembl_gene_id
  "chr",              # chromosome_name
  "start",            # start_position
  "end"               # end_position
)

## Keep autosomes + X/Y only
gene_annot <- gene_annot %>%
  filter(chr %in% c(as.character(1:22), "X", "Y"))

## Join back to predictors
predictors_annot <- predictors %>%
  inner_join(gene_annot, by = "gene_symbol")

nrow(predictors_annot)
head(predictors_annot)
# -------
## ================================
## 3. Load Belenguez GWAS summary
## ================================
gwas <- fread(gwas_file)

## Harmonize column names
gwas <- gwas %>%
  mutate(
    chr = as.character(chromosome),
    pos = as.integer(base_pair_location),
    snp_pos = paste0("chr", chr, ":", pos)
  )

## ---- Deduplicate GWAS by chr:pos to avoid duplicated SNPs ----
# Convert to data.table
setDT(gwas)

# Create snp_pos
gwas[, snp_pos := paste0("chr", chromosome, ":", base_pair_location)]

# Deduplicate keeping the lowest p-value
setorder(gwas, snp_pos, p_value)   # fast sort
gwas <- gwas[!duplicated(snp_pos)]

## verify
cat("GWAS snp_pos check:\n")
print(head(gwas$snp_pos))


## Basic sanity check on expected cols
stopifnot(all(c("chromosome", "base_pair_location",
                "effect_allele", "other_allele",
                "p_value") %in% names(gwas)))

## Derive beta if odds_ratio present
if ("odds_ratio" %in% names(gwas) && !("beta" %in% names(gwas))) {
  gwas[, beta := log(odds_ratio)]
}

if (!"standard_error" %in% names(gwas)) {
  stop("GWAS file missing 'standard_error' column; please check format.")
}

gwas <- gwas %>%
  filter(!is.na(chromosome),
         chromosome %in% as.character(1:22)) %>%
  mutate(
    chr  = as.character(chromosome),
    pos  = as.integer(base_pair_location),
    A1   = toupper(effect_allele),
    A2   = toupper(other_allele),
    maf  = if ("effect_allele_frequency" %in% names(.))
      pmin(effect_allele_frequency, 1 - effect_allele_frequency)
    else NA_real_,
    snp_pos = paste0("chr", chr, ":", pos)   ## <-- match key
  )

N_case    <- 111326
N_control <- 677663
N_total   <- N_case + N_control
prop_case <- N_case / N_total
# --------
## ================================
## 4. Load GTEx v10 eQTLs (parquet)
## ================================
library(arrow)
library(data.table)
library(dplyr)

eqtl_list <- list()
eqtl_N    <- numeric(length(tissues))
names(eqtl_N) <- tissues

for (t in tissues) {
  message("Loading eQTLs for ", t, " from ", eqtl_files[[t]])
  
  dt <- if (grepl("\\.parquet$", eqtl_files[[t]])) {
    as.data.table(arrow::read_parquet(eqtl_files[[t]]))
  } else {
    as.data.table(fread(eqtl_files[[t]]))
  }
  
  if (!all(c("gene_id", "variant_id", "slope", "slope_se") %in% names(dt))) {
    stop("eQTL file for ", t, " missing one of: gene_id, variant_id, slope, slope_se")
  }
  
  ## drop version from gene_id: ENSG00000123456.3 -> ENSG00000123456
  dt[, gene_id := gsub("\\.[0-9]+$", "", gene_id)]
  
  ## parse variant_id: chr_pos_ref_alt_b38
  tmp <- tstrsplit(dt$variant_id, "_", fixed=TRUE)
  
  dt[, chr := gsub("^chr","", tmp[[1]]) ]
  dt[, pos := as.integer(tmp[[2]]) ]
  dt[, ref := toupper(tmp[[3]]) ]
  dt[, alt := toupper(tmp[[4]]) ]
  
  dt[, snp_pos := paste0("chr", chr, ":", pos)]
  
  cat("GTEx snp_pos preview:\n")
  print(head(dt$snp_pos))
  
  ## ---- Deduplicate eQTL SNPs by chr:pos to avoid duplicated SNPs ----
  ## ---- Deduplicate eQTL SNPs by chr:pos (FAST + SAFE) ----
  setDT(dt)   # ensure dt is a data.table
  
  # fast sort
  setorder(dt, snp_pos, pval_nominal)
  
  # keep top row per snp_pos
  dt <- dt[!duplicated(snp_pos)]
  
  # MAF
  if ("af" %in% names(dt)) {
    dt[, maf := af]
  } else if (!"maf" %in% names(dt)) {
    dt[, maf := NA_real_]
  }
  
  # ensure ma_samples exists
  if (!"ma_samples" %in% names(dt)) {
    dt[, ma_samples := max(100, .N)]
  }
  
  eqtl_list[[t]] <- dt
  eqtl_N[t]      <- max(dt$ma_samples, na.rm = TRUE)
}

str(eqtl_list[[1]], max.level = 1)
eqtl_N

# --------
## ================================
## 5. Helper functions for coloc
## ================================
region_padding <- 250000L

make_region_limits <- function(chr, start, end, pad = region_padding) {
  list(
    chr   = as.character(chr),
    start = max(1L, as.integer(start) - pad),
    end   = as.integer(end) + pad
  )
}

make_gwas_dataset <- function(chr, start, end) {
  df <- gwas %>%
    filter(chr == !!chr,
           pos >= start,
           pos <= end)
  
  if (nrow(df) < 10L) return(NULL)
  
  list(
    snp      = df$snp_pos,
    beta     = df$beta,
    varbeta  = df$standard_error^2,
    MAF      = df$maf,
    N        = N_total,
    type     = "cc",
    s        = prop_case
  )
}

make_eqtl_dataset <- function(tissue, ensembl_gene_id, chr, start, end) {
  dt <- eqtl_list[[tissue]] %>%
    filter(gene_id == !!ensembl_gene_id,
           chr == !!chr,
           pos >= start,
           pos <= end)
  
  if (nrow(dt) < 10L) return(NULL)
  
  list(
    snp     = dt$snp_pos,
    beta    = dt$slope,
    varbeta = dt$slope_se^2,
    MAF     = dt$maf,
    N       = eqtl_N[[tissue]],
    type    = "quant"
  )
}
# --------
## ================================
## 6. Run coloc for each gene / CT / tissue
## ================================
results_list <- list()
idx <- 1L

for (tissue in tissues) {
  message("=== Tissue: ", tissue, " ===")
  
  for (i in seq_len(nrow(predictors_annot))) {
    row <- predictors_annot[i, ]
    message("  ", tissue, " | ", row$cell_type, " | ", row$gene_symbol)
    
    reg <- make_region_limits(row$chr, row$start, row$end)
    
    D1 <- make_gwas_dataset(reg$chr, reg$start, reg$end)
    if (is.null(D1)) {
      message("    Skipping: too few GWAS SNPs in region")
      next
    }
    
    D2 <- make_eqtl_dataset(tissue, row$ensembl_gene_id, reg$chr, reg$start, reg$end)
    if (is.null(D2)) {
      message("    Skipping: no eQTL SNPs for this gene/tissue in region")
      next
    }
    
    common <- intersect(D1$snp, D2$snp)
    if (length(common) < 10L) {
      message("    Skipping: <10 shared SNPs between GWAS and eQTL")
      next
    }
    
    keep1 <- D1$snp %in% common
    keep2 <- D2$snp %in% common
    
    D1_sub <- lapply(D1, function(x) if (length(x) == length(D1$snp)) x[keep1] else x)
    D2_sub <- lapply(D2, function(x) if (length(x) == length(D2$snp)) x[keep2] else x)
    
    coloc_res <- coloc.abf(D1_sub, D2_sub)
    
    summary_pp <- as.list(coloc_res$summary)
    summary_pp$tissue        <- tissue
    summary_pp$cell_type     <- row$cell_type
    summary_pp$gene_symbol   <- row$gene_symbol
    summary_pp$ensembl_gene  <- row$ensembl_gene_id
    summary_pp$chr           <- reg$chr
    summary_pp$region_start  <- reg$start
    summary_pp$region_end    <- reg$end
    summary_pp$nsnps         <- coloc_res$nsnps
    
    results_list[[idx]] <- summary_pp
    idx <- idx + 1L
  }
}

coloc_summary <- bind_rows(results_list)

out_file <- file.path(base_dir, "coloc_results_AD_GWAS_GTEx_v10_Cortex_and_BA9_all_brain.csv")
fwrite(coloc_summary, out_file)
message("Saved coloc summary to: ", out_file)

head(coloc_summary)

