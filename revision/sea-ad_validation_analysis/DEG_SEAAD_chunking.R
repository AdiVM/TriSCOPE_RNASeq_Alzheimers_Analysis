## ============================================================
## DEG PIPELINE FOR SEA-AD DATA (MTX-based, cell-level)
## ============================================================

library(dplyr)
library(Matrix)
library(edgeR)
library(RUVSeq)
library(lme4)
library(glmmTMB)
library(Seurat)

## -------------------------------
## 0. Command-line argument
## -------------------------------
args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 3) {
#   stop("Usage: Rscript run_DEG_SEAAD.R <celltype> <chunk_id> <n_chunks>")
# }

# celltype  <- args[1]
# chunk_id  <- as.integer(args[2])
# n_chunks  <- as.integer(args[3])

if (length(args) < 3 || length(args) > 4) {
  stop("Usage: Rscript run_DEG_SEAAD.R <celltype> <chunk_id> <n_chunks> [set_id]")
}

celltype  <- args[1]
chunk_id  <- as.integer(args[2])
n_chunks  <- as.integer(args[3])
set_id    <- if (length(args) == 4) as.integer(args[4]) else NULL

cat("Running DEG for cell type:", celltype, "\n")
cat("Chunk", chunk_id, "of", n_chunks, "\n")

## -------------------------------
## 1. Paths
## -------------------------------
base_dir <- "/n/scratch/users/a/adm808/SEAAD_DEG_input"

out_dir <- "/n/groups/patel/adithya/Alz_Outputs/Final_Outputs/seaad_DEG_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## -------------------------------
## 2. Resolve MTX file(s) per cell type
## -------------------------------
get_mtx_paths <- function(celltype) {
  switch(
    celltype,
    "Ast" = list(c(
      mtx       = file.path(base_dir, "Ast_counts.mtx"),
      barcodes  = file.path(base_dir, "Ast_barcodes.csv"),
      metadata  = file.path(base_dir, "Ast_metadata.csv")
    )),
    "Inh" = list(c(
      mtx       = file.path(base_dir, "Inh_counts.mtx"),
      barcodes  = file.path(base_dir, "Inh_barcodes.csv"),
      metadata  = file.path(base_dir, "Inh_metadata.csv")
    )),
    "Oli" = list(c(
      mtx       = file.path(base_dir, "Oli_counts.mtx"),
      barcodes  = file.path(base_dir, "Oli_barcodes.csv"),
      metadata  = file.path(base_dir, "Oli_metadata.csv")
    )),
    "Opc" = list(c(
      mtx       = file.path(base_dir, "Opc_counts.mtx"),
      barcodes  = file.path(base_dir, "Opc_barcodes.csv"),
      metadata  = file.path(base_dir, "Opc_metadata.csv")
    )),
    "Mic" = list(c(
      mtx       = file.path(base_dir, "Mic_counts.mtx"),
      barcodes  = file.path(base_dir, "Mic_barcodes.csv"),
      metadata  = file.path(base_dir, "Mic_metadata.csv")
    )),
    "Ex"  = list(
      c(
        mtx       = file.path(base_dir, "Ex_set1_counts.mtx"),
        barcodes  = file.path(base_dir, "Ex_set1_barcodes.csv"),
        metadata  = file.path(base_dir, "Ex_set1_metadata.csv")
      ),
      c(
        mtx       = file.path(base_dir, "Ex_set2_counts.mtx"),
        barcodes  = file.path(base_dir, "Ex_set2_barcodes.csv"),
        metadata  = file.path(base_dir, "Ex_set2_metadata.csv")
      ),
      c(
        mtx       = file.path(base_dir, "Ex_set3_counts.mtx"),
        barcodes  = file.path(base_dir, "Ex_set3_barcodes.csv"),
        metadata  = file.path(base_dir, "Ex_set3_metadata.csv")
      )
    ),
    stop("Unknown cell type: ", celltype)
  )
}

mtx_files <- get_mtx_paths(celltype)

## -------------------------------
## 3. Load MTX + barcodes + metadata (per set)
## -------------------------------
gene_path <- file.path(base_dir, "genes.csv")
gene_names <- read.csv(gene_path, stringsAsFactors = FALSE)$gene

load_set <- function(paths) {
  cat("Reading MTX:", paths[["mtx"]], "\n")
  mat <- readMM(paths[["mtx"]])  # cells x genes (from extraction script)
  mat <- t(mat)                  # transpose -> genes x cells (Seurat convention)
  mat <- as(mat, "CsparseMatrix")  # match original Seurat dgCMatrix class

  barcodes <- read.csv(paths[["barcodes"]], stringsAsFactors = FALSE)$barcode
  meta     <- read.csv(paths[["metadata"]], stringsAsFactors = FALSE,
                       check.names = FALSE)

  rownames(mat) <- gene_names
  colnames(mat) <- barcodes

  list(counts = mat, meta = meta)
}

if (celltype == "Ex") {
  if (!is.null(set_id)) {
    obj_list <- list(load_set(mtx_files[[set_id]]))
  } else {
    obj_list <- lapply(mtx_files, load_set)
  }
} else {
  obj_list <- lapply(mtx_files, load_set)
}

# for (obj_idx in seq_along(obj_list)) {

for (i in seq_along(obj_list)) {
  obj_idx <- if (!is.null(set_id)) set_id else i

  # obj <- obj_list[[obj_idx]]
  obj <- obj_list[[i]]

  cat("Processing object", obj_idx, "of", length(obj_list), "\n")

## Extract raw counts + metadata
expr_sparse <- obj$counts
meta <- obj$meta
meta$barcode <- as.character(meta$barcode)

stopifnot(all(colnames(expr_sparse) == meta$barcode))

cat("Loaded cells:", ncol(expr_sparse), "\n")
cat("Loaded genes:", nrow(expr_sparse), "\n")

# SEA-AD donor column is "Donor ID"; normalize to projid for code parity
names(meta)[names(meta) == "Donor ID"] <- "projid"
meta$projid <- as.character(meta$projid)

cat("Unique donors (projid):", length(unique(meta$projid)), "\n")

genes <- rownames(expr_sparse)

## -------------------------------
## 4. Clinical metadata
## -------------------------------
# SEA-AD metadata already contains AD status + donor covariates;
# no external clinical join required.

## -------------------------------
## 5. Gene filtering
## -------------------------------
expr_mic <- expr_sparse

## Remove zero-variance genes
rs  <- rowSums(expr_mic)
rs2 <- rowSums(expr_mic^2)
nc  <- ncol(expr_mic)
gene_var <- (rs2 / nc) - (rs / nc)^2

expr_mic_novar <- expr_mic[gene_var > 0, , drop = FALSE]
cat("After removing zero-variance genes:", nrow(expr_mic_novar), "\n")

## 20% detection filter (unchanged)
prop_detect <- rowSums(expr_mic_novar > 0) / ncol(expr_mic_novar)
gene_filter <- prop_detect > 0.20
expr_mic_filt <- expr_mic_novar[gene_filter, , drop = FALSE]

cat("After 20% detection filter:", nrow(expr_mic_filt), "\n")

expr_mic_final <- expr_mic_filt

## -------------------------------
## 6. Metadata prep
## -------------------------------
meta$projid <- as.factor(meta$projid)

meta$AD <- as.factor(as.integer(meta$alzheimers_or_control))

meta <- meta[match(colnames(expr_mic_final), meta$barcode), ]

## -------------------------------
## 7. RUVSeq (donor-level)
## -------------------------------
cell_tags  <- colnames(expr_mic_final)
projid_vec <- meta$projid

f_projid   <- factor(projid_vec)
donor_ids  <- as.integer(f_projid)
donor_lvls <- levels(f_projid)
ndonors    <- length(donor_lvls)

M <- sparseMatrix(
  i    = seq_along(donor_ids),
  j    = donor_ids,
  x    = 1,
  dims = c(length(donor_ids), ndonors)
)

data_ind <- expr_mic_final %*% M
data_ind <- as.matrix(data_ind)
colnames(data_ind) <- donor_lvls

ad_by_donor <- tapply(
  as.numeric(as.character(meta$AD)),
  meta$projid,
  function(x) unique(x[!is.na(x)])
)
ad_by_donor <- ad_by_donor[colnames(data_ind)]
ad_factor <- factor(ad_by_donor)

d_e <- DGEList(counts = data_ind, genes = rownames(data_ind))
keep <- rowSums(cpm(d_e) > 1) >= 3
d_e <- d_e[keep, , keep.lib.sizes = FALSE]

d_e <- calcNormFactors(d_e, method = "TMM")
design_ind <- model.matrix(~ ad_factor)

d_e <- estimateGLMCommonDisp(d_e, design_ind)
d_e <- estimateGLMTagwiseDisp(d_e, design_ind)

fit1 <- glmFit(d_e, design_ind)
res1 <- residuals(fit1, type = "deviance")

ruvn <- 10

set_ruv <- newSeqExpressionSet(
  as.matrix(round(d_e$counts)),
  phenoData = data.frame(
    row.names = colnames(d_e$counts),
    ad = ad_factor[colnames(d_e$counts)]
  )
)

set_ruv <- RUVr(
  set_ruv,
  cIdx = rownames(d_e$counts),
  k = ruvn,
  res = res1
)

W <- pData(set_ruv)[, grep("^W_", colnames(pData(set_ruv))), drop = FALSE]
ruv_df <- data.frame(projid = rownames(W), W)

meta <- meta %>%
  left_join(ruv_df, by = "projid")

w_cols <- grep("^W_", colnames(meta), value = TRUE)

## -------------------------------
## 8. Cell-level Poisson GLMM DE
## -------------------------------
results  <- list()
libsize  <- colSums(expr_mic_final)
genes_ct <- rownames(expr_mic_final)
# -------------------------------
# Gene chunking
# -------------------------------
ngenes <- length(genes_ct)

chunk_size <- ceiling(ngenes / n_chunks)
start_idx  <- (chunk_id - 1) * chunk_size + 1
end_idx    <- min(chunk_id * chunk_size, ngenes)

genes_ct <- genes_ct[start_idx:end_idx]

cat("Processing genes", start_idx, "to", end_idx, "out of", ngenes, "\n")


counter <- 0
for (gene in genes_ct) {
  counter <- counter + 1
  if (counter %% 1000 == 0) cat(counter, "genes processed\n")

  y <- as.numeric(expr_mic_final[gene, ])

  dt <- meta
  dt$y <- y
  dt$libsize <- libsize[dt$barcode]

  to_keep <- !is.na(dt$AD) & !is.na(dt$projid) & !is.na(y)
  for (wc in w_cols) to_keep <- to_keep & !is.na(dt[[wc]])

  dt2 <- dt[to_keep, ]
  dt2$projid <- as.factor(dt2$projid)
  dt2$AD <- as.factor(dt2$AD)

  ruv_term <- if (length(w_cols) > 0) paste(w_cols, collapse = " + ") else NULL
  formula_str <- if (!is.null(ruv_term)) {
    paste0("y ~ AD + ", ruv_term, " + offset(log(libsize)) + (1 | projid)")
  } else {
    "y ~ AD + offset(log(libsize)) + (1 | projid)"
  }

  fit_simple <- tryCatch(
    glmer(
      as.formula(formula_str),
      data = dt2,
      family = poisson(link = "log"),
      nAGQ = 10,
      control = glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) NULL
  )

  if (!is.null(fit_simple)) {
    coefs <- summary(fit_simple)$coefficients
    if ("AD1" %in% rownames(coefs)) {
      results[[gene]] <- data.frame(
        gene        = gene,
        estimate_AD = coefs["AD1", "Estimate"],
        se_AD       = coefs["AD1", "Std. Error"],
        pval_AD     = coefs["AD1", "Pr(>|z|)"],
        n_cells     = nrow(dt2)
      )
    }
  }
}

final_results <- do.call(rbind, results)
final_results$p_adj <- p.adjust(final_results$pval_AD, method = "fdr")
final_results$log2FC <- final_results$estimate_AD / log(2)
final_results$DEG <- ifelse(
  final_results$p_adj < 0.05 & abs(final_results$log2FC) > 0.25,
  "True", "False"
)

# out_path <- file.path(out_dir, paste0("poisson_DE_results_SEAAD_", celltype, ".csv"))

out_path <- file.path(
  out_dir,
  paste0("poisson_DE_results_SEAAD_", celltype,
       "_set", obj_idx,
       "_chunk", chunk_id, "_of_", n_chunks, ".csv")
)

write.csv(final_results, out_path, row.names = FALSE)

cat("DE results written to:", out_path, "\n")

}
