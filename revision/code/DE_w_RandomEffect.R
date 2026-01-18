## ============================================================
## DEG PIPELINE FOR NEW PFC DATA (major_cell_type)
## ============================================================

library(arrow)
library(dplyr)
library(Matrix)
library(edgeR)
library(RUVSeq)
library(lme4)
library(glmmTMB)

## -------------------------------
## 0. Paths
## -------------------------------
cells_path   <- "/Users/adithyamadduri/Downloads/PFC_cells.tsv"
genes_path   <- "/Users/adithyamadduri/Downloads/PFC_genes.tsv"
sparse_path  <- "/Users/adithyamadduri/Downloads/PFC_matrix_sparse.parquet"
meta_path    <- "/Users/adithyamadduri/Downloads/PFC_metadata_with_clinical.csv"

out_dir <- "/Users/adithyamadduri/Desktop/Projects/Patel_Lab/Revision/New_Data_DEG_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## -------------------------------
## 1. Load barcodes and genes
## -------------------------------
barcodes <- read.delim(
  cells_path,
  header = FALSE,
  stringsAsFactors = FALSE
)[, 1]

genes <- read.delim(
  genes_path,
  header = FALSE,
  stringsAsFactors = FALSE
)[, 1]

cat("Loaded barcodes:", length(barcodes), "\n")
cat("Loaded genes:", length(genes), "\n")

## -------------------------------
## 2. Load sparse triplet matrix
## -------------------------------
tbl_sparse <- read_parquet(sparse_path)
df_sparse  <- as.data.frame(tbl_sparse)

## Expect columns: cell, gene, value
## cell and gene are 0-based indices into barcodes / genes
if (!all(c("cell", "gene", "value") %in% colnames(df_sparse))) {
  stop("Expected columns 'cell', 'gene', 'value' in PFC_matrix_sparse.parquet")
}

cat("Sparse matrix loaded:", nrow(df_sparse), "non-zero entries\n")

## Attach barcode to each entry (for easier per-cell-type filtering)
df_sparse$barcode <- barcodes[df_sparse$cell + 1L]

## -------------------------------
## 3. Load metadata with clinical info
## -------------------------------
meta <- read.csv(meta_path, stringsAsFactors = FALSE)

## Expect at least: barcode, projid, major_cell_type, dcfdx_lv, ceradsc
req_cols <- c("barcode", "projid", "major_cell_type", "dcfdx_lv", "ceradsc")
missing_cols <- setdiff(req_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop(paste("Missing columns in metadata:", paste(missing_cols, collapse = ", ")))
}

## Keep only cells present in expression
meta <- meta %>%
  filter(barcode %in% barcodes)

cat("Metadata rows after restricting to barcodes in expression:", nrow(meta), "\n")

## -------------------------------
## 4. Define cell types for loop
## -------------------------------
## Based on your counts: Exc, Oli, Inh, Ast, OPC, Mic
cell_types <- c("Exc")
# "Mic", "Ast", "Oli", "Inh", "OPC", 

## ============================================================
## MAIN LOOP OVER CELL TYPES
## ============================================================

for (celltype in cell_types) {
  cat("\n=============================\n")
  cat("Processing cell type:", celltype, "\n")
  cat("=============================\n")
  
  ## --- Subset metadata for this cell type ---
  meta_mic <- meta %>%
    filter(major_cell_type == celltype)
  
  if (nrow(meta_mic) == 0) {
    cat("No cells found for cell type", celltype, "- skipping.\n")
    next
  }
  
  ## Ensure barcodes are valid
  meta_mic <- meta_mic %>%
    filter(barcode %in% barcodes)
  
  if (nrow(meta_mic) == 0) {
    cat("No matching barcodes in expression for cell type", celltype, "- skipping.\n")
    next
  }
  
  ## Row names = barcode (analogous to TAG before)
  rownames(meta_mic) <- meta_mic$barcode
  
  ## List of barcodes for this cell type (used as expression columns)
  cell_barcodes <- meta_mic$barcode
  
  cat("Number of cells in", celltype, ":", length(cell_barcodes), "\n")
  
  ## --- Subset sparse triplets to these barcodes only (per-cell-type matrix) ---
  df_ct <- df_sparse[df_sparse$barcode %in% cell_barcodes, ]
  
  if (nrow(df_ct) == 0) {
    cat("No non-zero counts for cell type", celltype, "- skipping.\n")
    next
  }
  
  ## Map barcodes -> column indices 1..ncells for this cell type
  ## Columns will be ordered as meta_mic$barcode
  cell_index_map <- setNames(seq_along(cell_barcodes), cell_barcodes)
  col_index <- cell_index_map[df_ct$barcode]
  
  ## Row indices: gene index (0-based) + 1
  row_index <- df_ct$gene + 1L
  
  ## Safety checks
  if (any(is.na(col_index))) {
    stop("NA in col_index mapping for cell type: ", celltype)
  }
  if (max(row_index) > length(genes)) {
    stop("Row index exceeds gene vector length for cell type: ", celltype)
  }
  
  ## --- Construct sparse gene x cell matrix for this cell type ---
  ngenes <- length(genes)
  ncells <- length(cell_barcodes)
  
  expr_sparse <- sparseMatrix(
    i    = row_index,
    j    = col_index,
    x    = df_ct$value,
    dims = c(ngenes, ncells),
    dimnames = list(genes, cell_barcodes)
  )
  
  ## Convert to dense matrix for downstream steps (as in original pipeline)
  # expr_mic <- as.matrix(expr_sparse)
  
  # cat("Initial gene count:", nrow(expr_mic), "\n")
  
  ## -------------------------------
  ## 1. Remove zero variance genes
  ## -------------------------------
  # gene_var <- apply(expr_mic, 1, var)
  expr_mic <- expr_sparse   # keep sparse
  
  # compute variance in sparse matrix
  # var(x) = mean(x^2) - mean(x)^2
  rs <- rowSums(expr_mic)
  rs2 <- rowSums(expr_mic^2)
  nc <- ncol(expr_mic)
  gene_var <- (rs2 / nc) - (rs / nc)^2
  
  expr_mic_novar <- expr_mic[gene_var > 0, , drop = FALSE]
  cat("After removing zero-variance genes:", nrow(expr_mic_novar), "\n")
  
  ## -------------------------------
  ## 2. CPM filter
  ## -------------------------------
  #lib_size_millions <- colSums(expr_mic_novar) / 1e6
  #expr_cpm <- sweep(expr_mic_novar, 2, lib_size_millions, FUN = "/")
  #gene_filter <- rowSums(expr_cpm > 1) >= 10
  #expr_mic_filt <- expr_mic_novar[gene_filter, , drop = FALSE]
  #cat("After CPM filter:", nrow(expr_mic_filt), "\n")
  
  ## -------------------------------
  ## 2. 20% detection filter (Mathys-style)
  ## -------------------------------
  # prop_detect <- rowMeans(expr_mic_novar > 0)
  prop_detect <- rowSums(expr_mic_novar > 0) / ncol(expr_mic_novar)
  gene_filter <- prop_detect > 0.20     # gene expressed in >20% of cells
  
  expr_mic_filt <- expr_mic_novar[gene_filter, , drop = FALSE]
  # expr_mic_filt <- expr_mic_novar[gene_filter, , drop = FALSE]
  cat("After 20% detection filter:", nrow(expr_mic_filt), "\n")
  
  ## -------------------------------
  ## 3. Metadata prep
  ## -------------------------------
  meta_mic$projid   <- as.factor(meta_mic$projid)
  meta_mic$dcfdx_lv <- as.numeric(meta_mic$dcfdx_lv)
  
  ## AD status: 1 if dcfdx_lv in {4,5}, else 0
  meta_mic$AD <- as.factor(ifelse(meta_mic$dcfdx_lv %in% c(4, 5), 1, 0))
  
  ## Align metadata to columns of expr_mic_filt
  meta_mic <- meta_mic[match(colnames(expr_mic_filt), meta_mic$barcode), ]
  
  expr_mic_final <- expr_mic_filt
  cat("Using all remaining genes:", nrow(expr_mic_final), "\n")
  
  ## =====================================================================
  ## 3b. RUVSeq step (donor-level counts -> residuals -> RUVr covariates)
  ## =====================================================================
  ## Build donor-level count matrix data_ind: genes x donors
  cell_tags  <- colnames(expr_mic_final)  ## barcodes
  projid_vec <- meta_mic$projid[match(cell_tags, meta_mic$barcode)]
  
  cell_index_list <- split(seq_along(cell_tags), projid_vec)
  
  data_ind <- sapply(cell_index_list, function(idx) {
    rowSums(expr_mic_final[, idx, drop = FALSE])
  })
  
  data_ind <- as.matrix(data_ind)
  colnames(data_ind) <- names(cell_index_list)
  
  ## =====================================================================
  ## 3b. RUVSeq step (donor-level counts -> residuals -> RUVr covariates)
  ## =====================================================================
  ## Build donor-level count matrix data_ind: genes x donors
  #cell_tags  <- colnames(expr_mic_final)  ## barcodes
  #projid_vec <- meta_mic$projid[match(cell_tags, meta_mic$barcode)]
  
  ## projid_vec may be a factor; enforce a clean factor
  #f_projid   <- factor(projid_vec)
  #donor_ids  <- as.integer(f_projid)           # 1..ndonors
  #donor_lvls <- levels(f_projid)
  #ndonors    <- length(donor_lvls)
  
  ## Sparse cell x donor matrix: cells (rows) -> donors (cols)
  #M <- sparseMatrix(
  #i    = seq_along(donor_ids),
  #j    = donor_ids,
  #x    = 1,
  #dims = c(length(donor_ids), ndonors)
  #)
  
  ## =====================================================================
  ## 3b. RUVSeq step (donor-level counts -> residuals -> RUVr covariates)
  ## =====================================================================
  ## Build donor-level count matrix data_ind: genes x donors
  cell_tags  <- colnames(expr_mic_final)  ## barcodes
  projid_vec <- meta_mic$projid[match(cell_tags, meta_mic$barcode)]
  
  ## projid_vec may be a factor; enforce a clean factor
  f_projid   <- factor(projid_vec)
  donor_ids  <- as.integer(f_projid)           # 1..ndonors
  donor_lvls <- levels(f_projid)
  ndonors    <- length(donor_lvls)
  
  ## Sparse cell x donor matrix: cells (rows) -> donors (cols)
  M <- sparseMatrix(
    i    = seq_along(donor_ids),
    j    = donor_ids,
    x    = 1,
    dims = c(length(donor_ids), ndonors)
  )
  
  ## genes x cells  %*%  cells x donors  =  genes x donors
  data_ind <- expr_mic_final %*% M
  data_ind <- as.matrix(data_ind)              # small enough: genes x donors
  colnames(data_ind) <- donor_lvls
  
  ## Donor-level AD (0/1) aligned to columns of data_ind
  ad_by_donor <- tapply(
    as.numeric(as.character(meta_mic$AD)),
    meta_mic$projid,
    function(x) unique(x[!is.na(x)])
  )
  ad_by_donor <- ad_by_donor[colnames(data_ind)]
  ad_factor <- factor(ad_by_donor)
  
  ## edgeR GLM on donor-level counts
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
  
  ## Build SeqExpressionSet for RUVr
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
  ruv_df <- data.frame(
    projid = rownames(W),
    W,
    row.names = NULL
  )
  
  ## Merge RUV covariates back to cell-level metadata
  meta_mic <- meta_mic %>%
    left_join(ruv_df, by = "projid")
  
  ## Identify RUV covariate columns
  w_cols <- grep("^W_", colnames(meta_mic), value = TRUE)
  
  ## -------------------------------
  ## 4. Run DE per gene
  ## -------------------------------
  results  <- list()
  libsize  <- colSums(expr_mic_final)
  genes_ct <- rownames(expr_mic_final)
  
  counter <- 0
  for (gene in genes_ct) {
    counter <- counter + 1
    if (counter %% 1000 == 0) cat(counter, "genes processed\n")
    
    y <- as.numeric(expr_mic_final[gene, ])
    
    dt <- meta_mic
    dt$y <- y
    dt$libsize <- libsize[dt$barcode]
    
    ## Only require AD, projid, y, library size, and non-missing RUV covariates
    to_keep <- !is.na(dt$AD) & !is.na(dt$projid) & !is.na(y)
    if (length(w_cols) > 0) {
      for (wc in w_cols) {
        to_keep <- to_keep & !is.na(dt[[wc]])
      }
    }
    dt2 <- dt[to_keep, ]
    
    dt2$projid <- as.factor(dt2$projid)
    dt2$AD     <- as.factor(dt2$AD)
    
    ## Build formula with RUV covariates
    if (length(w_cols) > 0) {
      ruv_term <- paste(w_cols, collapse = " + ")
      formula_str <- paste0(
        "y ~ AD + ", ruv_term,
        " + offset(log(libsize)) + (1 | projid)"
      )
    } else {
      formula_str <- "y ~ AD + offset(log(libsize)) + (1 | projid)"
    }
    
    ## Poisson GLMM with donor random intercept + RUV covariates
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
      summ  <- summary(fit_simple)
      coefs <- summ$coefficients
      
      if (!("AD1" %in% rownames(coefs))) {
        next
      }
      
      res <- data.frame(
        gene        = gene,
        estimate_AD = coefs["AD1", "Estimate"],
        se_AD       = coefs["AD1", "Std. Error"],
        pval_AD     = coefs["AD1", "Pr(>|z|)"],
        n_cells     = nrow(dt2)
      )
      results[[gene]] <- res
    }
  }
  
  if (length(results) == 0) {
    cat("No successfully fitted genes for cell type", celltype, "\n")
    next
  }
  
  final_results <- do.call(rbind, results)
  final_results$p_adj <- p.adjust(final_results$pval_AD, method = "fdr")
  final_results$log2FC <- final_results$estimate_AD / log(2)
  final_results$DEG <- ifelse(
    (final_results$p_adj < 0.05) & (abs(final_results$log2FC) > 0.25),
    "True",
    "False"
  )
  
  out_path <- file.path(
    out_dir,
    sprintf("poisson_DE_results_PFC_%s.csv", celltype)
  )
  write.csv(final_results, out_path, row.names = FALSE)
  cat(sprintf("DE results written to %s\n", out_path))
}