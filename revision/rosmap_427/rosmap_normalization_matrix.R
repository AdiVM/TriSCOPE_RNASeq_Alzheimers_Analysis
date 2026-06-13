## ============================================================
## rosmap_normalization_matrix.R
##
## Builds, for the Mathys PFC427 data, parquet files in the SAME
## schema the SEA-AD predictive pipeline consumes, so the seaad
## ein / non-ein Multirun scripts run on 427 with only path changes.
##
## Source object: the SAME per-cell-type Seurat RDS that
## DEG_PFC_Maximal_chunking.R used (raw counts = obj@assays$RNA@counts).
##
## Normalization: IDENTICAL to sea_ad_normalization_matrix.py
##     CPM    = counts / total_counts_per_cell * 1e6
##     X      = log2(CPM + 1)            (float32)
##
## Outputs (under OUT_DIR):
##   ROSMAP427_Matrix_<ct>/part-XXX.parquet   (cells x genes, 'TAG' column, float32)
##   ROSMAP427_Meta_<ct>.parquet              (per-cell metadata, SEA-AD column names)
##   ROSMAP427_CellMetadata.parquet           (all cell types concatenated)
##
## Cell types: Ex, Inh, Ast, Oli, Mic, Opc   (SEA-AD codes)
## Ex is processed set-by-set (never merged at cell level) to bound RAM.
##
## Run via submit_rosmap_normalization.sh (module R/4.4.2 + ~/R/library).
## ============================================================

suppressMessages({
  library(Seurat)
  library(Matrix)
  library(arrow)
})

## ── Paths (same data the DEG used) ─────────────────────────────────────────────
base_dir  <- "/n/groups/patel/randy/single_cell/syn52293417/Data/Gene Expression/Gene Expression (snRNAseq)/Gene Expression (snRNAseq - 10x)/Gene Expression (snRNAseq - 10x) processed"
clin_path <- "/n/groups/patel/randy/single_cell/syn18485175/metadata/ROSMAP_clinical.csv"

OUT_DIR  <- "/n/groups/patel/adithya/ROSMAP427_Outputs"
META_OUT <- file.path(OUT_DIR, "ROSMAP427_CellMetadata.parquet")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## Cells per sub-chunk when densifying (genes x SUBCHUNK float64 transient).
## ~33k genes x 20k cells x 8 bytes ~= 5.3 GB (x2 transient for transpose).
SUBCHUNK_SIZE <- 20000L

## ── Cell-type -> RDS file(s).  SEA-AD codes (Inh, Mic) as used by the ein. ──────
rds_map <- list(
  Ast = c("Astrocytes.rds"),
  Inh = c("Inhibitory_neurons.rds"),
  Oli = c("Oligodendrocytes.rds"),
  Opc = c("OPCs.rds"),
  Mic = c("Immune_cells.rds"),
  Ex  = c("Excitatory_neurons_set1.rds",
          "Excitatory_neurons_set2.rds",
          "Excitatory_neurons_set3.rds")
)
# Process small cell types first (fail-fast end-to-end validation + early donor
# count), with the large Ex object last. Order does not affect output.
CELL_TYPES <- c("Mic", "Opc", "Ast", "Oli", "Inh", "Ex")

## ── Clinical metadata (donor-level) ────────────────────────────────────────────
clin <- read.csv(clin_path, stringsAsFactors = FALSE)
clin$projid <- as.character(clin$projid)
# Keep first row per donor (clinical is donor-level; guard against dups)
clin <- clin[!duplicated(clin$projid), ]
rownames(clin) <- clin$projid

## Age at Death: ROSMAP codes oldest as "90+" -> 90 (per original ROSMAP pipeline)
norm_age <- function(x) {
  x <- as.character(x)
  x[x %in% c("90+", "90 +")] <- "90"
  suppressWarnings(as.numeric(x))
}

## Build the SEA-AD-schema metadata frame for a set of cell barcodes + projids.
build_meta <- function(barcodes, projid, ct) {
  projid <- as.character(projid)
  cl <- clin[projid, , drop = FALSE]
  dcfdx <- suppressWarnings(as.numeric(cl$dcfdx_lv))
  data.frame(
    TAG                    = barcodes,
    `Donor ID`             = projid,
    Sex                    = ifelse(cl$msex == 1, "Male",
                              ifelse(cl$msex == 0, "Female", NA_character_)),
    `Age at Death`         = norm_age(cl$age_death),
    PMI                    = suppressWarnings(as.numeric(cl$pmi)),
    `APOE Genotype`        = as.character(cl$apoe_genotype),
    Education              = suppressWarnings(as.numeric(cl$educ)),
    alzheimers_or_control  = as.integer(dcfdx %in% c(4, 5)),  # NA -> FALSE -> 0
    broad_cell_type        = ct,
    check.names            = FALSE,
    stringsAsFactors       = FALSE
  )
}

## CPM -> log2(x+1), IDENTICAL math to sea_ad_normalization_matrix.py.
## Input: sparse genes x cells (sub-chunk). Output: dense cells x genes (double).
normalize_cpm_log2_t <- function(counts_sub) {
  total <- Matrix::colSums(counts_sub)          # per-cell total counts
  total[total == 0] <- 1                         # guard divide-by-zero
  # scale each cell (column) by 1e6/total  -> CPM, still sparse
  cpm <- counts_sub %*% Matrix::Diagonal(x = 1e6 / total)
  dense <- as.matrix(cpm)                         # genes x cells
  dense <- log2(dense + 1)                        # log2(CPM+1)
  t(dense)                                         # cells x genes
}

meta_parts <- list()

for (ct in CELL_TYPES) {
  cat("\n==================== cell type:", ct, "====================\n")
  mat_dir   <- file.path(OUT_DIR, paste0("ROSMAP427_Matrix_", ct))
  meta_file <- file.path(OUT_DIR, paste0("ROSMAP427_Meta_", ct, ".parquet"))

  if (dir.exists(mat_dir) && length(list.files(mat_dir, pattern = "\\.parquet$")) > 0 &&
      file.exists(meta_file)) {
    cat("  Already done, loading meta part and skipping:", ct, "\n")
    meta_parts[[ct]] <- as.data.frame(arrow::read_parquet(meta_file))
    next
  }
  dir.create(mat_dir, recursive = TRUE, showWarnings = FALSE)

  files     <- rds_map[[ct]]
  ref_genes <- NULL          # gene order fixed by first file of this cell type
  schema    <- NULL
  ct_meta   <- list()
  part_idx  <- 0L

  for (fi in seq_along(files)) {
    f <- file.path(base_dir, files[fi])
    cat("\n  --- loading:", files[fi], "---\n")
    obj  <- readRDS(f)
    cnt  <- obj@assays$RNA@counts          # genes x cells, sparse, raw integer counts
    bc   <- colnames(cnt)
    proj <- as.character(obj@meta.data$projid)
    stopifnot(all(bc == rownames(obj@meta.data)))
    rm(obj); gc()

    genes <- rownames(cnt)
    if (is.null(ref_genes)) {
      ref_genes <- genes
      # float32 schema: TAG (string) + one float32 field per gene
      fields <- c(list(arrow::field("TAG", arrow::utf8())),
                  lapply(ref_genes, function(g) arrow::field(g, arrow::float32())))
      schema <- do.call(arrow::schema, fields)
    } else {
      # within a cell type all files must share the gene set/order
      if (!identical(genes, ref_genes)) {
        stopifnot(setequal(genes, ref_genes))
        cat("    (reordering genes to match reference order)\n")
        cnt <- cnt[ref_genes, , drop = FALSE]
      }
    }

    cat("    cells:", ncol(cnt), " genes:", nrow(cnt),
        " donors:", length(unique(proj)), "\n")

    ## metadata for these cells
    ct_meta[[fi]] <- build_meta(bc, proj, ct)

    ## stream normalize in sub-chunks
    ncells <- ncol(cnt)
    nsub   <- ceiling(ncells / SUBCHUNK_SIZE)
    for (s in seq_len(nsub)) {
      a <- (s - 1) * SUBCHUNK_SIZE + 1
      b <- min(s * SUBCHUNK_SIZE, ncells)
      sub <- cnt[, a:b, drop = FALSE]
      norm_t <- normalize_cpm_log2_t(sub)          # cells x genes (double)
      df <- as.data.frame(norm_t, check.names = FALSE)
      colnames(df) <- ref_genes
      df <- cbind(TAG = bc[a:b], df)
      tbl <- arrow::as_arrow_table(df)$cast(schema) # -> float32 + utf8
      part_path <- file.path(mat_dir, sprintf("part-%04d.parquet", part_idx))
      arrow::write_parquet(tbl, part_path, compression = "snappy")
      part_idx <- part_idx + 1L
      cat("    wrote", basename(part_path), "(cells", a, "-", b, ")\n")
      rm(sub, norm_t, df, tbl); gc()
    }
    rm(cnt); gc()
  }

  ct_meta_df <- do.call(rbind, ct_meta)
  stopifnot(!any(duplicated(ct_meta_df$TAG)))      # barcodes unique within cell type
  arrow::write_parquet(arrow::as_arrow_table(ct_meta_df), meta_file)
  meta_parts[[ct]] <- ct_meta_df
  cat("  ", ct, "done:", nrow(ct_meta_df), "cells,",
      length(unique(ct_meta_df$`Donor ID`)), "donors\n")
}

## ── Assemble full CellMetadata ─────────────────────────────────────────────────
cat("\nAssembling full CellMetadata parquet...\n")
full_meta <- do.call(rbind, meta_parts[CELL_TYPES])
# TAGs are unique WITHIN each cell type (asserted per-type above), which is all
# the modeling needs (it filters to one cell type before every join/merge). The
# raw 10x barcode string can repeat ACROSS cell types; that is harmless here, so
# warn rather than fail if the global set is not unique.
ndup <- sum(duplicated(full_meta$TAG))
if (ndup > 0) cat("NOTE:", ndup, "TAG strings repeat across cell types",
                  "(harmless; joins are per-cell-type).\n")
arrow::write_parquet(arrow::as_arrow_table(full_meta), META_OUT)

cat("\n=====================================================\n")
cat("DONE.\n")
cat("  Metadata:", META_OUT, "  (", nrow(full_meta), "cells )\n")
cat("  Donors (overall):", length(unique(full_meta$`Donor ID`)), "\n")
cat("  AD (1):",  sum(tapply(full_meta$alzheimers_or_control, full_meta$`Donor ID`, function(x) x[1])),
    " donors;  Control (0):",
    sum(tapply(full_meta$alzheimers_or_control, full_meta$`Donor ID`, function(x) x[1]) == 0), "\n")
print(table(full_meta$broad_cell_type))
cat("=====================================================\n")
