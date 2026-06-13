## ============================================================
## inspect_rosmap427.R
## Read-only inspection of the PFC427 per-cell-type Seurat RDS
## objects used by DEG_PFC_Maximal_chunking.R. Confirms, per cell
## type: n cells, n genes, n unique donors (projid), barcode format,
## gene-name format, and (for Ex) whether the 3 sets share an
## identical gene set. Loads one object at a time with gc() between
## to stay within memory. DOES NOT WRITE ANY DATA FILES.
## ============================================================

suppressMessages({
  library(Seurat)
  library(Matrix)
})

base_dir <- "/n/groups/patel/randy/single_cell/syn52293417/Data/Gene Expression/Gene Expression (snRNAseq)/Gene Expression (snRNAseq - 10x)/Gene Expression (snRNAseq - 10x) processed"
clin_path <- "/n/groups/patel/randy/single_cell/syn18485175/metadata/ROSMAP_clinical.csv"

# celltype code -> rds file(s)   (same mapping as the DEG script; Inh uses "Inh")
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

clin <- read.csv(clin_path, stringsAsFactors = FALSE)
clin$projid <- as.character(clin$projid)
cat("Clinical rows:", nrow(clin), " unique projid:", length(unique(clin$projid)), "\n")

all_donors <- character(0)

for (ct in names(rds_map)) {
  cat("\n==================== ", ct, " ====================\n")
  files <- rds_map[[ct]]
  ct_donors <- character(0)
  ex_gene_sets <- list()

  for (fi in seq_along(files)) {
    f <- file.path(base_dir, files[fi])
    cat("\n--- file:", files[fi], "---\n")
    obj <- readRDS(f)

    cnt <- obj@assays$RNA@counts
    meta <- obj@meta.data

    cat("class(counts):", class(cnt)[1],
        " | dim (genes x cells):", nrow(cnt), "x", ncol(cnt), "\n")
    cat("counts integer-valued?:",
        isTRUE(all(cnt@x == round(cnt@x))), "\n")
    cat("meta.data columns:\n"); print(colnames(meta))
    cat("projid in meta?:", "projid" %in% colnames(meta), "\n")

    proj <- as.character(meta$projid)
    cat("n cells:", ncol(cnt),
        " | n unique projid (this file):", length(unique(proj)), "\n")

    cat("barcode (rownames meta) head:\n"); print(head(rownames(meta), 3))
    cat("counts colnames == meta barcodes?:",
        isTRUE(all(colnames(cnt) == rownames(meta))), "\n")
    cat("gene (rownames counts) head:\n"); print(head(rownames(cnt), 5))

    # clinical-join coverage
    in_clin <- unique(proj) %in% clin$projid
    cat("donors covered by clinical:", sum(in_clin), "/", length(in_clin), "\n")

    ct_donors <- union(ct_donors, unique(proj))
    if (ct == "Ex") ex_gene_sets[[fi]] <- rownames(cnt)

    rm(obj, cnt, meta); gc()
  }

  cat("\n>>>", ct, "TOTAL unique donors across file(s):",
      length(ct_donors), "\n")
  all_donors <- union(all_donors, ct_donors)

  if (ct == "Ex" && length(ex_gene_sets) == 3) {
    same12 <- identical(ex_gene_sets[[1]], ex_gene_sets[[2]])
    same13 <- identical(ex_gene_sets[[1]], ex_gene_sets[[3]])
    cat(">>> Ex gene rownames identical across sets? set1==set2:",
        same12, " set1==set3:", same13, "\n")
    cat(">>> Ex gene counts:",
        sapply(ex_gene_sets, length), "\n")
  }
}

cat("\n==================================================\n")
cat("GRAND TOTAL unique donors across ALL cell types:",
    length(all_donors), "\n")
cat("==================================================\n")
