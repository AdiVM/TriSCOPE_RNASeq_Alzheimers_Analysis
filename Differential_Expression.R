library(arrow)
library(lme4)
library(dplyr)
library(glmmTMB)

# --- Load Data ---
expr_path <- "/Users/adithyamadduri/Desktop/Projects/Patel_Lab/CellMatrix_with_genenames.parquet"
meta_path <- "/Users/adithyamadduri/Desktop/Projects/Patel_Lab/New_CellMetadataSyn1848517_with_rnaBatch.parquet"

expr <- read_parquet(expr_path)     # Genes x cells
meta <- read_parquet(meta_path)

cell_types <- c("Mic", "Ast", "In", "Opc", "Oli", "Ex")

meta$apoe_genotype <- as.character(meta$apoe_genotype)
meta$apoe_genotype <- gsub("\\.0", "", meta$apoe_genotype)

# Mapping: 33 -> e3/e3, 34 -> e3/e4, etc.
meta$apoe_genotype <- recode(meta$apoe_genotype,
                             "23" = "e2/e3",
                             "33" = "e3/e3",
                             "34" = "e3/e4",
                             "44" = "e4/e4",
                             "24" = "e2/e4",
                             "22" = "e2/e2")

# Convert to factor and set reference level to "e3/e3"
meta$apoe_genotype <- factor(meta$apoe_genotype)
meta$apoe_genotype <- relevel(meta$apoe_genotype, ref = "e3/e3")

for (celltype in cell_types) {
  
  # --- Subset Mic Cells ---
  meta_mic <- meta %>% filter(broad.cell.type == celltype)
  rownames(meta_mic) <- meta_mic$TAG
  
  expr_mat <- as.data.frame(expr)
  rownames(expr_mat) <- expr_mat$index
  expr_mat$index <- NULL
  #meta_mic <- meta_mic[as.character(meta_mic$projid) != '21189544',]
  expr_mic <- expr_mat[, meta_mic$TAG, drop=FALSE]
  table(meta_mic$projid)
  
  
  cat("Initial gene count:", nrow(expr_mic), "\n")
  
  lib <- colSums(expr_mic)
  
  # --- 1. Remove genes with zero variance ---
  gene_var <- apply(expr_mic, 1, var)
  expr_mic_novar <- expr_mic[gene_var > 0, ]
  cat("After removing zero-variance genes:", nrow(expr_mic_novar), "\n")
  
  # --- 2. Filter on CPM > 1 in >=3 cells ---
  lib_size_millions <- colSums(expr_mic_novar) / 1e6
  expr_cpm <- sweep(expr_mic_novar, 2, lib_size_millions, FUN="/")
  gene_filter <- rowSums(expr_cpm > 1) >= 10
  expr_mic_filt <- expr_mic_novar[gene_filter, ]
  cat("After CPM filter:", nrow(expr_mic_filt), "\n")
  
  # --- 3. Remove genes expressed only in AD or only in controls ---
  meta_mic$projid <- as.factor(meta_mic$projid)
  meta_mic$pmi <- as.numeric(meta_mic$pmi)
  
  # Handling age at death
  # Replace "90+" (or any entry containing "90+") with 90
  meta_mic$age_death <- as.character(meta_mic$age_death)
  meta_mic$age_death[meta_mic$age_death == "90+"] <- "90"
  
  # Convert to numeric
  meta_mic$age_death <- as.numeric(meta_mic$age_death)
  
  # Round to nearest integer
  meta_mic$age_death <- round(meta_mic$age_death, 0)
  
  meta_mic$msex <- as.factor(meta_mic$msex)
  meta_mic$AD <- as.factor(ifelse(meta_mic$ceradsc %in% c(1, 2), 1, 0))
  print(table(meta_mic$AD, meta_mic$projid) %>% rowSums() %>% table())
  # meta_mic$AD <- as.factor(ifelse(!is.na(meta_mic$age_first_ad_dx), 1, 0))
  
  expr_mic_final <- expr_mic_filt
  cat("Skipping AD/control expression filter; using all remaining genes:", nrow(expr_mic_final), "\n")
  
  # --- 4/5. Loop through all genes for DE analysis ---
  results <- list()
  libsize <- colSums(expr_mic_final)
  
  counter <- 0
  for (gene in rownames(expr_mic_final)) {
    counter <- counter + 1
    if(counter%%1000 == 0){
      print(counter)
    }
    y <- as.numeric(expr_mic_final[gene, ])
    dt <- meta_mic
    dt$y <- y
    dt$libsize <- libsize[dt$TAG]
    
    # Only keep cells where all covariates are non-missing
    to_keep <- !is.na(dt$pmi) & !is.na(dt$age_death) & !is.na(dt$msex) & !is.na(dt$AD) & !is.na(y) & !is.na(dt$apoe_genotype) & !is.na(dt$sequencingBatch)
    dt2 <- dt[to_keep, ]
    
    dt2$apoe_genotype <- factor(dt2$apoe_genotype)
    dt2$apoe_genotype <- relevel(dt2$apoe_genotype, ref = "e3/e3")
    y2 <- y[to_keep]
    
    #dt2$projid <- relevel(dt2$projid, ref = '10101291')
    
    lib2 <- lib[to_keep]
    
    
    
    fit_simple <- tryCatch(
      glmmTMB(y~AD + offset(log(lib2)) + pmi + age_death + msex + apoe_genotype + sequencingBatch,data=dt2,ziformula=~1,family=poisson),
      #glmer(y ~ AD + offset(log(lib2)) + (1|projid), data=dt2, family=poisson, nAGQ=10),
      error = function(e) NULL
    )
    
    
    
    #if (!is.null(fit_simple) && !isSingular(fit_simple)) {
    # Store result
    #}
    
    if (!is.null(fit_simple)) {
      summ <- summary(fit_simple)
      
      res <- data.frame(
        gene = gene,
        estimate_AD = summ$coefficients$cond["AD1", "Estimate"],
        se_AD = summ$coefficients$cond["AD1", "Std. Error"],
        pval_AD = summ$coefficients$cond["AD1", "Pr(>|z|)"],
        n_cells = length(y2)
      )
      results[[gene]] <- res
    }
  }
  
  summary(fit_simple)
  
  final_results <- do.call(rbind, results)
  final_results$p_adj <- p.adjust(final_results$pval_AD, method = "fdr")
  final_results$log2FC <- final_results$estimate_AD / log(2)
  final_results$DEG <- ifelse((final_results$p_adj < 0.05) & (abs(final_results$estimate_AD / log(2)) > 0.25), "True", "False")
  
  out_path <- sprintf("/Users/adithyamadduri/Desktop/Projects/Patel_Lab/Cerad_batch_DE_Outputs/poisson_DE_results_%s.csv", celltype)
  write.csv(final_results, out_path, row.names = FALSE)
  cat(sprintf("DE results written to %s\n", out_path))
  
}
