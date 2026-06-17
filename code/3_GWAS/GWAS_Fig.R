# =============================================================================
# Author: Haodong Tian
# Description: Generates GWAS visualizations — heritability vs. independent
#              SNP count scatter plot, and Manhattan plot for the most heritable
#              liver MRI feature (firstorder_TotalEnergy).
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

library(data.table)
library(qqman)


# The 59 retained features — features passing heritability and correlation QC
retained_names <- readRDS("/path/to/project/phenotype data/features with other images/retained_names.rds")

# Univariate LDSC heritability results, subsetted to retained features
LDSC <- read.table('/path/to/server/project/GWAS_regenie/LDSC_new/heritability_summary.tsv',
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
LDSC <- LDSC[LDSC$Feature %in% retained_names, ]
LDSC_res <- LDSC
LDSC_res$Z <- LDSC_res$h2 / LDSC_res$h2_se  # univariate LDSC Z score


# Count independent lead SNPs per feature at p < 1e-6 (from UVMR GWAS clumped results)
tsv_dir <- "/path/to/project/new_UVMR/retained_UVMR_summary_1e6_new"

snp_counts <- sapply(retained_names, function(feat) {
  fpath <- file.path(tsv_dir, paste0(feat, ".tsv"))
  if (!file.exists(fpath)) {
    return(0L)
  } else {
    dt <- fread(fpath)
    return(nrow(dt))
  }
})
snp_df <- data.frame(Feature   = retained_names,
                     snp_count = as.integer(snp_counts),
                     stringsAsFactors = FALSE)

snp_counts_5e8 <- sapply(retained_names, \(feat)
  if (file.exists(f <- file.path(tsv_dir, paste0(feat, ".tsv")))) sum(fread(f)$P < 5e-8) else 0L)
snp_df_5e8 <- data.frame(Feature = retained_names, snp_count_5e8 = snp_counts_5e8)

LDSC_final <- Reduce(function(x, y) merge(x, y, by = "Feature", all.x = TRUE),
                     list(LDSC_res, snp_df, snp_df_5e8))
# firstorder_TotalEnergy has the largest h2, # of SNPs (1e-6), and # of SNPs (5e-8)

plot(LDSC_final$snp_count_5e8, LDSC_final$h2)



### Manhattan plot for firstorder_TotalEnergy (the most heritable feature) ----

gwas <- fread(
  "/path/to/server/project/GWAS_regenie/MR_GWAS_new/firstorder_TotalEnergy.regenie",
  select = c("CHR", "BP", "P", "SNP")
)

gwas[, CHR := as.integer(CHR)]

# Retain all significant hits; subsample background points to limit plot size
sig <- gwas[P < 1e-4]
bg  <- gwas[P >= 1e-4][sample(.N, min(5e5, .N))]
gwas_plot <- rbind(sig, bg)
dim(gwas_plot)

manhattan(gwas_plot,
          suggestiveline = FALSE,
          genomewideline = -log10(5e-8),
          cex      = 0.4,
          cex.axis = 0.5,
          col      = c("#7CAAE2", "#183B5A"),
          ylim     = c(0, 80)
)
