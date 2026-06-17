# =============================================================================
# Author: Haodong Tian
# Description: Pre-processing for MR analysis — selects and aligns instrumental
#              variables from GWAS summary statistics.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================


### Pre-MR analysis steps:
# (1) Q-clustering of features based on local GWAS variants (for visualization)
#     Note: Q-clustering does not require the covariance matrix of genetic associations.
# (2) Individual-level rank-normal transformed (rint) and residualized phenotype data
#     (required for the covariance matrix of genetic associations in robust MVMR estimation)
# (3) LDSC heritability results (computed on HPC)


pathname <- '/path/to/project'


library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(MASS)          # for latent MR
library(LaplacesDemon) # for latent MR
library(caret)         # for quick correlation selection
library(MendelianRandomization)
library(susieR)
library(combinat)      # for MR-BMA
library(hash)          # for MR-BMA
library(pheatmap)
library(mr.raps)
library(data.table)


### Individual-level phenotype data (features must be rint-transformed, then residualized)
liver_data_all <- read.table(
  "/path/to/project/phenotype data/features with other images/rint_and_residualizing.txt",
  header = TRUE, sep = "", stringsAsFactors = FALSE)
dim(liver_data_all)  # 37791 x 76 = 17 covariates/ID columns + 59 retained features
hist(liver_data_all$shape_Flatness, n = 100)  # spot check distribution of one feature
apply(liver_data_all[, -(1:17)], 2, sd)
summary(apply(liver_data_all[, -(1:17)], 2, sd))


### --------------------------------------------------------------------------------------------------
### Q-clustering based on local GWAS variants -------------------------------------------------------
### --------------------------------------------------------------------------------------------------

##
## (1) LD clumping of the combined lead SNPs to obtain independent IVs
##

# Combined lead SNP coordinates for all 59 retained features
UVMRcorrds <- read.table(paste0(pathname, "/new_UVMR/UVMRcoords_new.txt"),
                         header = FALSE, sep = " ", stringsAsFactors = FALSE)
dim(UVMRcorrds)  # 461 x 3
rsID_list <- UVMRcorrds$V3

## LD clumping for the combined lead SNPs
dummy_sumstats <- data.table(SNP = rsID_list, P = 0.001)
dummy_file <- file.path(pathname, 'new_UVMR', 'PLINK_files', "dummy_sumstats.txt")
fwrite(dummy_sumstats, file = dummy_file, sep = "\t", quote = FALSE)

plink_path <- file.path(pathname, "PheWAS/plink_mac_20250615/plink")
bfile_path <- file.path(pathname, "PheWAS/1000G_QC")
out_prefix <- file.path(pathname, 'new_UVMR', 'PLINK_files', "rsID_clumped_results_for_Q")
plink_cmd  <- sprintf(
  '"%s" --bfile "%s" --clump "%s" --clump-p1 1 --clump-p2 1 --clump-r2 0.01 --clump-kb 1000 --out "%s"',
  plink_path, bfile_path, dummy_file, out_prefix)
system(plink_cmd, intern = TRUE)

clumped_results <- fread(paste0(out_prefix, ".clumped"))
final_rsID_list <- clumped_results$SNP  # independent SNPs after LD clumping
length(final_rsID_list)  # 282 independent SNPs each associated with >= 1 radiomics feature


##
## (2) Retrieve Bx and Bxse for all 59 features at the independent variants
##
folder <- paste0(pathname, "/new_UVMR/retained_UVMR_summary_1e6_new/")
feature_paths <- list.files(folder, pattern = "\\.tsv$", full.names = TRUE)
all_names <- c()
for (kk in 1:length(feature_paths)) {
  feature_path <- feature_paths[kk]
  feature_name <- sub("\\.tsv$", "", basename(feature_path))
  all_names    <- c(all_names, feature_name)
}
all_names  # the 59 retained feature names


# Pre-extracted local GWAS files are used to avoid reading full GWAS results
# files in R (which is substantially slower than on HPC).
gwas_dir <- "/path/to/server/project/GWAS_regenie/UVMRtoMVMR_new/local_MR_GWAS"

num_snps  <- length(final_rsID_list)
num_traits <- length(all_names)
num_snps; num_traits  # 282 independent SNPs; 59 retained features

Bx   <- matrix(NA_real_,      nrow = num_snps, ncol = num_traits,
               dimnames = list(final_rsID_list, all_names))
Bxse <- matrix(NA_real_,      nrow = num_snps, ncol = num_traits,
               dimnames = list(final_rsID_list, all_names))
EA   <- matrix(NA_character_,  nrow = num_snps, ncol = num_traits,
               dimnames = list(final_rsID_list, all_names))

for (trait in all_names) {
  file_path    <- file.path(gwas_dir, paste0(trait, ".regenie"))
  cat("Reading GWAS file:", file_path, "\n")
  gwas_data    <- fread(file_path, select = c("SNP", "ALLELE1", "BETA", "SE"))
  matched_data <- gwas_data[SNP %in% final_rsID_list]
  setkey(matched_data, SNP)
  matched_rsIDs <- intersect(final_rsID_list, matched_data$SNP)
  Bx[matched_rsIDs, trait]   <- matched_data[matched_rsIDs, BETA]
  Bxse[matched_rsIDs, trait] <- matched_data[matched_rsIDs, SE]
  EA[matched_rsIDs, trait]   <- matched_data[matched_rsIDs, ALLELE1]
  cat("Trait", trait, "done. Matched SNPs:", length(matched_rsIDs), "\n")
}
dim(Bx); dim(Bxse); dim(EA)  # all 282 x 59
View(Bx); View(Bxse); View(EA)


##
## (3) Q-clustering
##

# Exclude BMI-associated variants before clustering to reduce confounding
BMI_for_UVMR <- read.table(paste0(pathname, "/new_UVMR/", 'BMI', '.tsv'),
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
BMI_data_matched <- BMI_for_UVMR[match(final_rsID_list, BMI_for_UVMR$rsID), ]
rsID_vec <- as.vector(na.omit(
  BMI_data_matched$rsID[BMI_data_matched$pvalue > 5 * 10^(-4)]))  # key threshold
length(rsID_vec)  # 226

clustering_Bx_matrix   <- Bx[match(rsID_vec, rownames(Bx)), ]
clustering_Bxse_matrix <- Bxse[match(rsID_vec, rownames(Bxse)), ]
dim(clustering_Bx_matrix); dim(clustering_Bxse_matrix)  # 226 x 59


# Q statistic function: pairwise heterogeneity between two features using 2nd-order IVW
# Minimised numerically to give a symmetric dissimilarity measure for clustering.
get_Q <- function(k1, k2) {
  bl1   <- clustering_Bx_matrix[, k1];  bl1se <- clustering_Bxse_matrix[, k1]
  bl2   <- clustering_Bx_matrix[, k2];  bl2se <- clustering_Bxse_matrix[, k2]
  # Naive Q (ignoring sample correlation between features) is preferred because
  # its value is stable for highly correlated traits and the phenotypic correlation
  # matrix is not straightforward to obtain in all settings.
  the_Q_function <- function(b) {
    sum((bl2 - b * bl1)^2 / (bl2se^2 + b^2 * bl1se^2))
  }
  start_b_value <- mr_ivw(mr_input(bx = bl1, bxse = bl1se,
                                    by = bl2, byse = bl2se), weights = 'delta')@Estimate
  mr_plot(mr_input(bx = bl1, bxse = bl1se, by = bl2, byse = bl2se))
  b_lower <- start_b_value - 5 * mr_ivw(mr_input(bx = bl1, bxse = bl1se,
                                                    by = bl2, byse = bl2se), weights = 'delta')@StdError
  b_upper <- start_b_value + 5 * mr_ivw(mr_input(bx = bl1, bxse = bl1se,
                                                    by = bl2, byse = bl2se), weights = 'delta')@StdError
  numeric_res <- optim(par = start_b_value, fn = the_Q_function, method = "Brent",
                       lower = b_lower, upper = b_upper)$value
  return(numeric_res)
}

# Compute pairwise Q matrix (symmetric)
Q_matrix <- matrix(rep(0, ncol(clustering_Bx_matrix)^2),
                   ncol(clustering_Bx_matrix), ncol(clustering_Bx_matrix))
for (i in 1:(ncol(clustering_Bx_matrix) - 1)) {
  cat(paste0(i, '-'))
  for (j in (i + 1):ncol(clustering_Bx_matrix)) {
    Qval <- get_Q(k1 = i, k2 = j)
    Q_matrix[i, j] <- Qval; Q_matrix[j, i] <- Qval  # symmetric fill
  }
}
dim(Q_matrix)  # 59 x 59
colnames(Q_matrix) <- colnames(clustering_Bx_matrix)
rownames(Q_matrix) <- colnames(clustering_Bx_matrix)

Q_matrix <- readRDS(paste0(pathname, "/new_UVMR/Q_matrix.rds"))

# Diagnostics: check for negative Q values (should be 0)
n_neg <- sum(Q_matrix < 0, na.rm = TRUE)
message("Number of negative entries: ", n_neg)

Q_matrix_ <- Q_matrix

# Assess which transformation is more Gaussian for use as a dissimilarity
hist(Q_matrix_)
hist(log(1 + Q_matrix_ / (dim(clustering_Bx_matrix)[1] - 1)))
# Log-transformed Q (standardised by degrees of freedom) used for clustering
Q_matrix_used <- log(1 + Q_matrix_ / (dim(clustering_Bx_matrix)[1] - 1))


# Hierarchical clustering using average linkage on log-transformed Q dissimilarity
dist_mat <- as.dist(Q_matrix_used)
hc_cols  <- hclust(dist_mat, method = "average")
k        <- 10  # k=10 chosen by elbow method (see below); captures all prioritised features
clusters <- cutree(hc_cols, k = k)
table(clusters)
cluster_df <- data.frame(feature = names(clusters), cluster = as.integer(clusters),
                         row.names = NULL, stringsAsFactors = FALSE)
cluster_df <- cluster_df %>% arrange(cluster, feature)

pheatmap(Q_matrix_used,
         clustering_distance_rows = dist_mat, clustering_distance_cols = dist_mat,
         clustering_method        = "average",
         cutree_rows = k, cutree_cols = k,
         display_numbers = FALSE,
         number_format   = "%.1f",
         show_colnames   = FALSE,
         annotation_legend_pos = "left",
         treeheight_row = 0, treeheight_col = 50)

# Reference line: log(1 + chi2_0.05 / (n-1)) — Q values above this indicate
# significant heterogeneity between the two features (used for visual reference in heatmap)
log(1 + qchisq(1 - 0.05, nrow(clustering_Bx_matrix) - 1) / (dim(clustering_Bx_matrix)[1] - 1))


### Elbow method to select optimal k
d_mat <- as.matrix(dist_mat)
sse <- sapply(1:30, function(k) {
  cl    <- cutree(hc_cols, k = k)
  sse_k <- 0
  for (g in unique(cl)) {
    idx   <- which(cl == g)
    if (length(idx) <= 1) next
    d_sub <- d_mat[idx, idx]
    sse_k <- sse_k + sum(d_sub[upper.tri(d_sub)]^2)
  }
  sse_k
})
plot(1:30, sse, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters (k)", ylab = "Sum of Squared Distances",
     main = "Elbow Method for Hierarchical Clustering")
# k=10 is the selected value; provides good cluster separation while
# capturing all UVMR-highlighted features

View(cluster_df)
