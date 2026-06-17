# =============================================================================
# Author: Haodong Tian
# Description: Main UVMR and MVMR Mendelian Randomization PheWAS pipeline for
#              59 retained liver MRI radiomics features across multiple
#              cardiometabolic outcomes.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================


### UVMR + MVMR MR PheWAS analysis for the 59 retained liver features

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


OUTnames <- list.files(file.path(pathname, "new_UVMR"), pattern = "\\.tsv$", full.names = FALSE)
OUTnames <- OUTnames[OUTnames != "BMI.tsv"]  # BMI is used as a covariate, not an outcome
OUTnames <- sub("\\.tsv$", "", OUTnames)
OUTnames  # 12 outcomes of interest


## Read liver feature Q-clustering results (only for the 59 retained features)
cluster_df <- readRDS(paste0(pathname, "/new_UVMR/cluster_df.rds"))

## The 59 retained features — only these are included in MR analyses
retained_names <- readRDS("/path/to/project/phenotype data/features with other images/retained_names.rds")

## Read feature univariate LDSC heritability results
LDSC <- read.table('/path/to/server/project/GWAS_regenie/LDSC_new/heritability_summary.tsv',
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
LDSC <- LDSC[LDSC$Feature %in% retained_names, ]
LDSC_res <- LDSC
LDSC_res$Z <- LDSC_res$h2 / LDSC_res$h2_se  # univariate LDSC Z score


## Read individual-level phenotype data
# Required for estimating the covariance matrix of genetic associations with
# multiple exposures in MVMR-SuSiE (2nd-order multivariable IVW).
# Individual-level data is needed to account for the uncertainty (distribution)
# of the genetic association estimates with the exposures.
liver_data_all <- read.table(  # features must be rank-normal transformed then residualized
  "/path/to/project/phenotype data/features with other images/rint_and_residualizing.txt",
  header = TRUE, sep = "", stringsAsFactors = FALSE)
dim(liver_data_all)  # 37791 x 76 = ID + 16 covariates + 59 features


## Get file paths and names for all 59 retained liver MRI features
folder <- paste0(pathname, "/new_UVMR/retained_UVMR_summary_1e6_new/")
feature_paths <- list.files(folder, pattern = "\\.tsv$", full.names = TRUE)
all_names <- c()  # all liver feature names
for (kk in 1:length(feature_paths)) {
  feature_path  <- feature_paths[kk]
  feature_name  <- sub("\\.tsv$", "", basename(feature_path))
  all_names     <- c(all_names, feature_name)
}
all_names  # the 59 retained features


do_Steiger <- TRUE  # use Steiger filtering to remove SNPs more likely to first affect the outcome

for (OUTname in OUTnames) {

  OUT_UVMRmatrix   <- read.table(paste0(pathname, "/new_UVMR/", OUTname, ".tsv"),
                                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  OUT_UVMRsummary  <- OUT_UVMRmatrix
  print('---------------------------------------------')
  print(paste0('Current phenotype: ', OUTname,
               ' (', nrow(OUT_UVMRsummary), ' matched SNPs in this trait GWAS)'))

  fpath <- paste0(pathname, "/new_UVMR/UVMRpipeline_res/UVMR_highlighted_", OUTname, ".rds")
  if (file.exists(fpath)) {
    print('Results already exist for this phenotype; skipping')
    next
  }

  UVMR_RES <- data.frame(
    feature           = rep(NA, length(feature_paths)),
    cluster           = rep(NA, length(feature_paths)),
    UVMR_ID           = rep(NA, length(feature_paths)),
    h2                = rep(NA, length(feature_paths)),
    MVMR_SNP_number   = rep(NA, length(feature_paths)),
    UVMR_SNP_number   = rep(NA, length(feature_paths)),
    MVMR_BMI_est      = rep(NA, length(feature_paths)),
    MVMR_BMI_p        = rep(NA, length(feature_paths)),
    IVW_est           = rep(NA, length(feature_paths)),
    IVW_pvalue        = rep(NA, length(feature_paths)),
    Q_pvalue          = rep(NA, length(feature_paths)),
    Fvalue            = rep(NA, length(feature_paths)),
    weighted_median   = rep(NA, length(feature_paths)),
    weighted_mode     = rep(NA, length(feature_paths)),
    MR_RAPS           = rep(NA, length(feature_paths)),
    P_adj_bonf        = rep(NA, length(feature_paths))
  )

  cat('Feature-wide UVMR: ')
  # Run UVMR separately for each of the 59 retained features under the current outcome
  for (kk in 1:length(feature_paths)) {
    cat(paste0(kk, '-'))

    feature_path <- feature_paths[kk]
    feature_name <- sub("\\.tsv$", "", basename(feature_path))

    # Load feature-specific IV summary data (p-value threshold 1e-6 for IV selection)
    dat <- read.table(feature_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # Record feature metadata
    UVMR_RES$feature[kk] <- feature_name
    UVMR_RES$cluster[kk] <- cluster_df$cluster[match(feature_name, cluster_df$feature)]
    UVMR_RES$UVMR_ID[kk] <- kk
    UVMR_RES$h2[kk]      <- LDSC_res$h2[match(feature_name, LDSC_res$Feature)]

    # Remove IVs absent from the outcome GWAS summary data
    if (sum(is.na(match(dat$SNP, OUT_UVMRsummary$rsID))) != 0) {
      dat <- dat[!is.na(match(dat$SNP, OUT_UVMRsummary$rsID)), ]
    }

    # Remove IVs absent from the BMI GWAS summary data
    BMI_for_UVMR <- read.table(paste0(pathname, "/new_UVMR/", 'BMI', '.tsv'),
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (sum(is.na(match(dat$SNP, BMI_for_UVMR$rsID))) != 0) {
      dat <- dat[!is.na(match(dat$SNP, BMI_for_UVMR$rsID)), ]
    }


    # Strand check (run only for the first feature since all MRI exposures share the same IV format)
    if (kk == 1) {
      Data1 <- dat
      Data2 <- OUT_UVMRsummary[match(dat$SNP, OUT_UVMRsummary$rsID), ]
      Data3 <- BMI_for_UVMR[match(dat$SNP, BMI_for_UVMR$rsID), ]
      # In the outcome GWAS, proxy SNPs have non-empty proxyID; exclude them for strand checking
      # since strand assignment requires true EA/nonEA from the original SNP
      cat('\n   ==== same-strand check between [MRI exposure dat] and [outcome data]   \n')
      is_empty <- function(x) is.na(x) | x == ""
      Data2 <- Data2[is_empty(Data2$proxyID), ]
      twoData_strand_check(Data1, Data2)  # function defined in twoData_strand_check.R
      cat('\n   ==== same-strand check between [MRI exposure dat] and [BMI data]   \n')
      Data3 <- Data3[is_empty(Data3$proxyID), ]
      twoData_strand_check(Data1, Data3)
    }
    # Result: all SNPs are on the same strand (same-strand rate = 100%); no ambiguous palindromic
    # SNPs (0.45 < EAF < 0.55); MR harmonization via effect allele alone is therefore safe.


    # Steiger filtering:
    # Steiger filtering removes SNPs that explain more variance in the outcome than in the exposure,
    # which are likely invalid IVs acting via a reverse pathway.
    # Note: because radiomics features may themselves be surrogate exposures (i.e. G -> liver
    # condition -> radiomics feature), the Steiger model may not strictly hold for all SNPs.
    # However, applying Steiger filtering still improves MR reliability by removing obvious
    # reverse-causal SNPs, while robust MR methods (MR-Median, MR-Mode, MR-RAPS) provide
    # additional safeguards.
    if (do_Steiger & (nrow(dat) != 0)) {
      # Outcome effective sample size (use median across SNPs for meta-analysed GWAS such as CAD)
      clean_GWAS_data_path <- list.files(
        "/path/to/project/LDSC/clean_GWAS",
        pattern = paste0("^", OUTname, "\\.tsv\\.gz$"), full.names = TRUE)
      N2 <- median(fread(clean_GWAS_data_path)$effective_sample_size, na.rm = TRUE)

      # Compute squared correlation with the outcome for each IV
      local_OUT_UVMRsummary <- OUT_UVMRsummary[match(dat$SNP, OUT_UVMRsummary$rsID), ]
      ref_effect_allele <- dat$ALLELE1
      by_vector  <- local_OUT_UVMRsummary$beta *
        (-1 + 2 * (local_OUT_UVMRsummary$trueEA == ref_effect_allele))
      byse_vector <- local_OUT_UVMRsummary$se
      t_scores   <- by_vector / byse_vector
      abs_r2     <- sqrt(t_scores^2 / (t_scores^2 + N2 - 2))

      # Compute squared correlation with the exposure for each IV
      N1       <- nrow(liver_data_all)  # 37791 — radiomics feature GWAS sample size
      t_scores <- dat$BETA / dat$SE
      abs_r1   <- sqrt(t_scores^2 / (t_scores^2 + N1 - 2))

      # Fisher's z transformation for the one-sided Steiger test
      z1      <- 1/2 * log((1 + abs_r1) / (1 - abs_r1))
      z2      <- 1/2 * log((1 + abs_r2) / (1 - abs_r2))
      final_z <- (z1 - z2) / sqrt(1 / (N1 - 3) + 1 / (N2 - 3))

      # One-sided test with Bonferroni correction
      # H0: abs_r1 <= abs_r2 (SNP explains more variance in outcome than exposure)
      # Retain SNPs that do NOT significantly reject H0 (i.e. those consistent with
      # the exposure-first direction); this conservative approach preserves power while
      # robust MR methods handle any remaining invalid IVs.
      one_side_pvalues <- pnorm(final_z, 0, 1)
      p_threhold <- 0.05 / nrow(dat)
      dat <- dat[one_side_pvalues >= p_threhold, ]
    }


    # Proceed only if IVs remain after Steiger filtering
    if (nrow(dat) != 0) {

      if (nrow(dat) >= 2) {  # MVMR requires SNPs >= number of exposures
        ###### MVMR with BMI adjustment -----------------------------------------------
        UVMR_RES$MVMR_SNP_number[kk] <- nrow(dat)
        # SNP effect harmonization
        local_OUT_UVMRsummary <- OUT_UVMRsummary[match(dat$SNP, OUT_UVMRsummary$rsID), ]
        ref_effect_allele <- dat$ALLELE1
        by_vector   <- local_OUT_UVMRsummary$beta *
          (-1 + 2 * (local_OUT_UVMRsummary$trueEA == ref_effect_allele))
        byse_vector <- local_OUT_UVMRsummary$se

        BMI_data_matched <- BMI_for_UVMR[match(dat$SNP, BMI_for_UVMR$rsID), ]
        BMI_harmonized   <- BMI_data_matched$beta *
          (-1 + 2 * (BMI_data_matched$trueEA == dat$ALLELE1))

        MVMRres <- mr_mvivw(mr_mvinput(
          bx    = cbind(dat$BETA, BMI_harmonized),
          bxse  = cbind(dat$SE, BMI_data_matched$se),
          by    = by_vector,
          byse  = byse_vector))

        UVMR_RES$MVMR_BMI_est[kk] <- MVMRres@Estimate[1]  # estimate for liver feature (index 1)
        UVMR_RES$MVMR_BMI_p[kk]   <- MVMRres@Pvalue[1]
      }

      ####### UVMR with BMI-variant removal -------------------------------------------
      # Remove BMI-associated variants (p < 5e-4) before running UVMR, to reduce
      # horizontal pleiotropy via the adiposity pathway
      BMI_pvalues <- (BMI_for_UVMR$pvalue)[match(dat$SNP, BMI_for_UVMR$rsID)]
      dat <- dat[BMI_pvalues > 5 * 10^(-4), ]  # key threshold for BMI-associated variant removal

      UVMR_RES$UVMR_SNP_number[kk] <- nrow(dat)

      # Outcome harmonization
      local_OUT_UVMRsummary <- OUT_UVMRsummary[match(dat$SNP, OUT_UVMRsummary$rsID), ]
      ref_effect_allele <- dat$ALLELE1
      by_vector   <- local_OUT_UVMRsummary$beta *
        (-1 + 2 * (local_OUT_UVMRsummary$trueEA == ref_effect_allele))
      byse_vector <- local_OUT_UVMRsummary$se

      if (nrow(dat) != 0) {
        # Primary UVMR: IVW (random-effects)
        MRres <- mr_ivw(mr_input(bx = dat$BETA, bxse = dat$SE,
                                 by = by_vector, byse = byse_vector))
        mr_plot(mr_input(bx = dat$BETA, bxse = dat$SE, by = by_vector, byse = byse_vector),
                orientate = TRUE)
        dat$snp      <- paste0('snp_', 1:nrow(dat))
        dat$proxysnp <- local_OUT_UVMRsummary$proxyID[match(dat$SNP, local_OUT_UVMRsummary$rsID)]

        UVMR_RES$IVW_est[kk]    <- MRres@Estimate
        UVMR_RES$IVW_pvalue[kk] <- MRres@Pvalue
        UVMR_RES$Q_pvalue[kk]   <- MRres@Heter.Stat[2]  # Cochran Q heterogeneity p-value
        UVMR_RES$Fvalue[kk]     <- MRres@Fstat           # F-statistic for weak instrument check
      }

      ### Robust UVMR sensitivity analyses: MR-Median, MR-Mode, MR-RAPS ---
      if (nrow(dat) >= 3) {
        MRmedian <- mr_median(mr_input(bx = dat$BETA, bxse = dat$SE,
                                       by = by_vector, byse = byse_vector),
                              weighting = 'simple')
        UVMR_RES$weighted_median[kk] <- MRmedian@Pvalue

        MRmode <- mr_mbe(mr_input(bx = dat$BETA, bxse = dat$SE,
                                  by = by_vector, byse = byse_vector))
        UVMR_RES$weighted_mode[kk] <- MRmode@Pvalue

        MRRAPS <- mr.raps(data.frame(beta.exposure = dat$BETA, beta.outcome = by_vector,
                                     se.exposure = dat$SE, se.outcome = byse_vector),
                          diagnostics = FALSE, over.dispersion = TRUE)
        UVMR_RES$MR_RAPS[kk] <- 2 * pnorm(-abs(MRRAPS$beta.hat / MRRAPS$beta.se))
      }
    }

  }

  UVMR_RES$P_adj_bonf <- round(p.adjust(as.numeric(UVMR_RES$IVW_pvalue), method = "bonferroni"), 3)
  UVMR_RES$P_adj_BH   <- round(p.adjust(as.numeric(UVMR_RES$IVW_pvalue), method = "BH"), 3)

  ### Weak instrument check: all F > 10 is required
  if (sum((UVMR_RES$Fvalue < 10)[!is.na(UVMR_RES$Fvalue)]) == 0) {
    print('All F > 10; good')
  } else {
    print('Caution: some features have F < 10 — potential weak instrument bias')
  }

  saveRDS(UVMR_RES, file = paste0(pathname, "/new_UVMR/UVMRpipeline_res/UVMR_results_", OUTname, ".rds"))


  ### Highlighted feature selection ---------------------------------------------------

  ## Strict filter (UVMR_superhighlighted):
  # (1) UVMR BH-FDR p-value < 0.05 (primary criterion; BMI-associated SNPs already removed)
  # (2) All robust MR p-values < 0.05 (MR-Median, MR-Mode, MR-RAPS)
  # (3) MVMR-BMI p-value < 0.05 (BMI-adjusted MVMR confirms the association)
  # (4) MVMR estimate direction consistent with UVMR estimate
  UVMR_superhighlighted <- na.omit(UVMR_RES[
    UVMR_RES$P_adj_BH < 0.05 &
      pmax(UVMR_RES$weighted_median, UVMR_RES$weighted_mode, UVMR_RES$MR_RAPS) < 0.05 &
      UVMR_RES$MVMR_BMI_p < 0.05 &
      (UVMR_RES$MVMR_BMI_est * UVMR_RES$IVW_est) > 0, ])


  ## Lenient filter (UVMR_highlighted):
  # Requires criterion (1) above, plus at least 2 out of 3 additional criteria:
  # (2) All robust MR p-values < 0.05
  # (3) MVMR-BMI p-value < 0.05
  # (4) Consistent direction between MVMR and UVMR estimates
  # At least 3 SNPs are implicitly required for robust MR methods to be computed.
  cond3    <- pmax(UVMR_RES$weighted_median, UVMR_RES$weighted_mode, UVMR_RES$MR_RAPS) < 0.05
  cond4    <- UVMR_RES$MVMR_BMI_p < 0.05
  cond5    <- (UVMR_RES$MVMR_BMI_est * UVMR_RES$IVW_est) > 0
  cond_sum <- cond3 + cond4 + cond5
  UVMR_highlighted <- na.omit(UVMR_RES[
    UVMR_RES$P_adj_BH < 0.05 &   # criterion (1) mandatory
      cond_sum >= 2,              # at least 2 of the 3 additional criteria satisfied
    ])

  if (nrow(UVMR_highlighted) == 0) {
    print('No features pass UVMR (lenient version); skipping this phenotype')
    next
  }

  print(paste0("Number of superhighlighted features: ", nrow(UVMR_superhighlighted)))
  print(paste0("Number of highlighted features: ", nrow(UVMR_highlighted)))
  saveRDS(UVMR_superhighlighted,
          file = paste0(pathname, "/new_UVMR/UVMRpipeline_res/UVMR_superhighlighted_", OUTname, ".rds"))
  saveRDS(UVMR_highlighted,
          file = paste0(pathname, "/new_UVMR/UVMRpipeline_res/UVMR_highlighted_", OUTname, ".rds"))


  ## Cluster enrichment analysis: are the highlighted features enriched in any Q-cluster?
  all_features    <- all_names
  cluster_labels  <- cluster_df$cluster[match(all_features, cluster_df$feature)]
  names(cluster_labels) <- all_features
  selected_features <- UVMR_highlighted$feature

  clusters <- unique(cluster_labels)
  results  <- data.frame(cluster = clusters, OR = NA, p = NA, selected_count = NA)
  for (i in seq_along(clusters)) {
    cl       <- clusters[i]
    in_cl    <- names(cluster_labels)[cluster_labels == cl]
    not_in_cl <- setdiff(names(cluster_labels), in_cl)
    a <- sum(in_cl %in% selected_features)
    b <- length(in_cl) - a
    c <- sum(not_in_cl %in% selected_features)
    d <- length(not_in_cl) - c
    # One-sided Fisher's exact test for over-representation
    fisher_res <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    results$OR[i]             <- fisher_res$estimate
    results$p[i]              <- fisher_res$p.value
    results$selected_count[i] <- a
  }
  saveRDS(results,
          file = paste0(pathname, "/new_UVMR/UVMRpipeline_res/enrichment_", OUTname, ".rds"))


  if (nrow(UVMR_highlighted) == 1) {
    print('Only one feature highlighted; MVMR variable selection requires >= 2 features; skipping')
    next
  }

  ### MVMR variable selection via MVMR-SuSiE =============================================

  ## Only UVMR-highlighted features are included

  ### STEP 0: Collect all IVs and run PLINK clumping to obtain independent SNPs ----------
  UVMR_highlighted_names <- UVMR_highlighted$feature
  UVMR_highlighted_kk    <- UVMR_highlighted$UVMR_ID
  length(UVMR_highlighted_kk)

  rsID_list <- c()
  for (kk in UVMR_highlighted_kk) {
    feature_path <- feature_paths[kk]
    dat <- read.table(feature_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    rsID_list <- c(rsID_list, dat$SNP)
  }
  rsID_list <- unique(rsID_list)

  # Create a dummy GWAS summary file for PLINK clumping (uniform p-value avoids PLINK issues)
  dummy_sumstats <- data.table(SNP = rsID_list, P = 0.001)
  dummy_file <- file.path(pathname, 'new_UVMR', 'PLINK_files', 'forMVMRSuSiE_dummy_sumstats.txt')
  fwrite(dummy_sumstats, file = dummy_file, sep = "\t", quote = FALSE)

  plink_path <- file.path(pathname, "PheWAS/plink_mac_20250615/plink")
  bfile_path <- file.path(pathname, "PheWAS/1000G_QC")
  out_prefix <- file.path(pathname, 'new_UVMR', 'PLINK_files', "rsID_clumped_results")
  plink_cmd  <- sprintf(
    '"%s" --bfile "%s" --clump "%s" --clump-p1 1 --clump-p2 1 --clump-r2 0.01 --clump-kb 1000 --out "%s"',
    plink_path, bfile_path, dummy_file, out_prefix)
  system(plink_cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)

  clumped_results  <- fread(paste0(out_prefix, ".clumped"))
  final_rsID_list  <- clumped_results$SNP

  print(paste0("MVMR | total unique SNPs: ", length(rsID_list),
               " | after-clumping independent SNPs: ", length(final_rsID_list)))


  ### STEP 1: Retrieve summary data for all post-clumping independent SNPs ---------------
  # Pre-extracted local GWAS files are used because reading full GWAS files from the server
  # in R is substantially slower than on HPC; local files are pre-subset on the HPC.
  gwas_dir <- "/path/to/server/project/GWAS_regenie/UVMRtoMVMR_new/local_MR_GWAS"

  num_snps  <- length(final_rsID_list)
  num_traits <- length(UVMR_highlighted_names)
  Bx   <- matrix(NA_real_,      nrow = num_snps, ncol = num_traits,
                 dimnames = list(final_rsID_list, UVMR_highlighted_names))
  Bxse <- matrix(NA_real_,      nrow = num_snps, ncol = num_traits,
                 dimnames = list(final_rsID_list, UVMR_highlighted_names))
  EA   <- matrix(NA_character_,  nrow = num_snps, ncol = num_traits,
                 dimnames = list(final_rsID_list, UVMR_highlighted_names))

  cat("Reading GWAS files from server and matching SNPs: ")
  for (trait in UVMR_highlighted_names) {
    cat(paste0(which(UVMR_highlighted_names == trait), '-'))
    file_path    <- file.path(gwas_dir, paste0(trait, ".regenie"))
    gwas_data    <- fread(file_path, select = c("SNP", "ALLELE1", "BETA", "SE"))
    matched_data <- gwas_data[SNP %in% final_rsID_list]
    setkey(matched_data, SNP)
    matched_rsIDs <- intersect(final_rsID_list, matched_data$SNP)
    Bx[matched_rsIDs, trait]   <- matched_data[matched_rsIDs, BETA]
    Bxse[matched_rsIDs, trait] <- matched_data[matched_rsIDs, SE]
    EA[matched_rsIDs, trait]   <- matched_data[matched_rsIDs, ALLELE1]
  }

  print(paste0('MVMR summary matrix dimensions (Bx): ', nrow(Bx), ' x ', ncol(Bx)))


  ### STEP 2: BMI-variant removal -------------------------------------------------------
  BMI_data_matched <- BMI_for_UVMR[match(final_rsID_list, BMI_for_UVMR$rsID), ]
  rsID_vec <- as.vector(na.omit(
    BMI_data_matched$rsID[BMI_data_matched$pvalue > 5 * 10^(-4)]))  # key threshold

  # Ensure all remaining SNPs also have outcome GWAS data
  rsID_vec <- rsID_vec[rsID_vec %in% OUT_UVMRsummary$rsID]

  print(paste0('After BMI-variant removal, SNP count: ', length(rsID_vec)))


  ### STEP 3: MR-cML-BIC to exclude outlier/pleiotropic SNPs ---------------------------
  # MR-cML-BIC is recommended by the MVMR-SuSiE paper for outlier removal prior to MVMR.
  # It is applied per-exposure to flag SNPs with evidence of horizontal pleiotropy.
  cat(paste0('MR-cML-BIC outlier removal for MVMR: '))
  undesired_rsID <- c()
  for (kkk in 1:ncol(Bx)) {

    # Outcome data harmonization
    local_OUT_UVMRsummary <- OUT_UVMRsummary[match(rsID_vec, OUT_UVMRsummary$rsID), ]
    ref_effect_allele <- EA[match(rsID_vec, rownames(EA)), 1]
    by_vector   <- local_OUT_UVMRsummary$beta *
      (-1 + 2 * (local_OUT_UVMRsummary$trueEA == ref_effect_allele))
    byse_vector <- local_OUT_UVMRsummary$se

    # MR-cML-BIC for this exposure
    bx_vector   <- Bx[match(rsID_vec, rownames(Bx)), kkk]
    bxse_vector <- Bxse[match(rsID_vec, rownames(Bxse)), kkk]
    mr_plot(mr_input(bx = bx_vector, bxse = bxse_vector,
                     by = by_vector, byse = byse_vector), orientate = TRUE)

    MR_cML_res <- mr_cML(mr_input(bx = bx_vector, bxse = bxse_vector,
                                  by = by_vector, byse = byse_vector),
                         MA = FALSE, DP = FALSE, n = 41743)
    cat(paste0(kkk, '(', length(MR_cML_res@BIC_invalid), ')-'))
    undesired_rsID <- c(undesired_rsID, rsID_vec[MR_cML_res@BIC_invalid])

    # Scatter plot coloring valid vs. invalid SNPs
    df <- data.frame(bx = bx_vector, bxse = bxse_vector,
                     by = by_vector, byse = byse_vector)
    df$group <- "Valid"
    df$group[MR_cML_res@BIC_invalid] <- "Invalid"
    color_map <- c("Valid" = "black", "Invalid" = "red")
    df_plot <- transform(df, bx = abs(bx), by = ifelse(bx < 0, -by, by))
    ggplot(df_plot, aes(x = bx, y = by, color = group)) +
      geom_point() +
      geom_errorbar(aes(ymin = by - 1.96 * byse, ymax = by + 1.96 * byse), width = 0) +
      geom_errorbarh(aes(xmin = bx - 1.96 * bxse, xmax = bx + 1.96 * bxse), height = 0) +
      scale_color_manual(values = color_map) +
      xlab("SNP effect on exposure (bx)") + ylab("SNP effect on outcome (by)") +
      theme_minimal() + theme(legend.title = element_blank())
  }
  undesired_rsID <- unique(undesired_rsID)


  ### STEP 4: MVMR-SuSiE fitting -------------------------------------------------------
  MMVR_rsID_vec <- rsID_vec[!rsID_vec %in% undesired_rsID]

  print(paste0('After removing ', length(undesired_rsID),
               ' outlier SNPs (MR-cML-BIC), final SNP count: ', length(MMVR_rsID_vec)))

  ## Prepare outcome data
  local_OUT_UVMRsummary <- OUT_UVMRsummary[match(MMVR_rsID_vec, OUT_UVMRsummary$rsID), ]
  ref_effect_allele <- EA[match(MMVR_rsID_vec, rownames(EA)), 1]
  by_vector   <- local_OUT_UVMRsummary$beta *
    (-1 + 2 * (local_OUT_UVMRsummary$trueEA == ref_effect_allele))
  byse_vector <- local_OUT_UVMRsummary$se

  ## Prepare exposure summary matrices
  Bx_matrix_used   <- Bx[match(MMVR_rsID_vec, rownames(Bx)), ]
  Bxse_matrix_used <- Bxse[match(MMVR_rsID_vec, rownames(Bxse)), ]
  dim(Bx_matrix_used); dim(Bxse_matrix_used)

  ## Covariance matrix of the G-X association estimates (required for 2nd-order MVMR-IVW)
  # The phenotypic correlation matrix from individual-level data is used to approximate
  # the sampling covariance between SNP-exposure effect estimates across correlated exposures.
  Corr <- cor(liver_data_all[, match(colnames(Bx), colnames(liver_data_all)), drop = FALSE])
  Sigma_exposures_list <- list()
  for (jj in 1:length(MMVR_rsID_vec)) {
    Sigma_exposures <- diag(Bxse_matrix_used[jj, ]) %*% Corr %*% diag(Bxse_matrix_used[jj, ])
    Sigma_exposures_list[[jj]] <- Sigma_exposures
  }

  ## Iterative MVMR-SuSiE: update the 2nd-order standard error iteratively
  # The advanced (2nd-order) SE accounts for uncertainty in the IV-exposure estimator,
  # equivalent to a profile likelihood correction. Iteration continues until convergence.
  beta0      <- rep(0, ncol(Bx_matrix_used))
  beta0_prev <- rep(Inf, length(beta0))
  tol        <- 1e-6; max_iter <- 50
  iter       <- 0
  while (max(abs(beta0 - beta0_prev)) > tol && iter < max_iter) {
    iter <- iter + 1
    beta0_prev <- beta0

    # Update 2nd-order error term (profile-likelihood equivalent)
    byse_vector_adv <- sapply(1:length(MMVR_rsID_vec), function(jj) {
      sqrt(t(beta0) %*% Sigma_exposures_list[[jj]] %*% beta0 + byse_vector[jj]^2)
    })

    # Weighted regression to standardise error terms
    Sigma_inv_0.5     <- diag(1 / byse_vector_adv)
    weighted_responses  <- Sigma_inv_0.5 %*% by_vector
    weighted_Covariates <- Sigma_inv_0.5 %*% as.matrix(Bx_matrix_used)

    # SuSiE for sparse variable selection in MVMR
    res <- susie(weighted_Covariates, weighted_responses, L = 10, intercept = FALSE,
                 standardize = FALSE, max_iter = 500,
                 estimate_residual_variance = FALSE, residual_variance = 1)

    cs_list   <- susie_get_cs(res)
    PIPmatrix <- res$alpha[1:length(cs_list$cs), , drop = FALSE]
    PMImatrix <- res$mu[1:length(cs_list$cs), , drop = FALSE]
    beta0     <- apply(PIPmatrix * PMImatrix, 2, sum)  # posterior mean -> updated causal effect

    cat("Iteration:", iter, "max diff:", max(abs(beta0 - beta0_prev)), "\n")
  }


  ### STEP 5: Store MVMR-SuSiE results -------------------------------------------------
  cs_len <- ifelse(is.null(cs_list$cs) || length(cs_list$cs) == 0, 0, length(cs_list$cs))
  if (cs_len > 0) {
    PIP_vals <- res$alpha[1:cs_len, , drop = FALSE]
  } else {
    PIP_vals <- rep(NA, ncol(res$alpha))
  }

  print(paste0('MVMR-SuSiE done | number of credible sets: ', cs_len))
  exposure_names <- colnames(Bx_matrix_used)
  MVMR_SuSiE_res <- data.frame(
    Feature = exposure_names,
    PIP     = t(PIP_vals),
    Cluster = cluster_df$cluster[match(exposure_names, cluster_df$feature)])

  saveRDS(MVMR_SuSiE_res,
          file = paste0(pathname, "/new_UVMR/UVMRpipeline_res/MVMRsusie_", OUTname, ".rds"))
}

# Note: warnings() are mostly from MR-RAPS (expected behaviour)




### ----------------------------------------------------------------------
### Visualization --------------------------------------------------------

OUTnames  # phenotype / outcome names
all_names # feature / exposure names


files <- file.path(pathname, "new_UVMR/UVMRpipeline_res",
                   paste0("UVMR_highlighted_", OUTnames, ".rds"))

mr_pwas_df <- do.call(rbind, lapply(seq_along(OUTnames), function(i) {
  f <- files[i]; if (!file.exists(f)) return(NULL)
  x <- readRDS(f); if (is.null(x) || nrow(x) == 0) return(NULL)
  x <- subset(x, feature %in% all_names)
  if (nrow(x) == 0) return(NULL)
  data.frame(
    outcome = OUTnames[i], feature = x$feature,
    sign    = ifelse(x$IVW_est > 0, "pos", ifelse(x$IVW_est < 0, "neg", NA_character_)),
    stringsAsFactors = FALSE)
}))
if (is.null(mr_pwas_df))
  mr_pwas_df <- data.frame(outcome = character(), feature = character(), sign = character())

## Order rows by cluster membership (same cluster grouped together)
cf <- merge(data.frame(feature = all_names, idx = seq_along(all_names)),
            cluster_df, by = "feature", all.x = TRUE)
cf <- cf[order(cf$cluster, cf$idx), ]
feature_levels <- cf$feature

## Complete grid (all OUTnames x all features shown even if NA)
grid    <- expand.grid(outcome = OUTnames, feature = feature_levels,
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
plot_df <- merge(grid, mr_pwas_df, by = c("outcome", "feature"), all.x = TRUE, sort = FALSE)

## Cluster bounding boxes for the heatmap
pos <- setNames(seq_along(feature_levels), feature_levels)
cb  <- do.call(rbind, lapply(split(cf$feature, cf$cluster, drop = TRUE), function(v) {
  r <- range(pos[v])
  data.frame(xmin = 0.5, xmax = length(OUTnames) + 0.5,
             ymin = r[1] - 0.5, ymax = r[2] + 0.5)
}))

## Cluster label positions (midpoint of each cluster group)
rng    <- tapply(pos[cf$feature], cf$cluster, range)
lab_df <- do.call(rbind, lapply(names(rng), function(k) {
  r <- rng[[k]]
  data.frame(cluster = k, y_mid_idx = mean(r), stringsAsFactors = FALSE)
}))
lab_df$y_lab <- factor(feature_levels[round(lab_df$y_mid_idx)], levels = feature_levels)
lab_df$x_lab <- length(OUTnames) + 0.6
lab_df$label <- paste0("C", lab_df$cluster)


## Read MVMR-SuSiE credible set membership for overlay on heatmap
mvmr_files <- file.path(pathname, "new_UVMR", "UVMRpipeline_res",
                        paste0("MVMRsusie_", OUTnames, ".rds"))
mvmr_marks <- do.call(rbind, lapply(seq_along(OUTnames), function(i) {
  f <- mvmr_files[i]
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f)
  if (is.null(x) || nrow(x) == 0) return(NULL)
  pip_cols <- grep("^PIP(\\.\\d+)?$", names(x), value = TRUE)
  if (length(pip_cols) == 0) return(NULL)
  pip_cols <- if ("PIP" %in% pip_cols) {
    c("PIP", setdiff(pip_cols, "PIP"))
  } else {
    pip_cols[order(as.integer(sub("PIP\\.(\\d+)", "\\1", pip_cols)))]
  }
  do.call(rbind, lapply(seq_along(pip_cols), function(k) {
    pk <- x[, c("Feature", pip_cols[k])]
    colnames(pk) <- c("feature", "PIP")
    pk <- pk[pk$feature %in% feature_levels & is.finite(pk$PIP) & pk$PIP > 0, , drop = FALSE]
    if (nrow(pk) == 0) return(NULL)
    pk <- pk[order(-pk$PIP), ]
    k_cut <- which(cumsum(pk$PIP) >= 0.5)[1]
    if (is.na(k_cut)) k_cut <- nrow(pk)
    data.frame(outcome = OUTnames[i],
               feature = pk$feature[seq_len(k_cut)],
               cs      = paste0("CS", k),
               stringsAsFactors = FALSE)
  }))
}))


# Manually reordered trait names for display
OUTnames_reorder <- c('TG', 'HDL_C', 'TG_HDL_C', 'Lp_a', 'LDL_C', 'CAD',
                      'HbA1c', 'T2D', 'liver_fat', 'liver_iron', 'Cirrhosis', 'liver_cancer')
length(OUTnames_reorder) == length(OUTnames)


# Identify outcomes with a single highlighted feature (for special annotation)
dt_tmp       <- as.data.table(plot_df)
single_marks <- dt_tmp[!is.na(sign), .SD[.N == 1], by = outcome][, .(outcome, feature)]

# Grid line data for heatmap
grid_h <- data.frame(y  = seq(0.5, length(feature_levels) + 0.5, by = 1),
                     x1 = 0.5, x2 = length(OUTnames_reorder) + 0.5)
grid_v <- data.frame(x  = seq(0.5, length(OUTnames_reorder) + 0.5, by = 1),
                     y1 = 0.5, y2 = length(feature_levels) + 0.5)

# Rename outcomes for display (liver iron -> liver cT1, TG_HDL_C -> TG:HDL_C)
map <- c(liver_iron = "liver_cT1", TG_HDL_C = "TG:HDL_C")
ren <- function(x) { x <- as.character(x); i <- x %in% names(map); x[i] <- map[x[i]]; x }
plot_df$outcome     <- ren(plot_df$outcome)
OUTnames_reorder    <- ren(OUTnames_reorder)
OUTnames            <- ren(OUTnames)
if (exists("mvmr_marks"))   mvmr_marks$outcome   <- ren(mvmr_marks$outcome)
if (exists("single_marks")) single_marks$outcome <- ren(single_marks$outcome)

outcome_labels <- c(
  "TG"       = "TG",       "HDL_C"    = "HDL-C",    "TG:HDL_C" = "TG:HDL-C",
  "Lp_a"     = "Lp(a)",    "LDL_C"    = "LDL-C",    "CAD"       = "CAD",
  "HbA1c"    = "HbA1c",    "T2D"      = "T2D",
  "liver_fat"= "Liver fat","liver_iron"="Liver iron",
  "Cirrhosis"= "Cirrhosis","liver_cancer"="Liver cancer","liver_cT1"="Liver cT1"
)


## Heatmap: star (CS1) and triangle (CS2) symbols mark MVMR-SuSiE credible set membership
p <- ggplot(plot_df,
            aes(x = factor(outcome, levels = OUTnames_reorder),
                y = factor(feature, levels = feature_levels),
                fill = sign)) +
  geom_tile(width = 0.9, height = 0.9) +
  geom_segment(data = grid_h,
               aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, color = "grey85", linewidth = 0.2) +
  geom_segment(data = grid_v,
               aes(x = x, xend = x, y = y1, yend = y2),
               inherit.aes = FALSE, color = "grey85", linewidth = 0.2) +
  {if (!is.null(cb) && nrow(cb) > 0)
    geom_rect(data = cb, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA, color = "grey50", linewidth = 0.4)} +
  {if (!is.null(lab_df))
    geom_text(data = lab_df, aes(x = x_lab, y = y_lab, label = label),
              inherit.aes = FALSE, hjust = 0, size = 3.3, color = "grey20")} +
  {if (!is.null(mvmr_marks) && nrow(mvmr_marks) > 0)
    geom_point(data = mvmr_marks,
               aes(x = factor(outcome, levels = OUTnames),
                   y = factor(feature, levels = feature_levels),
                   shape = cs),
               inherit.aes = FALSE, size = 1.8, color = "black")} +
  {if (nrow(single_marks) > 0)
    geom_point(data = single_marks,
               aes(x = factor(outcome, levels = OUTnames_reorder),
                   y = factor(feature, levels = feature_levels)),
               inherit.aes = FALSE, shape = 8, size = 1.8, color = "black")} +
  scale_fill_manual(values = c(pos = "pink", neg = "lightblue"),
                    na.value = NA, guide = "none") +
  scale_shape_manual(values = c(CS1 = 8, CS2 = 17), guide = "none") +  # 8=star, 17=triangle
  scale_x_discrete(drop = FALSE, labels = outcome_labels) +
  scale_y_discrete(drop = FALSE) +
  coord_cartesian(xlim = c(0.5, length(OUTnames) + 1.8), clip = "off") +
  theme_minimal(base_size = 11) +
  theme(panel.grid     = element_blank(),
        axis.title     = element_blank(),
        axis.text.x    = element_text(angle = 45, hjust = 1),
        axis.text.y    = element_text(size = 8),
        plot.margin    = margin(5, 30, 5, 5))

p


### Generate summary MR results table (highlighted features only) ----------------------
OUTnames <- list.files(file.path(pathname, "new_UVMR"), pattern = "\\.tsv$", full.names = FALSE)
OUTnames <- OUTnames[OUTnames != "BMI.tsv"]
OUTnames <- sub("\\.tsv$", "", OUTnames)

uvmr_files <- file.path(pathname, "new_UVMR/UVMRpipeline_res",
                        paste0("UVMR_highlighted_", OUTnames, ".rds"))
mvmr_files <- file.path(pathname, "new_UVMR", "UVMRpipeline_res",
                        paste0("MVMRsusie_", OUTnames, ".rds"))

MRres_table <- do.call(rbind, lapply(seq_along(OUTnames), function(i) {
  UVMRres <- readRDS(uvmr_files[i])

  MVMRres <- readRDS(mvmr_files[i])
  MVMRres$Cluster <- NULL
  pip_cols <- grep("^PIP", names(MVMRres), value = TRUE)
  fmt      <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.3f", as.numeric(x)))
  MVMRres$PIPres <- if (length(pip_cols) == 1) {
    fmt(MVMRres[[pip_cols]])
  } else {
    apply(MVMRres[, pip_cols, drop = FALSE], 1,
          function(r) paste(fmt(r), collapse = ";"))
  }

  UVMRres$PIPres  <- MVMRres$PIPres[match(UVMRres$feature, MVMRres$Feature)]
  UVMRres$Outcome <- OUTnames[i]
  UVMRres         <- UVMRres[, c("Outcome", setdiff(names(UVMRres), "Outcome"))]
  UVMRres
}))

MRres_table$Outcome[MRres_table$Outcome == 'TG_HDL_C']  <- 'TG:HDL_C'
MRres_table$Outcome[MRres_table$Outcome == 'liver_iron'] <- 'liver_cT1'

saveRDS(MRres_table, file = "/path/to/project/new_UVMR/MRres_table.rds")
