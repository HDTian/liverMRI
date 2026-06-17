# =============================================================================
# Author: Haodong Tian
# Description: Multivariable MR using factor analysis (MVMR-FA/SuSiE) to
#              estimate direct causal effects of liver MRI latent factors on
#              outcomes.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================


### MVMR-FA

# Based on the UVMR-highlighted MRI features:
# In addition to MVMR-SuSiE, MVMR-FA is conducted as a sensitivity analysis.
# MVMR-FA identifies groups of MRI features via factor analysis applied to local
# genetic information (i.e. in the MR context), and then performs MVMR using the
# resulting factors. This approach accounts for the high correlation among
# radiomics features while still estimating direct causal effects.


source('/path/to/project/new_UVMR/MVMR_FA_functions.R')


### Data collection ==============================================================================

pathname <- '/path/to/project'
OUTnames  # all trait names
OUTname <- OUTnames[12]
OUTname <- 'T2D'


### MVMR-FA analysis =================================================================================
MVMR_FA_RES <- list()
for (OUTname in OUTnames) {

  cat(paste0('================================= \n'))
  cat(paste0('Current outcome: ', OUTname, ' \n'))

  UVMR_highlighted <- readRDS(
    paste0(pathname, "/new_UVMR/UVMRpipeline_res/UVMR_highlighted_", OUTname, ".rds"))


  ## Restrict analysis to UVMR-highlighted features

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

  clumped_results <- fread(paste0(out_prefix, ".clumped"))
  final_rsID_list <- clumped_results$SNP

  print(paste0("MVMR | total unique SNPs: ", length(rsID_list),
               " | after-clumping independent SNPs: ", length(final_rsID_list)))


  ### STEP 1: Retrieve summary data for all post-clumping independent SNPs ---------------
  # Pre-extracted local GWAS files are used (see MRpipeline.R for rationale)
  gwas_dir <- "/path/to/server/project/GWAS_regenie/UVMRtoMVMR_new/local_MR_GWAS"

  num_snps   <- length(final_rsID_list)
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

  rsID_vec <- rsID_vec[rsID_vec %in% OUT_UVMRsummary$rsID]

  print(paste0('After BMI-variant removal, SNP count: ', length(rsID_vec)))


  ### STEP 3: MR-cML-BIC to exclude outlier/pleiotropic SNPs ---------------------------
  # See MRpipeline.R for the rationale; same procedure applied here for MVMR-FA.
  cat(paste0('MR-cML-BIC outlier removal for MVMR: '))
  undesired_rsID <- c()
  for (kkk in 1:ncol(Bx)) {

    local_OUT_UVMRsummary <- OUT_UVMRsummary[match(rsID_vec, OUT_UVMRsummary$rsID), ]
    ref_effect_allele <- EA[match(rsID_vec, rownames(EA)), 1]
    by_vector   <- local_OUT_UVMRsummary$beta *
      (-1 + 2 * (local_OUT_UVMRsummary$trueEA == ref_effect_allele))
    byse_vector <- local_OUT_UVMRsummary$se

    bx_vector   <- Bx[match(rsID_vec, rownames(Bx)), kkk]
    bxse_vector <- Bxse[match(rsID_vec, rownames(Bxse)), kkk]
    mr_plot(mr_input(bx = bx_vector, bxse = bxse_vector,
                     by = by_vector, byse = byse_vector), orientate = TRUE)

    MR_cML_res <- mr_cML(mr_input(bx = bx_vector, bxse = bxse_vector,
                                  by = by_vector, byse = byse_vector),
                         MA = FALSE, DP = FALSE, n = 41743)
    cat(paste0(kkk, '(', length(MR_cML_res@BIC_invalid), ')-'))
    undesired_rsID <- c(undesired_rsID, rsID_vec[MR_cML_res@BIC_invalid])

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
  undesired_rsID  <- unique(undesired_rsID)
  MMVR_rsID_vec   <- rsID_vec[!rsID_vec %in% undesired_rsID]

  print(paste0('After removing ', length(undesired_rsID),
               ' outlier SNPs (MR-cML-BIC), final SNP count: ', length(MMVR_rsID_vec)))


  ### STEP 4: Prepare final summary data for MVMR-FA -----------------------------------

  ## Outcome data and harmonization
  local_OUT_UVMRsummary <- OUT_UVMRsummary[match(MMVR_rsID_vec, OUT_UVMRsummary$rsID), ]
  ref_effect_allele <- EA[match(MMVR_rsID_vec, rownames(EA)), 1]
  by_vector   <- local_OUT_UVMRsummary$beta *
    (-1 + 2 * (local_OUT_UVMRsummary$trueEA == ref_effect_allele))
  byse_vector <- local_OUT_UVMRsummary$se

  ## Exposure summary matrices
  Bx_matrix_used   <- Bx[match(MMVR_rsID_vec, rownames(Bx)), ]
  Bxse_matrix_used <- Bxse[match(MMVR_rsID_vec, rownames(Bxse)), ]
  dim(Bx_matrix_used); dim(Bxse_matrix_used)


  ### STEP 5: MVMR-FA fitting (Bayesian Factor Analysis) --------------------------------
  # MVMR_BFA performs sparse Bayesian factor analysis on the exposure matrix, extracts
  # latent factor scores, and runs multivariable IVW using those scores as the exposures.
  # Sparsity on the factor loadings is enforced via BayesFM priors (Kmax = 2 factors here).

  MVMR_FA_res <- tryCatch({
    MVMR_BFA_res <- MVMR_BFA(
      exposure_beta    = Bx_matrix_used,
      outcome_beta     = by_vector,
      outcome_se       = byse_vector,
      Kmax             = 2, model = "fixed", Nid = 2,
      burnin           = 2000, iter = 10000,
      nu0 = 7, kappa = 0.15, kappa0 = 3, xi0 = 0.8,
      seed             = 725,
      plot_diagnostics = TRUE,
      generate_report  = TRUE
    )
    list(
      loadings           = MVMR_BFA_res$bfa_summary$alpha,
      loading_plot       = MVMR_BFA_res$plots$loadings,
      factor_correlation = MVMR_BFA_res$diagnostics$factor_correlations,
      MVMRres            = MVMR_BFA_res$mvmr_results$model
    )
  }, error = function(e) {
    message("MVMR_BFA failed: ", conditionMessage(e))
    NULL
  })

  MVMR_FA_RES[[OUTname]] <- MVMR_FA_res

}
# Analysis complete

saveRDS(MVMR_FA_RES, file = paste0(pathname, "/new_UVMR/MVMR_FA/MVMR_FA_RES.rds"))


### Summary of results by outcome
# A good MVMR-FA result: significant effect of at least one latent factor on the outcome,
# with factor loadings whose direction matches the UVMR heatmap (after adjusting for
# MVMR effect sign).

# Lipid outcomes
# LDL-C, Lp(a), CAD: MVMR-BFA failed (insufficient support for 2 latent factors)
MVMR_FA_RES$TG       # (near-)null MVMR result; loading pattern inconsistent
MVMR_FA_RES$HDL_C    # significant MVMR result; loading pattern inconsistent
MVMR_FA_RES$TG_HDL_C # null MVMR result

# Glycaemic outcomes
MVMR_FA_RES$T2D    # significant MVMR result; loading pattern consistent with UVMR
MVMR_FA_RES$HbA1c  # significant MVMR result; loading pattern consistent with UVMR

### Liver-related outcomes
MVMR_FA_RES$liver_fat    # significant MVMR result; loading pattern consistent with UVMR
MVMR_FA_RES$liver_iron   # significant MVMR result; loading pattern consistent with UVMR
MVMR_FA_RES$Cirrhosis    # significant MVMR result; loading pattern consistent with UVMR
MVMR_FA_RES$liver_cancer # significant MVMR result; loading pattern consistent with UVMR


### Conclusion:
# For T2D/HbA1c and liver-related outcomes, MVMR-FA results are consistent with
# both the UVMR PheWAS and MVMR-SuSiE analyses (significant effect + consistent direction).
# MVMR-FA additionally characterises the factor structure, indicating which sets of
# correlated MRI features share the same genetic pathway to the nominated outcome.


### Export MVMR-FA results table (loadings, factor correlations, MVMR estimates) ------
library(data.table)

obj2mat <- function(x) { m <- as.matrix(x); storage.mode(m) <- "character"; m }
pad     <- function(m, k) { if (ncol(m) < k) cbind(m, matrix("", nrow(m), k - ncol(m))) else m }

MR_FA_table <- do.call(rbind, lapply(names(MVMR_FA_RES), function(out) {
  obj  <- MVMR_FA_RES[[out]]$MVMRres
  tab1 <- obj2mat(rbind(
    colnames(cbind(Est = obj@Estimate, SE = obj@StdError,
                   CI_lower = obj@CILower, CI_upper = obj@CIUpper)),
    cbind(Est = obj@Estimate, SE = obj@StdError,
          CI_lower = obj@CILower, CI_upper = obj@CIUpper)))
  tab2 <- obj2mat(MVMR_FA_RES[[out]]$factor_correlation)
  tab3 <- obj2mat(rbind(colnames(MVMR_FA_RES[[out]]$loadings), MVMR_FA_RES[[out]]$loadings))
  k    <- max(ncol(tab1), ncol(tab2), ncol(tab3))
  rbind(cbind(out,                 matrix("", 1, k - 1)),
        cbind('MVMR results',      matrix("", 1, k - 1)),
        pad(tab1, k),
        cbind('Factor correlation',matrix("", 1, k - 1)),
        pad(tab2, k),
        cbind('Loading info',      matrix("", 1, k - 1)),
        pad(tab3, k),
        matrix("", 1, k))
}))
rownames(MR_FA_table) <- NULL

saveRDS(MR_FA_table, file = "/path/to/project/new_UVMR/MR_FA_table.rds")
