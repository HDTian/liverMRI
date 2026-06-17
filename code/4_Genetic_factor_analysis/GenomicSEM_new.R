# =============================================================================
# Author: Haodong Tian
# Description: Genomic SEM factor analysis on liver MRI radiomics features —
#   runs LDSC to estimate the genetic covariance matrix (S) and its sampling
#   covariance (V), then performs EFA (even chromosomes) and CFA (odd
#   chromosomes) to identify latent genetic factors, with final fitting using
#   all chromosomes.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

### Genomic SEM — overview
# Starting from 59 retained liver MRI radiomics features (high r_g among all
# 200 features makes numeric inversion unstable), we focus on category-defined
# lead features: for each imaging category we pick the 4 most heritable
# features, yielding ~28 liver MRI lead features.
# Pipeline:
#   Step 0: select category-representative lead features
#   Step 1: munge GWAS summary statistics (done on Broad Server)
#   Step 2: LDSC — estimate genetic covariance (S) and sampling covariance (V)
#            separately for EVEN and ODD chromosomes, then for ALL chromosomes
#   Step 3: Genomic SEM factor analysis
#            EFA (even-chr only): factanal() to determine the number of factors
#                                 and their loading patterns
#            CFA (odd-chr only):  usermodel() to test model fit (CFI, SRMR,
#                                 AIC) and estimate parameters
#            Final model:         refit chosen factor structure using ALL chr



pathname <- '/path/to/project'

library(readxl)
library(GenomicSEM) # packageVersion("GenomicSEM")  # 0.0.5
library(data.table)
library(ggplot2)
library(Matrix)  # contains nearPD() for smoothing near-non-positive-definite matrices
library(stats)   # contains factanal() for naive factor analysis (ignores SE of the input covariance matrix)
library(semPlot) # path diagram visualisation
library(lavaan)




### Step 0: get the category lead liver MRI feature names ---------------------------------------
### ---------------------------------------------------------------------------------------------


retained_names <- readRDS("/path/to/project/phenotype data/features with other images/retained_names.rds")

LDSC <- read.table('/path/to/server/project/GWAS_regenie/LDSC_new/heritability_summary.tsv', header=TRUE, sep="\t", stringsAsFactors=FALSE)
LDSC <- LDSC[ LDSC$Feature %in% retained_names,  ]
dim(LDSC)  # 59 6 = Feature h2 se(h2) lambda_GC LDSC_intercept se(LDSC_intercept)

DT <- as.data.table(LDSC)
DT[, category := sub("_.*", "", Feature)]   # extract category prefix (everything before the first "_")
name_vec <- DT[order(-h2), head(.SD, 4), by = category]$Feature
name_vec

# Check the genetic correlation structure; avoid pairs that are too large or near zero
df <- fread("/path/to/server/project/GWAS_regenie/LDSC_new/genetic_correlation.tsv", sep = "\t")
dim(df) #19900     6
df_retained <- df[  (df$Feature1 %in% name_vec) & (df$Feature2 %in% name_vec)  , ]
dim(df_retained) # 91  6
View(df_retained)

cate_features <- name_vec  # these are used for LDSC analysis; further exclusion based on h^2 will be applied later




### Step 1: munged files ------------------------------------------------------------------------
### ---------------------------------------------------------------------------------------------
# Munged files can be taken directly from Broad Server results.




### Step 2: LDSC fitting ------------------------------------------------------------------------
### ---------------------------------------------------------------------------------------------
# Run multivariable LD-Score regression to obtain:
#   S — the genetic covariance matrix
#   V — the corresponding sampling covariance matrix (via Jackknife)
#
# Must use the ldsc() function from GenomicSEM because off-diagonal elements of
# the V matrix (sampling covariances between entries of S) cannot be derived
# from univariate LDSC runs alone.


### Non-BMI-adjusted GWAS — munged files folder
liver_features_munged_folder <- '/path/to/server/project/GWAS_regenie/LDSC_new/munged_files'
# Focus only on the ~14 category-defined lead features.


# Liver feature files
liver_feature_files <- list.files(liver_features_munged_folder, pattern = "\\.sumstats\\.gz$", full.names = TRUE, recursive = FALSE)
liver_feature_files  # vector of munged file paths
liver_feature_names <- sub("\\.sumstats\\.gz$", "", basename(liver_feature_files)); liver_feature_names  # corresponding feature name vector

# Keep only the category-representative lead features
cate_features  # the most heritable lead liver MRI feature within each category
keep_idx <- liver_feature_names %in% cate_features
liver_feature_files <- liver_feature_files[keep_idx]; liver_feature_names <- liver_feature_names[keep_idx]
liver_feature_files   # munged file paths for category-lead liver MRI features
liver_feature_names   # corresponding feature names

# Additional traits included in the LDSC run:
# [liver features] liver_fat, liver_iron


liver_feature_files_plus_liver_fat_iron <- c(liver_feature_files,
                                             '/path/to/project/LDSC/munged_GWAS/liver_fat.sumstats.gz',
                                             '/path/to/project/LDSC/munged_GWAS/liver_iron.sumstats.gz')
liver_feature_names_plus_liver_fat_iron <- c(liver_feature_names, 'liver_fat', 'liver_iron')



# LDSC fitting with ODD and EVEN chromosomes to obtain S_odd/V_odd and S_even/V_even.
# Splitting by chromosome parity prevents overfitting: EFA structure learned on
# EVEN chromosomes is tested independently on ODD chromosomes in the CFA step.

# EVEN chromosome -------------
setwd( file.path(pathname, "GenomicSEM/LDSC_cate_liverMRI_EVEN") )  # LDSC log output directory
LDSC_path <- file.path(pathname, 'LDSC', 'eur_w_ld_chr')
LDSCoutput_EVEN <- ldsc(
  traits          = liver_feature_files_plus_liver_fat_iron,
  trait.names     = liver_feature_names_plus_liver_fat_iron,
  sample.prev     = rep(NA, length(liver_feature_files_plus_liver_fat_iron)),   # NA for continuous traits
  population.prev = rep(NA, length(liver_feature_files_plus_liver_fat_iron)),   # only relevant for liability-scale h^2 of binary traits
  ld              = LDSC_path,
  wld             = LDSC_path,
  stand           = TRUE,   # also returns the genetic correlation matrix (easier to inspect r_g and for heatmap visualisation)
  select          = "EVEN"
)

dim( LDSCoutput_EVEN$S )  # 16*16 — valid entries: 16*(16-1)/2 + 16 = 136
dim( LDSCoutput_EVEN$V )  # 136*136

save(LDSCoutput_EVEN,
     file = "/path/to/project/GenomicSEM/LDSCoutput_cate_liverMRI_EVEN.RData")



# ODD chromosome ---------------
setwd( file.path(pathname, "GenomicSEM/LDSC_cate_liverMRI_ODD") )  # LDSC log output directory
LDSC_path <- file.path(pathname, 'LDSC', 'eur_w_ld_chr')
LDSCoutput_ODD <- ldsc(
  traits          = liver_feature_files_plus_liver_fat_iron,
  trait.names     = liver_feature_names_plus_liver_fat_iron,
  sample.prev     = rep(NA, length(liver_feature_files_plus_liver_fat_iron)),   # NA for continuous traits
  population.prev = rep(NA, length(liver_feature_files_plus_liver_fat_iron)),   # only relevant for liability-scale h^2 of binary traits
  ld              = LDSC_path,
  wld             = LDSC_path,
  stand           = TRUE,   # also returns the genetic correlation matrix
  select          = "ODD"
)

dim( LDSCoutput_ODD$S )  # 16*16 — valid entries: 16*(16-1)/2 + 16 = 136
dim( LDSCoutput_ODD$V )  # 136*136

save(LDSCoutput_ODD,
     file = "/path/to/project/GenomicSEM/LDSCoutput_cate_liverMRI_ODD.RData")




### Concordance check: are ODD- and EVEN-chromosome estimates consistent? ---------------------
# In a clean LDSC analysis we expect ODD and EVEN chromosome genetic correlation
# estimates to agree within sampling uncertainty.
load("/path/to/project/GenomicSEM/LDSCoutput_cate_liverMRI_EVEN.RData")
load("/path/to/project/GenomicSEM/LDSCoutput_cate_liverMRI_ODD.RData")

vech <- function(M) M[lower.tri(M, diag = TRUE)]
S_odd  <- vech(LDSCoutput_ODD$S);  S_even  <- vech(LDSCoutput_EVEN$S)
Var_odd <- diag(LDSCoutput_ODD$V); Var_even <- diag(LDSCoutput_EVEN$V)
Z_diff <- (S_odd - S_even) / sqrt(Var_odd + Var_even)
p_diff <- 2 * pnorm(-abs(Z_diff))  # length(p_diff)  # 465 = 30*29/2 + 30

# QQ-plot of p-values for ODD vs EVEN differences
p <- p_diff[is.finite(p_diff) & p_diff > 0 & p_diff <= 1]
n <- length(p)
obs <- -log10(sort(p))
exp <- -log10(ppoints(n))
plot(exp, obs, xlab="Expected -log10(p)", ylab="Observed -log10(p)", pch=16)
abline(0, 1, col="red")

library(ggplot2)
p <- sort(p[p>0 & p<=1]); p[p==0] <- .Machine$double.xmin; n <- length(p)
exp <- -log10((1:n - 0.5)/n); obs <- -log10(p); k <- 1:n
ci_l <- -log10(qbeta(0.975, k, n-k+1)); ci_u <- -log10(qbeta(0.025, k, n-k+1))
ggplot(data.frame(exp,obs,ci_l,ci_u), aes(exp,obs)) +
  geom_ribbon(aes(ymin=ci_l,ymax=ci_u), fill="grey80", alpha=1.0) +
  geom_point(size=.8) + geom_abline(slope=1,intercept=0,col="red") +
  labs(x="Expected -log10(p)", y="Observed -log10(p)") + theme_classic()

# QQ plots of Z statistics comparing odd- and even-chromosome estimates showed
# approximately linear behaviour without evidence of tail inflation, indicating
# that differences were consistent with sampling variability.




# ALL chromosomes: ODD + EVEN ---------------
setwd( file.path(pathname, "GenomicSEM/LDSC_cate_liverMRI_ALL") )  # LDSC log output directory
LDSC_path <- file.path(pathname, 'LDSC', 'eur_w_ld_chr')
LDSCoutput_ALL <- ldsc(
  traits          = liver_feature_files_plus_liver_fat_iron,
  trait.names     = liver_feature_names_plus_liver_fat_iron,
  sample.prev     = rep(NA, length(liver_feature_files_plus_liver_fat_iron)),   # NA for continuous traits
  population.prev = rep(NA, length(liver_feature_files_plus_liver_fat_iron)),   # only relevant for liability-scale h^2 of binary traits
  ld              = LDSC_path,
  wld             = LDSC_path,
  stand           = TRUE,
  select          = FALSE
)

dim( LDSCoutput_ALL$S )  # 30*30 — valid entries: 30*(30-1)/2 + 30 = 465
dim( LDSCoutput_ALL$V )  # 465*465

save(LDSCoutput_ALL,
     file = "/path/to/project/GenomicSEM/LDSCoutput_cate_liverMRI_ALL.RData")

load("/path/to/project/GenomicSEM/LDSCoutput_cate_liverMRI_ALL.RData")




### Step 3: SEM factor analysis with LDSC-derived S and V --------------------------------------------------------------------
### --------------------------------------------------------------------------------------------------------------------------


### Trait exclusion based on h^2 quality ----------------------
# Traits with low heritability (h^2 not significantly different from zero)
# degrade Genomic SEM fitting: nearPD smoothing introduces large distortions
# when the genetic signal is weak. We therefore exclude traits whose h^2
# Z-score falls below a threshold of 7.

liver_feature_names_plus_liver_fat_iron  # the 30 trait names

liver_h2 <- DT[order(-h2), head(.SD, 4), by = category]  # h^2 for the 28 category-lead liver MRI features
large_h2_features <- liver_h2$Feature[liver_h2$h2/liver_h2$h2_se > 7]  # features with sufficiently significant h^2

large_h2_features_plus_liver_fat_iron <- c(large_h2_features, 'liver_fat')
large_h2_features_plus_liver_fat_iron  # 10 traits (9 lead MRI features + liver_fat)
# Note: liver_iron is excluded because its h^2 is too low (h^2 Z-score < 7),
# which causes poor SEM fitting.
# For reference: liver_fat h^2 = 0.1642 (SE 0.0188), Z = 8.75;
#                liver_iron h^2 was not sufficiently significant.




### EFA -------------------------------------------------------
rownames(LDSCoutput_EVEN$S) <- colnames(LDSCoutput_EVEN$S)
keep <- intersect(large_h2_features_plus_liver_fat_iron, colnames(LDSCoutput_EVEN$S))
EVEN_S_sub <- LDSCoutput_EVEN$S[keep, keep]

Ssmooth <- as.matrix((nearPD(EVEN_S_sub, corr = FALSE))$mat)  # project to nearest positive-definite matrix
max(Ssmooth - EVEN_S_sub)  # 0.001003099 — record the smoothing difference introduced by nearPD


set.seed(1123)
EFA <- factanal(covmat = Ssmooth, factors = 3, rotation = "promax", start = NULL)  # naive factor analysis with a given number of latent factors

print(EFA$loadings, cutoff = 0.3)
# When a liver MRI feature has no loading above the threshold, assign it to the
# factor with the largest absolute loading to ensure every indicator belongs to
# exactly one factor.
L <- as.matrix(EFA$loadings)
L_keep <- t(apply(L, 1, function(x){
  a <- abs(x); keep <- which(a > 0.3)
  if (!length(keep)) keep <- which.max(a)
  y <- rep(NA_real_, length(x)); y[keep] <- x[keep]
  y
}))
colnames(L_keep) <- colnames(L); rownames(L_keep) <- rownames(L)
L_keep




### CFA -------------------------------------------------------
make_lavaan_syntax <- function(L_keep, scale = TRUE){
  stopifnot(is.matrix(L_keep))
  syntax <- c()
  for (f in colnames(L_keep)) {
    feats <- rownames(L_keep)[!is.na(L_keep[, f])]
    if (length(feats) == 0) next
    if (scale) { feats[1] <- paste0("NA*", feats[1]) }  # free the first loading for scale identification (recommended for GenomicSEM)
    line <- paste0(f, " =~ ", paste(feats, collapse = " + "))
    syntax <- c(syntax, line)
  }
  paste(syntax, collapse = "\n")
}

make_lavaan_syntax(L_keep)

genomicSEMresult <- usermodel(LDSCoutput_ODD,           # use odd-chromosome genetic covariance for independent CFA validation
                              estimation = "DWLS",      # Diagonally Weighted Least Squares
                              model = make_lavaan_syntax(L_keep),  # negative residual variances will be constrained: trait ~~ a*trait, a > 0.0001
                              CFIcalc = TRUE,
                              std.lv = TRUE,            # unit variance identification: constrain latent factor variances to 1
                              imp_cov = FALSE)
current_result <- genomicSEMresult$results
neg_trait_name <- current_result$lhs[ (current_result$lhs %in% rownames(L_keep)) & (current_result$Unstand_Est < 0) ]  # traits with negative residual variance

# Construct constrained model syntax to handle negative residual variances
labs        <- letters[seq_along(neg_trait_name)]
model_solve_neg <- paste0( make_lavaan_syntax(L_keep),
                           sprintf("\n%s ~~ %s*%s\n%s > .001", neg_trait_name, labs, neg_trait_name, labs) )
genomicSEMresult <- usermodel(LDSCoutput_ODD,
                              estimation = "DWLS",
                              model = model_solve_neg,
                              CFIcalc = TRUE,
                              std.lv = TRUE,
                              imp_cov = FALSE)




### Formal model selection iteration -----------------------------------------------------------------------------
### ---------------------------------------------------------------------------------------------------------------
# large_h2_features was filtered using h^2 Z-score > 7 (see trait exclusion section above).

# Final trait set for Genomic SEM fitting
large_h2_features_plus_liver_fat_iron <- c(large_h2_features, 'liver_fat')
# liver_iron is excluded: its h^2 is too low (see trait exclusion note above).

saveRDS(large_h2_features_plus_liver_fat_iron,
        "/path/to/project/GenomicSEM/final_GenomicSEM_names.rds")

# GenomicSEM_names <- readRDS("/path/to/project/GenomicSEM/final_GenomicSEM_names.rds")


# Reference heritability values used for quality filtering:
# liver_fat:  Total Observed Scale h2 = 0.1642 (SE 0.0188),  Z = 8.75
# liver_iron: Total Observed Scale h2 = 0.1282 (SE 0.0222),  Z = 5.77
# Genetic Correlation between liver_fat and liver_iron: 0.6752 (0.1369)



rownames(LDSCoutput_EVEN$S) <- colnames(LDSCoutput_EVEN$S)
keep <- intersect(large_h2_features_plus_liver_fat_iron, colnames(LDSCoutput_EVEN$S))
EVEN_S_sub <- LDSCoutput_EVEN$S[keep, keep]

Ssmooth <- as.matrix((nearPD(EVEN_S_sub, corr = FALSE))$mat)  # project to nearest positive-definite matrix
max(Ssmooth - EVEN_S_sub)  # 8.61421e-05 — smoothing difference introduced by nearPD
dim(Ssmooth)  # 11 11



GenomicSEM_RES <- list()
for(num_of_factors in 1:5){  # models with >5 factors fail to converge

  cat( paste0('current factor number: ', num_of_factors, ' -----------------\n') )

  ### Step 1: EFA (based on EVEN chr) ----------------------------
  # Using even chromosomes for EFA prevents overfitting: the factor structure
  # learned here is validated independently on odd chromosomes in the CFA step.
  EFA <- factanal(covmat = Ssmooth, factors = num_of_factors, rotation = "promax", start = NULL)
  L <- as.matrix(EFA$loadings)
  L_keep <- matrix(NA_real_, nrow(L), ncol(L), dimnames = dimnames(L))
  athreshold <- 0.3  # loading threshold following the GenomicSEM Nature paper convention
  for (i in seq_len(nrow(L))) { x <- L[i,]; a <- abs(x); keep <- which(a > athreshold)
    if (!length(keep)) keep <- which.max(a)
    L_keep[i, keep] <- x[keep] }


  ### Step 2: CFA (based on ODD chr) ------------------------------
  # Using odd chromosomes for CFA provides an independent test of the factor
  # structure identified by EFA on even chromosomes, avoiding overfitting.
  model_now <- make_lavaan_syntax(L_keep)
  for (iter in 1:10) {  # iterate to resolve negative residual variances (at most dim(Ssmooth) traits can be affected)
    fit <- usermodel(LDSCoutput_ODD, estimation = "DWLS", model = model_now, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
    res <- fit$results
    neg_trait_name <- res$lhs[res$lhs %in% rownames(L_keep) & res$Unstand_Est < 0]

    if (length(neg_trait_name) == 0) break  # no more negative residual variances; accept current model

    labs <- paste0("p", iter, "_", seq_along(neg_trait_name))
    add  <- paste(sprintf("%s ~~ %s*%s\n%s > .001", neg_trait_name, labs, neg_trait_name, labs), collapse = "\n")
    model_now <- paste(model_now, add, sep = "\n")
  }
  genomicSEMfit <- fit

  ### Store the i-th Genomic SEM result
  GenomicSEM_RES[[num_of_factors]] <- genomicSEMfit
}

fit_table <- rbindlist(lapply( 1:5, \(i) as.data.table(cbind(model_id=i, GenomicSEM_RES[[i]]$modelfit)) ), fill=TRUE)




GenomicSEM_varimax_RES <- list()
for(num_of_factors in 1:5){  # models with >5 factors fail to converge

  cat( paste0('current factor number: ', num_of_factors, ' -----------------\n') )

  ### Step 1: EFA (based on EVEN chr) ----------------------------
  EFA <- factanal(covmat = Ssmooth, factors = num_of_factors, rotation = "varimax", start = NULL)
  L <- as.matrix(EFA$loadings)
  L_keep <- matrix(NA_real_, nrow(L), ncol(L), dimnames = dimnames(L))
  for (i in seq_len(nrow(L))) { x <- L[i,]; a <- abs(x); keep <- which(a > 0.3)
    if (!length(keep)) keep <- which.max(a)
    L_keep[i, keep] <- x[keep] }


  ### Step 2: CFA (based on ODD chr) ------------------------------
  model_now <- make_lavaan_syntax(L_keep)
  for (iter in 1:10) {
    fit <- usermodel(LDSCoutput_ODD, estimation = "DWLS", model = model_now, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
    res <- fit$results
    neg_trait_name <- res$lhs[res$lhs %in% rownames(L_keep) & res$Unstand_Est < 0]

    if (length(neg_trait_name) == 0) break

    labs <- paste0("p", iter, "_", seq_along(neg_trait_name))
    add  <- paste(sprintf("%s ~~ %s*%s\n%s > .001", neg_trait_name, labs, neg_trait_name, labs), collapse = "\n")
    model_now <- paste(model_now, add, sep = "\n")
  }
  genomicSEMfit <- fit

  ### Store the i-th Genomic SEM result
  GenomicSEM_varimax_RES[[num_of_factors]] <- genomicSEMfit
}

fit_varimax_table <- rbindlist(lapply( 1:5, \(i) as.data.table(cbind(model_id=i, GenomicSEM_varimax_RES[[i]]$modelfit)) ), fill=TRUE)
# Varimax rotation model fit summary:
# model_id      chisq    df       p_chisq        AIC       CFI       SRMR
#        1 14075.8358    35  0.000000e+00 14115.8358 0.7892530 0.13349695
#        2 12086.3870    30  0.000000e+00 12136.3870 0.8190387 0.12300741
#        3   902.0262    27 1.213123e-172   958.0262 0.9868662 0.08594221
#        4   446.0240    21  2.634567e-81   514.0240 0.9936206 0.06105651
#        5   307.5032    15  1.539422e-56   387.5032 0.9956097 0.04922065
fit_table
# Promax rotation model fit summary:
# model_id      chisq    df       p_chisq        AIC       CFI       SRMR
#        1 14075.8358    35  0.000000e+00 14115.8358 0.7892530 0.13349695
#        2 12719.5199    33  0.000000e+00 12763.5199 0.8095807 0.13047558
#        3   668.0284    29 4.620671e-122   720.0284 0.9904085 0.09549580
#        4   476.1628    24  1.465025e-85   538.1628 0.9932132 0.07556298
#        5   197.7763    19  9.421731e-32   269.7763 0.9973166 0.05601299  # <-- chosen model for final visualisation

GenomicSEM_fit_comparision <- rbind(fit_varimax_table, fit_table)
varimax_row <- data.table(model_id = "Varimax", chisq = NA, df = NA, p_chisq = NA, AIC = NA, CFI = NA, SRMR = NA)
promax_row  <- data.table(model_id = "Promax",  chisq = NA, df = NA, p_chisq = NA, AIC = NA, CFI = NA, SRMR = NA)
result <- rbindlist(list(varimax_row, GenomicSEM_fit_comparision[1:5], promax_row, GenomicSEM_fit_comparision[6:10]), fill = TRUE)
GenomicSEM_fit_comparision <- result
saveRDS(GenomicSEM_fit_comparision, file = "/path/to/project/GenomicSEM/GenomicSEM_fit_comparision.rds")


### Model selection conclusion:
# 5-factor promax model (correlated factors) achieves the best performance by
# AIC and CFI on the independent validation set (ODD chromosomes).
GenomicSEM_RES[[5]]
semPlotModel_object <- create_semPlotModel( GenomicSEM_RES[[5]] )
semPaths(semPlotModel_object, layout="tree2", whatLabels = "est")




### ------------------------------------------------------------------------------------------------------
################# Final best Genomic model (refitted on *ALL* chromosomes; strictly speaking, autosomes) -----
num_of_factors <- 5

### Step 1: EFA (still based on EVEN chr — must match the structure used for model selection) ---
# Ssmooth was already computed above using nearPD on the EVEN-chromosome sub-matrix.
EFA <- factanal(covmat = Ssmooth, factors = num_of_factors, rotation = "promax", start = NULL)
L <- as.matrix(EFA$loadings)
L_keep <- matrix(NA_real_, nrow(L), ncol(L), dimnames = dimnames(L))
athreshold <- 0.3  # loading threshold following the GenomicSEM Nature paper convention
for (i in seq_len(nrow(L))) { x <- L[i,]; a <- abs(x); keep <- which(a > athreshold)
  if (!length(keep)) keep <- which.max(a)
  L_keep[i, keep] <- x[keep] }
L_keep


### Step 2: Final CFA (now based on ALL chromosomes) ------------------------------
model_now <- make_lavaan_syntax(L_keep)
for (iter in 1:10) {
  fit <- usermodel(LDSCoutput_ALL, estimation = "DWLS", model = model_now, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
  res <- fit$results
  neg_trait_name <- res$lhs[res$lhs %in% rownames(L_keep) & res$Unstand_Est < 0]

  if (length(neg_trait_name) == 0) break

  labs <- paste0("p", iter, "_", seq_along(neg_trait_name))
  add  <- paste(sprintf("%s ~~ %s*%s\n%s > .001", neg_trait_name, labs, neg_trait_name, labs), collapse = "\n")
  model_now <- paste(model_now, add, sep = "\n")
}
genomicSEMfit_ALL <- fit
genomicSEMfit_ALL  # AIC: 297.0126    CFI: 0.9977    SRMR: 0.05205
View( genomicSEMfit_ALL$results )
Final_GenomicSEM_fit <- genomicSEMfit_ALL$results
saveRDS(Final_GenomicSEM_fit, file = "/path/to/project/GenomicSEM/Final_GenomicSEM_fit.rds")


semPaths(create_semPlotModel( genomicSEMfit_ALL ), layout="tree2", whatLabels = "est")


# Extract standardised estimates and standard errors for reporting
summary(genomicSEMfit_ALL$results$STD_All - genomicSEMfit_ALL$results$STD_Genotype)  # verify STD_Genotype and STD_All are similar
est_se_info <- genomicSEMfit_ALL$results[, c('lhs','op','rhs','STD_Genotype','STD_Genotype_SE')]

myround <- function(x, n = 2) sprintf(paste0("%.", n, "f"), x)
est_se_info$est_se <- paste0( myround(est_se_info$STD_Genotype, 2), ' (', myround(as.numeric(est_se_info$STD_Genotype_SE), 2), ')' )

data.frame( est_se_info$lhs, est_se_info$op, est_se_info$rhs, est_se_info$est_se )




### Final visualisation -----------------------

### Pie plots: variance explained by each latent factor and the residual (error) term
# For features loading on multiple factors, only the Var(F) part is plotted for simplicity
# (the inter-factor covariance term is omitted).
model_for_plot <- genomicSEMfit_ALL$results[, c('lhs','op','rhs','STD_All','p_value')]

# Feature pie charts ---------------
for(current_feature in large_h2_features_plus_liver_fat_iron){
  # Residual variance component
  error_part         <- model_for_plot[ (model_for_plot$rhs==current_feature) & (model_for_plot$lhs==current_feature), ]
  error_term_variance <- error_part$STD_All
  # Factor loading components
  factor_part   <- model_for_plot[ (model_for_plot$rhs==current_feature) & (model_for_plot$lhs!=current_feature), ]
  factor_effect2 <- factor_part$STD_All^2  # for a single-factor feature: error_term_variance + factor_effect2 = 1

  # Proportions: slot 1 = error term; slots 2-6 = Factor 1-5
  proportion       <- c( error_term_variance, factor_effect2/sum(factor_effect2)*(1-error_term_variance) )
  final_proportion <- rep(0, 6)
  final_proportion[c(1, 1+as.integer(sub("Factor", "", factor_part$lhs)))] <- proportion  # insert in factor order so colours always match
  p1  <- final_proportion[1] / sum(final_proportion)
  ang <- 270 - 180 * p1  # position the error-term slice symmetrically around the bottom of the pie
  par(mar = c(5, 5, 0, 5))
  pie(final_proportion, labels = NA, init.angle = ang, radius=1)
  lab   <- gsub("_", "\n", current_feature)
  nline <- length(strsplit(lab, "\n", fixed = TRUE)[[1]])
  line_use <- 1.5 + (nline - 1) / 2
  mtext(lab, side = 1, line = line_use, cex = 2.0)

  # Save to PDF
  out_dir <- "/path/to/project/GenomicSEM/pie_plots"
  pdf(file = file.path(out_dir, paste0(current_feature, ".pdf")), width = 4, height = 4)
  par(mar = c(5, 5, 0, 5))
  pie(final_proportion, labels = NA, init.angle = ang, radius=1)
  lab   <- gsub("_", "\n", current_feature)
  nline <- length(strsplit(lab, "\n", fixed = TRUE)[[1]])
  line_use <- 1.2 + (nline - 1) / 2
  mtext(lab, side = 1, line = line_use, cex = 2.0)
  dev.off()
}




# Inter-factor correlation matrix (Factor1–Factor5)
factors     <- paste0('Factor', 1:5)
factor_part <- model_for_plot[ (model_for_plot$rhs %in% factors) & (model_for_plot$lhs %in% factors), ]

library(dplyr)

df <- factor_part
# Retain only off-diagonal correlations (lhs != rhs)
d <- df %>% filter(lhs != rhs) %>%
  mutate( z_abs = qnorm(1 - p_value/2),      # |Z| from two-sided p-value
          se    = abs(STD_All) / z_abs        # SE = |r| / |Z|
  )
facs <- c('Factor1','Factor5','Factor2','Factor3','Factor4')
# Correlation matrix
cor_mat <- matrix(NA_real_, length(facs), length(facs), dimnames=list(facs, facs))
diag(cor_mat) <- 1
idx <- cbind(match(d$lhs, facs), match(d$rhs, facs))
cor_mat[idx] <- d$STD_All; idx <- idx[,c(2,1)]; cor_mat[idx] <- d$STD_All  # fill symmetrically
# SE matrix
se_mat <- matrix(NA_real_, length(facs), length(facs), dimnames=list(facs, facs))
diag(se_mat) <- 0
se_mat[idx] <- d$se; idx <- idx[,c(2,1)]; se_mat[idx] <- d$se              # fill symmetrically
# Significance testing
z <- cor_mat / se_mat
p <- 2 * pnorm(-abs(z))
diag(p) <- 1

m    <- nrow(cor_mat) * (nrow(cor_mat) - 1) / 2
bonf <- 0.05 / m
corrplot::corrplot(cor_mat,
                   method="color", type="full", diag=TRUE,
                   p.mat=p, sig.level=bonf, insig="label_sig",
                   pch="*", pch.cex=1.2, pch.col="black",
                   addgrid.col="grey50", tl.col="black", tl.cex=0.8
)




# FDR-adjusted significance for path arrows
model_for_plot$FDRp    <- p.adjust(model_for_plot$p_value, method = "BH")
model_for_plot$FDRtrue <- model_for_plot$FDRp < 0.05
model_for_plot
# Note: Standard errors are not provided for the STD_All column because the
# sandwich-corrected SE components are not currently available in GenomicSEM.
# STD_All is expected to closely approximate STD_Genotype when std.lv = TRUE
# (i.e., latent factor variances are constrained to 1).




######### GenomicSEM column interpretation --------------------------------------------------
## Unstand_Est / Unstand_SE:
##   Neither the factor nor the observed variable is standardised.
## STD_Genotype / STD_Genotype_SE:
##   The observed variables are standardised (i.e., the genetic covariance matrix
##   is rescaled to a genetic correlation matrix with diagonal = 1). This is the
##   most commonly reported standardisation in GenomicSEM analyses.
## STD_All:
##   Both the factor and the observed variable are standardised (the classical
##   psychometric standardisation where latent factors also have unit variance).

## STD_Genotype and STD_Genotype_SE are produced by re-running the specified
## model using the standardised input (i.e., a genetic correlation and sampling
## correlation matrix). Results are therefore always standardised with respect
## to the genetic variance in the phenotypes.

## In the factor analysis interpretation (std.lv = TRUE):
##   "factor =~ variable"    — the loading is the factor effect on the indicator
##   "variable ~~ variable"  — the estimate is the residual (error) variance
##
## Under STD_All, for a single-factor indicator, variance decomposes as:
##   1 = [loading^2] + [residual variance]
##
## Pie chart proportions therefore use:
##   single-factor indicator: residual variance + loading^2
##   multi-factor indicator:  residual variance + Var(est1*F1 + ... + est5*F5)




########  GenomicSEM common problems and solutions:

### 1. Negative residual variance for some traits
### Solution: constrain the residual variance: trait ~~ a*trait, a > 0.0001

### 2. A smoothing difference > 0.025 for Z-statistics in the genetic covariance matrix
###    (pre- versus post-nearPD smoothing)
### Solution: exclude traits with non-significant h^2 (the problem originates
###    from weak genetic signal, not from the SEM model itself).
###    Note: this check occurs after trait selection according to the lavaan model
###    syntax, so the LDSC run does not need to be repeated.
S_LD  <- LDSCoutput_ODD$S
S_LDb <- LDSCoutput_ODD$S
smooth1  <- ifelse(eigen(S_LD)$values[nrow(S_LD)] <= 0, S_LD <- as.matrix((nearPD(S_LD, corr=FALSE))$mat), S_LD <- S_LD)
LD_sdiff <- max(abs(S_LD - S_LDb))  # S-matrix smoothing difference
V_LD  <- LDSCoutput_ODD$V
V_LDb <- LDSCoutput_ODD$V
smooth2   <- ifelse(eigen(V_LD)$values[nrow(V_LD)] <= 0, V_LD <- as.matrix((nearPD(V_LD, corr=FALSE))$mat), V_LD <- V_LD)
LD_sdiff2 <- max(abs(V_LD - V_LDb))  # V-matrix smoothing difference
k        <- nrow(S_LD)
SE_pre   <- matrix(0, k, k); SE_pre[lower.tri(SE_pre,  diag=TRUE)] <- sqrt(diag(V_LDb))
SE_post  <- matrix(0, k, k); SE_post[lower.tri(SE_post, diag=TRUE)] <- sqrt(diag(V_LD))
Z_pre  <- S_LDb / SE_pre   # Z-matrix using original S and V
Z_post <- S_LD  / SE_post  # Z-matrix using smoothed S and V
Z_diff <- Z_pre - Z_post
Z_diff[which(!is.finite(Z_diff))] <- 0
Z_diff <- max(Z_diff)  # maximum Z-statistic change due to smoothing
Z_diff  # 0.2083034

### 3. Covariance matrix of latent variables is not positive definite
### Solution: This indicates the fitted model implies a covariance matrix that
###    cannot exist (a valid covariance/correlation matrix must be positive
###    semidefinite). Consider simplifying the factor structure or adding
###    constraints.




####### Hierarchical (second-order / p-factor) model ------------------------------------------
# Fit a second-order factor model to characterise the inter-factor covariance
# structure among the five first-order genetic factors.

# Naive factor analysis on the first-order factor covariance matrix (1 further factor for initial inspection)
factor_covariance <- genomicSEMfit_ALL$results[ , c('lhs', 'rhs', 'STD_All') ]
factor_covariance <- factor_covariance[ factor_covariance$lhs %in% paste0('Factor', 1:5), ]
factor_covariance <- factor_covariance[ factor_covariance$rhs %in% paste0('Factor', 1:5), ]
factor_covariance  # covariance matrix of the 5 first-order latent factors

dt      <- as.data.table(factor_covariance)
dt      <- rbind(dt, dt[lhs != rhs, .(lhs = rhs, rhs = lhs, STD_All)], use.names = TRUE)
corr_mat <- as.matrix(dcast(dt, lhs ~ rhs, value.var = "STD_All")[, -1, with = FALSE], rownames = dcast(dt, lhs ~ rhs, value.var = "STD_All")$lhs)
eigen(corr_mat, symmetric=TRUE)$values

# Apply nearPD smoothing and run naive EFA on the first-order factor correlation matrix
FSsmooth <- as.matrix((nearPD(corr_mat, corr = TRUE))$mat)
FSsmooth - corr_mat  # verify smoothing difference is negligible
FAonFactors <- factanal(covmat = FSsmooth, factors = 1, rotation = "promax", start = NULL)
FAonFactors
# Result is rotation-invariant with a single further factor.
# Uniquenesses (residual variance of each first-order factor):
# Factor1 Factor2 Factor3 Factor4 Factor5
#   0.474   0.711   0.801   0.882   0.057
#
# Loadings on the single further factor:
#         further_Factor1
# Factor1   0.725          # 0.725^2 + 0.474 = 0.999 ≈ 1
# Factor2   0.537
# Factor3   0.446
# Factor4   0.344
# Factor5   0.971


# Identification note for correlation FA (diagonal fixed to 1) with K variables,
# using one further factor (K free loading parameters; error terms = 1 - loading^2):
# K=2: 1 observed covariance < 2 parameters -> under-identified
# K=3: 3 observed covariances = 3 parameters -> just-identified
# K=4: 6 observed covariances > 4 parameters -> over-identified
# K=5: 10 observed covariances > 5 parameters -> over-identified

# For exact identification with K variables: need (# of factors) = (K-1)/2,
# so for K=5 the exact-identification model uses 2 further factors.

FAonFactors <- factanal(covmat = FSsmooth, factors = 2, rotation = "varimax", start = NULL)
FAonFactors
# Loadings:
#         Factor1 Factor2
# Factor1  0.997             # Factor1 is dominated by further factor 1
# Factor2  0.585   0.551
# Factor3          0.275
# Factor4  0.123   0.990     # Factor4 is dominated by further factor 2
# Factor5  0.704   0.236
lambda       <- FAonFactors$loadings[, 1:2]
uniquenesses <- FAonFactors$uniquenesses
# Verify the FA-implied covariance matrix matches FSsmooth
tcrossprod(lambda) + diag(uniquenesses)  # should be close to FSsmooth
FSsmooth

uniq_mat <- matrix(c("Uniquenesses",rep(NA,4), names(FAonFactors$uniquenesses), round(FAonFactors$uniquenesses,3),NA), nrow=3, ncol=5, byrow=TRUE)
load_mat <- matrix(c("Loadings",rep(NA,4), NA,"Factor1","Factor2",NA,NA, "Factor1",round(FAonFactors$loadings[1,],3),NA,NA, "Factor2",round(FAonFactors$loadings[2,],3),NA,NA, "Factor3",round(FAonFactors$loadings[3,],3),NA,NA, "Factor4",round(FAonFactors$loadings[4,],3),NA,NA, "Factor5",round(FAonFactors$loadings[5,],3),NA,NA), nrow=8, ncol=5, byrow=TRUE)
fit_mat  <- matrix(c("Fit Statistics",rep(NA,4), NA,"Factor1","Factor2",NA,NA, "SS loadings",1.851,1.415,NA,NA, "Proportion Var",0.370,0.283,NA,NA, "Cumulative Var",0.370,0.653,NA,NA), nrow=5, ncol=5, byrow=TRUE)
FAres <- as.data.table(rbind(uniq_mat, rep(NA,5), load_mat, rep(NA,5), fit_mat))
View(FAres)
saveRDS(FAres, file = "/path/to/project/GenomicSEM/naiveFAres.rds")




# Second-order CFA within GenomicSEM (free-parameter approach)
# Factor1 and Factor4 each anchor one second-order factor;
# remaining factors load on both higher-order factors.
model_2nd_nofixed <- paste0(
  model_now,
  sprintf("
# ---------- second-order factor (free loadings; anchored by Factor1 and Factor4) ----------
Higher1 =~ 1*Factor1 + Factor2 + Factor3  + Factor5
Higher4 =~ 1*Factor4  + Factor2 + Factor3 + Factor5
Factor1 ~~ 0*Factor1   # residual variance set to zero (fully explained by Higher1)
Factor4 ~~ 0*Factor4

# ---------- identification ----------
Higher1 ~~ 1*Higher1
Higher4 ~~ 1*Higher4
Higher1 ~~ 0*Higher4   # second-order factors assumed orthogonal

") )

fit_2nd <- usermodel(LDSCoutput_ALL, estimation="DWLS", model=model_2nd_nofixed, std.lv=FALSE, CFIcalc=TRUE, imp_cov=FALSE)
res <- fit_2nd$results
neg_trait_name <- res$lhs[res$lhs %in% rownames(L_keep) & res$Unstand_Est < 0]
labs <- paste0("pp", 1, "_", seq_along(neg_trait_name))
add  <- paste(sprintf("%s ~~ %s*%s\n%s > .001", neg_trait_name, labs, neg_trait_name, labs), collapse = "\n")
model_2nd_nofixed <- paste(model_2nd_nofixed, add, sep = "\n")
fit_2nd <- usermodel(LDSCoutput_ALL, estimation="DWLS", model=model_2nd_nofixed, std.lv=FALSE, CFIcalc=TRUE, imp_cov=FALSE)

genomicSEMfit_ALL_second_order <- fit_2nd
genomicSEMfit_ALL_second_order  # AIC=545.3352  CFI=0.9949798  SRMR=0.07312367
# CFI and SRMR indicate acceptable fit; results are broadly consistent with the
# first-order model, though the additional second-order structure adds complexity.
Final_GenomicSEM_fit_second_order <- genomicSEMfit_ALL_second_order$results
saveRDS(Final_GenomicSEM_fit_second_order, file = "/path/to/project/GenomicSEM/Final_GenomicSEM_fit_second_order.rds")


semPaths(create_semPlotModel( genomicSEMfit_ALL_second_order ), layout="tree2", whatLabels = "est")
# Verification: F2 variance = F2_residual + Higher1->F2 loading^2 + Higher4->F2 loading^2
# = 0.40348940 + 0.53751593^2 + 0.55460547^2 ≈ 1


# Pie-chart visualisation for second-order model (each first-order factor as the "feature")
model_for_plot <- genomicSEMfit_ALL_second_order$results[, c('lhs','op','rhs','STD_All','p_value')]

for(current_feature in paste0('Factor', 1:5)){
  # Residual variance component
  error_part          <- model_for_plot[ (model_for_plot$rhs==current_feature) & (model_for_plot$lhs==current_feature), ]
  error_term_variance <- error_part$STD_All
  # Second-order factor components
  factor_part   <- model_for_plot[ (model_for_plot$rhs==current_feature) & (model_for_plot$lhs!=current_feature), ]
  factor_effect2 <- factor_part$STD_All^2

  proportion       <- c( error_term_variance, factor_effect2/sum(factor_effect2)*(1-error_term_variance) )
  final_proportion <- rep(0, length(1:2) + 5)  # 5 first-order factors + 2 second-order factors
  final_proportion[c( as.integer(sub("Factor", "", unique(factor_part$rhs))), 5+1:2 )] <- proportion
  p1  <- final_proportion[ as.integer(sub("Factor", "", unique(factor_part$rhs))) ] / sum(final_proportion)
  ang <- 270 - 180 * p1
  final_colors <- c("lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk", "lightblue", "lavender")
  par(mar = c(5, 5, 0, 5))
  pie(final_proportion, col=final_colors, labels = NA, init.angle = ang, radius=1)

  out_dir <- "/path/to/project/GenomicSEM/pie_plots"
  pdf(file = file.path(out_dir, paste0(current_feature, ".pdf")), width = 4, height = 4)
  final_colors <- c("lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk", "lightblue", "lavender")
  par(mar = c(5, 5, 0, 5))
  pie(final_proportion, col=final_colors, labels = NA, init.angle = ang, radius=1)
  dev.off()
}


# Factor pie charts (single-colour fill for each factor)
for(ff in 1:5){
  out_dir <- "/path/to/project/GenomicSEM/pie_plots"
  pdf(file = file.path(out_dir, paste0('F', ff, ".pdf")), width = 4, height = 4)
  final_proportion <- rep(0, 6)  # residual term + Factor 1-5
  final_proportion[1+ff] <- 1
  pie(final_proportion, labels = NA, border=0, lty=0)
  tt <- seq(0, 2*pi, length.out = 400)
  lines(r*cos(tt), r*sin(tt), col = "black", lwd = 1.5)
  dev.off()
}




### ----------------------------------------------------------------------------
# Deprecated: fixed-parameter second-order Genomic SEM approach
# (loadings and disturbances fixed to values derived from the naive FA above;
#  superseded by the free-parameter second-order model fitted above)
# ----------------------------------------------------------------------------
model_2nd_fixed <- paste0(
  model_now,
  sprintf("
# ---------- second-order factor (fixed loadings from FA) ----------
Higher =~ %0.6f*Factor1 + %0.6f*Factor2 + %0.6f*Factor3 +
          %0.6f*Factor4 + %0.6f*Factor5

# ---------- identification ----------
Higher ~~ 1*Higher


# ---------- fixed first-order disturbances (from FA uniquenesses) ----------
Factor1 ~~ %0.6f*Factor1
Factor2 ~~ %0.6f*Factor2
Factor3 ~~ %0.6f*Factor3
Factor4 ~~ %0.6f*Factor4
Factor5 ~~ %0.6f*Factor5
",
          lambda["Factor1"], lambda["Factor2"], lambda["Factor3"],
          lambda["Factor4"], lambda["Factor5"],
          uniquenesses["Factor1"], uniquenesses["Factor2"], uniquenesses["Factor3"],
          uniquenesses["Factor4"], uniquenesses["Factor5"]
  ) )


model_2nd_fixed <- paste0(
  model_now,
  sprintf("
# ---------- second-order factor (fixed loadings from FA) ----------
Higher1 =~ %0.6f*Factor1 + %0.6f*Factor2 + %0.6f*Factor3 +%0.6f*Factor4 + %0.6f*Factor5
Higher2 =~ %0.6f*Factor1 + %0.6f*Factor2 + %0.6f*Factor3 +%0.6f*Factor4 + %0.6f*Factor5

# ---------- identification ----------
Higher1 ~~ 1*Higher1
Higher2 ~~ 1*Higher2
Higher1 ~~ 0*Higher2

# ---------- fixed first-order disturbances (from FA uniquenesses) ----------
Factor1 ~~ %0.6f*Factor1
Factor2 ~~ %0.6f*Factor2
Factor3 ~~ %0.6f*Factor3
Factor4 ~~ %0.6f*Factor4
Factor5 ~~ %0.6f*Factor5
",
  lambda["Factor1",1], lambda["Factor2",1], lambda["Factor3",1], lambda["Factor4",1], lambda["Factor5",1],
  lambda["Factor1",2], lambda["Factor2",2], lambda["Factor3",2], lambda["Factor4",2], lambda["Factor5",2],
          uniquenesses["Factor1"], uniquenesses["Factor2"], uniquenesses["Factor3"],
          uniquenesses["Factor4"], uniquenesses["Factor5"]
  ) )


fit_2nd <- usermodel(LDSCoutput_ALL, estimation="DWLS", model=model_2nd_fixed, std.lv=FALSE, CFIcalc=TRUE, imp_cov=FALSE)
res <- fit_2nd$results
neg_trait_name <- res$lhs[res$lhs %in% rownames(L_keep) & res$Unstand_Est < 0]
labs <- paste0("pp", 1, "_", seq_along(neg_trait_name))
add  <- paste(sprintf("%s ~~ %s*%s\n%s > .001", neg_trait_name, labs, neg_trait_name, labs), collapse = "\n")
model_2nd_fixed <- paste(model_2nd_fixed, add, sep = "\n")
fit_2nd <- usermodel(LDSCoutput_ALL, estimation="DWLS", model=model_2nd_fixed, std.lv=FALSE, CFIcalc=TRUE, imp_cov=FALSE)

genomicSEMfit_ALL_second_order <- fit_2nd
genomicSEMfit_ALL_second_order
# CFI and SRMR indicate acceptable fit; broadly consistent with the first-order model.
semPaths(create_semPlotModel( genomicSEMfit_ALL_second_order ), layout="tree2", whatLabels = "est")


# Visualisation for fixed-parameter second-order model (factor-level pie charts)
model_for_plot <- genomicSEMfit_ALL_second_order$results[, c('lhs','op','rhs','STD_All','p_value')]

for(current_feature in paste0('Factor', 1:5)){
  error_part          <- model_for_plot[ (model_for_plot$rhs==current_feature) & (model_for_plot$lhs==current_feature), ]
  error_term_variance <- error_part$STD_All
  factor_part   <- model_for_plot[ (model_for_plot$rhs==current_feature) & (model_for_plot$lhs!=current_feature), ]
  factor_effect2 <- factor_part$STD_All^2

  proportion       <- c( error_term_variance, factor_effect2/sum(factor_effect2)*(1-error_term_variance) )
  final_proportion <- rep(0, 6)  # slot 1: residual; slots 2-6: Factor 1-5
  final_proportion[c( 1+as.integer(sub("Factor", "", factor_part$rhs)), 1 )] <- proportion
  p1  <- final_proportion[1] / sum(final_proportion)
  ang <- 90 - 180 * p1
  final_colors <- c("grey95", "lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk")
  par(mar = c(5, 5, 0, 5))
  pie(final_proportion, col=final_colors, labels = NA, init.angle = ang, radius=1)

  out_dir <- "/path/to/project/GenomicSEM/pie_plots"
  pdf(file = file.path(out_dir, paste0(current_feature, ".pdf")), width = 4, height = 4)
  final_colors <- c("grey95", "lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk")
  par(mar = c(5, 5, 0, 5))
  pie(final_proportion, col=final_colors, labels = NA, init.angle = ang, radius=1)
  dev.off()
}


# Pie colour reference: c("white", "lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk")
c("white", "lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk")
