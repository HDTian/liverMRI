# =============================================================================
# Author: Haodong Tian
# Description: Post-GWAS analysis and visualization: LDSC intercept/GC/QQ plots,
#              heritability summary, genetic correlation eigen-decomposition,
#              bivariate LDSC between MRI features and liver fat/cT1,
#              stratified S-LDSC cell-type enrichment, locus summary with
#              PheWAS annotation, and PDFF-lead variant effect scatter plots.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

library(data.table)
library(qqman)
library(CMplot)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggsignif)


### LDSC intercept, genomic inflation factor (Lambda_GC), QQ plots --------------
### -----------------------------------------------------------------------------

# Univariate LDSC heritability results for all 200 MRI features
LDSC <- read.table('/path/to/server/project/GWAS_regenie/LDSC_new/heritability_summary.tsv',
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dim(LDSC)
LDSC_res <- LDSC; LDSC_res$Z <- LDSC_res$h2 / LDSC_res$h2_se  # univariate LDSC Z score

summary(LDSC_res$Intercept)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9887  1.0122  1.0199  1.0180  1.0244  1.0343

summary(LDSC_res$Lambda_GC)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.022   1.053   1.071   1.076   1.096   1.152


# The 59 retained features — features passing heritability and correlation QC
retained_names <- readRDS("/path/to/project/phenotype data/features with other images/retained_names.rds")

LDSC_retained <- LDSC_res[LDSC_res$Feature %in% retained_names, ]
dim(LDSC_retained)
LDSC_retained <- LDSC_retained[order(-LDSC_retained$Z), ]
saveRDS(LDSC_retained, file = "/path/to/project/GWAS_res/LDSC_retained.rds")

# QQ plots for the 10 most heritable retained features
qqplot_names <- LDSC_retained$Feature[order(-LDSC_retained$h2)][1:10]

indir <- "/path/to/server/project/GWAS_regenie/MR_GWAS_new"
setwd('/path/to/project/GWAS_res/QQ_pplots')

for (nm in qqplot_names) {
  cat('==================================\n')
  cat(paste0(nm, '\n'))

  f <- file.path(indir, paste0(nm, ".regenie"))
  d <- fread(f, select = c("P"))

  d_for_plot <- data.frame(
    SNP     = paste0("rs", 1:length(d$P)),
    Chr     = 1,
    Pos     = 1:length(d$P),
    P_value = d$P
  )

  summary(d$P)
  lambda_GC      <- median(qchisq(1 - d$P, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
  LDSC_lambda_GC <- LDSC_retained$Lambda_GC[LDSC_retained$Feature == nm]
  # Note: LDSC Lambda_GC uses HapMap3 QC-filtered SNPs, giving a slightly different value
  cat(paste0('lambda_GC: ', round(lambda_GC, 3), ' (from full p-values) and ',
             round(LDSC_lambda_GC, 3), ' (LDSC)\n'))

  CMplot(d_for_plot, plot.type = "q",
         main      = paste0(nm, "\nlambda_GC = ", round(LDSC_lambda_GC, 3)),
         threshold = 5e-8,
         conf.int  = TRUE,
         box       = FALSE, file.output = TRUE,
         file.name = nm,
         file      = "jpg",
         dpi       = 300)
}


### 200x200 genetic correlation matrix -> 59 retained features -> eigen decomposition --------
### ------------------------------------------------------------------------------------------

df <- fread("/path/to/server/project/GWAS_regenie/LDSC_new/genetic_correlation.tsv", sep = "\t")
dim(df)  # 19900 = 200*199/2; lower triangle of the 200x200 pairwise rg matrix

summary(abs(df$Genetic_Correlation))
summary(df$Genetic_Correlation^2)
hist(abs(df$Genetic_Correlation))

rg_vectors <- df$Genetic_Correlation
rg_vectors[rg_vectors < -1] <- -1  # clip |rg| > 1 (occasional numerical overflow from LDSC)
rg_vectors[rg_vectors >  1] <-  1
summary(rg_vectors)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -1.0000 -0.3201  0.0225  0.0290  0.3863  1.0000
hist(rg_vectors, n = 20, probability = TRUE, col = "steelblue",
     xlab = 'Genetic correlation', main = '')


# Restrict to the 59 retained features
df_retained_names <- df[(df$Feature1 %in% retained_names) & (df$Feature2 %in% retained_names)]
dim(df_retained_names)  # 1711 = 59*58/2
rg_retained_vectors <- df_retained_names$Genetic_Correlation
summary(rg_retained_vectors)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -0.89840 -0.22385  0.02170  0.02652  0.27425  0.89310
hist(rg_retained_vectors, n = 20, probability = TRUE, xlab = 'Genetic correlation', main = '')


# Eigen-decomposition of the 59x59 genetic correlation matrix
traits <- sort(retained_names)
rg     <- matrix(1, length(traits), length(traits), dimnames = list(traits, traits))
dim(rg)  # 59 x 59
i1 <- match(df_retained_names$Feature1, traits); i2 <- match(df_retained_names$Feature2, traits)
rg[cbind(i1, i2)] <- df_retained_names$Genetic_Correlation
rg[cbind(i2, i1)] <- df_retained_names$Genetic_Correlation
rg[rg < -1] <- -1; rg[rg > 1] <- 1
e <- eigen(rg, symmetric = TRUE)

k           <- 20
total_var   <- sum(e$values[1:k])
df_eig      <- data.frame(idx = 1:k,
                           ev  = e$values[1:k],
                           cum = cumsum(e$values[1:k]) / total_var)
threshold   <- 0.8
cutoff_idx  <- min(which(df_eig$cum >= threshold))
scale_factor <- max(df_eig$ev)

ggplot(df_eig, aes(x = idx)) +
  geom_bar(aes(y = ev), stat = "identity", fill = "grey60", width = 0.7) +
  geom_line(aes(y = cum * scale_factor), color = "black", size = 0.8) +
  geom_point(aes(y = cum * scale_factor), shape = 21, fill = "white", size = 2) +
  annotate("segment", x = cutoff_idx, xend = cutoff_idx,
           y = 0, yend = df_eig$cum[cutoff_idx] * scale_factor,
           linetype = "dashed", color = "royalblue") +
  annotate("segment", x = 0, xend = cutoff_idx,
           y = df_eig$cum[cutoff_idx] * scale_factor,
           yend = df_eig$cum[cutoff_idx] * scale_factor,
           linetype = "dashed", color = "royalblue") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "firebrick", alpha = 0.8) +
  scale_y_continuous(
    name     = "Eigenvalue",
    sec.axis = sec_axis(~./scale_factor, name = "Cumulative Variance", labels = scales::percent),
    expand   = expansion(mult = c(0, 0.1))
  ) +
  scale_x_continuous(breaks = seq(1, k, by = 1)) +
  labs(x = "Dimension (Rank)") + theme_classic()



### Heritability of the 59 retained MRI features --------------------------------
### -----------------------------------------------------------------------------

summary(LDSC_retained$h2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.02230 0.06085 0.07660 0.09046 0.11405 0.21310

# Reference values from published LDSC analyses:
# liver fat h2:  0.1642 (SE 0.0188)
# liver iron h2: 0.1282 (SE 0.0222)

qqplot_names <- LDSC_retained$Feature[order(-LDSC_retained$h2)][1:10]


### Bivariate LDSC between 10 most heritable MRI features and liver fat / liver cT1 ----------

library(readxl)
library(GenomicSEM)

liver_features_munged_folder <- '/path/to/server/project/GWAS_regenie/LDSC_new/munged_files'
liver_feature_files  <- list.files(liver_features_munged_folder, pattern = "\\.sumstats\\.gz$",
                                    full.names = TRUE, recursive = FALSE)
liver_feature_names  <- sub("\\.sumstats\\.gz$", "", basename(liver_feature_files))

retained_names <- qqplot_names  # restrict to the 10 most heritable features for bLDSC
keep_idx <- liver_feature_names %in% retained_names
liver_feature_files <- liver_feature_files[keep_idx]
liver_feature_names <- liver_feature_names[keep_idx]

pathname  <- '/path/to/project'
LDSC_path <- file.path(pathname, 'LDSC', 'eur_w_ld_chr')

trait_info_table <- read_excel(paste0(pathname, '/trait_info.xlsx')); setDT(trait_info_table)
trait_info_table <- trait_info_table[!is.na(`file name`)]
trait_info_table <- trait_info_table[9:10]  # liver fat (row 9) and liver cT1/iron (row 10)

quiet <- function(expr) {
  out <- file(nullfile(), "wt"); msg <- file(nullfile(), "wt")
  on.exit({ try(sink(type = "message")); try(close(msg)); try(sink()); try(close(out)) }, add = TRUE)
  sink(out); sink(msg, type = "message")
  suppressWarnings(suppressMessages(force(expr)))
}

setwd(file.path(pathname, 'GWAS_res', 'bLDSC'))

for (kkk in 1:nrow(trait_info_table)) {
  current_trait <- sub("\\.(tsv|txt|csv)(\\.gz)?$", "", trait_info_table$`file name`[kkk], ignore.case = TRUE)
  out_file <- file.path(pathname, 'LDSC', "munged_GWAS", paste0(current_trait, ".sumstats.gz"))

  print('================================================')
  print(paste0('current phenotype:', current_trait))

  # Use sample_prev = 0.5 for binary traits as recommended by GenomicSEM documentation
  sample_prev_used <- if (trait_info_table$trait_type[kkk] == "binary") 0.5 else NA
  print(paste0('sample_prev_used: ', sample_prev_used))

  for (l in 1:length(liver_feature_names)) {
    cat(paste0(l, '-'))
    LDSCout <- quiet(ldsc(
      traits          = c(out_file, liver_feature_files[l]),
      trait.names     = c(current_trait, liver_feature_names[l]),
      sample.prev     = c(sample_prev_used, NA),
      population.prev = c(NA, NA),
      ld              = LDSC_path,
      wld             = LDSC_path,
      ldsc.log        = paste0(current_trait, '_vs_', liver_feature_names[l]),
      stand           = TRUE
    ))
  }
  cat('\n')
}

## Parse bLDSC log files to extract rg estimates and SE
logs_dir    <- file.path(pathname, 'GWAS_res', 'bLDSC')
trait_names <- sub("\\.(tsv|txt|csv)(\\.gz)?$", "", trait_info_table$`file name`, ignore.case = TRUE)
rg_mat <- matrix(NA_real_, length(liver_feature_names), length(trait_names),
                  dimnames = list(liver_feature_names, trait_names))
se_mat <- rg_mat

logs <- list.files(logs_dir, pattern = "_ldsc\\.log$", full.names = TRUE)
parse_log <- function(f) {
  x <- readLines(f, warn = FALSE); i <- grep("^Genetic Correlation between ", x)
  if (!length(i)) return(NULL)
  s <- x[tail(i, 1)]
  m <- regexec("^Genetic Correlation between\\s+(.+?)\\s+and\\s+(.+?):\\s*([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\\s*\\(([^)]+)\\)", s)
  r <- regmatches(s, m)[[1]]; if (length(r) < 5) return(NULL)
  list(a = trimws(r[2]), b = trimws(r[3]), rg = as.numeric(r[4]), se = as.numeric(r[5]))
}
who_is_who <- function(a, b, f) {
  if (a %in% liver_feature_names && b %in% trait_names) return(list(feat = a, trait = b))
  if (b %in% liver_feature_names && a %in% trait_names) return(list(feat = b, trait = a))
  z <- strsplit(sub("_ldsc\\.log$", "", basename(f)), "_vs_")[[1]]
  if (length(z) == 2) {
    if (z[1] %in% trait_names && z[2] %in% liver_feature_names) return(list(feat = z[2], trait = z[1]))
    if (z[2] %in% trait_names && z[1] %in% liver_feature_names) return(list(feat = z[1], trait = z[2]))
  }
  NULL
}
for (f in logs) {
  pr <- parse_log(f)
  if (is.null(pr) || any(!is.finite(c(pr$rg, pr$se)))) next
  id <- who_is_who(pr$a, pr$b, f)
  if (is.null(id)) next
  rg_mat[id$feat, id$trait] <- pr$rg
  se_mat[id$feat, id$trait] <- pr$se
}

dim(rg_mat); dim(se_mat)  # 10 x 2

LDSC_res_sub <- LDSC_res[LDSC_res$Feature %in% rownames(rg_mat), ]
out_tab <- data.table(
  MRI                 = LDSC_res_sub$Feature,
  heritability        = sprintf("%.3f", LDSC_res_sub$h2),
  `rg with liver fat` = sprintf("%.3f (%.3f)", rg_mat[LDSC_res_sub$Feature, "liver_fat"],
                                                se_mat[LDSC_res_sub$Feature, "liver_fat"]),
  `rg with liver cT1` = sprintf("%.3f (%.3f)", rg_mat[LDSC_res_sub$Feature, "liver_iron"],
                                                se_mat[LDSC_res_sub$Feature, "liver_iron"])
)
out_tab
saveRDS(out_tab, file = "/path/to/project/GWAS_res/the10heritableMRI.rds")



### Heritability barplot — 10 most heritable MRI features + liver fat / liver cT1 -----------

heritable_table <- data.frame(
  names  = qqplot_names,
  h2_est = LDSC_retained$h2[match(qqplot_names, LDSC_retained$Feature)],
  h2_se  = LDSC_retained$h2_se[match(qqplot_names, LDSC_retained$Feature)]
)
heritable_table <- rbind(heritable_table, c('liver_fat',  0.1642, 0.0188))
heritable_table <- rbind(heritable_table, c('liver_iron', 0.1282, 0.0222))
heritable_table$h2_est <- as.numeric(heritable_table$h2_est)
heritable_table$h2_se  <- as.numeric(heritable_table$h2_se)

d <- heritable_table[order(-heritable_table$h2_est), ]
d <- d[nrow(d):1, ]
d$names[d$names == 'liver_iron'] <- 'liver_cT1'

par(mar = c(4, 15, 2, 1))
bp <- barplot(d$h2_est,
              horiz     = TRUE,
              names.arg = d$names,
              las       = 1,
              xlim      = c(0, max(d$h2_est + 1.96*d$h2_se)),
              col       = "steelblue",
              border    = "black",
              lwd       = 0.6,
              xlab      = "SNP heritability (h²)")
arrows(d$h2_est - 1.96*d$h2_se, bp, d$h2_est + 1.96*d$h2_se, bp,
       angle = 90, code = 3, length = 0.04)
box()
par(mar = c(5.1, 4.1, 4.1, 2.1))



### Heritability vs. number of independent lead SNPs ----------------------------
### -----------------------------------------------------------------------------

tsv_dir <- "/path/to/project/new_UVMR/retained_UVMR_summary_1e6_new"

snp_counts <- sapply(retained_names, function(feat) {
  fpath <- file.path(tsv_dir, paste0(feat, ".tsv"))
  if (!file.exists(fpath)) return(0L) else { dt <- fread(fpath); return(nrow(dt)) }
})
snp_df <- data.frame(Feature = retained_names, snp_count = as.integer(snp_counts),
                      stringsAsFactors = FALSE)

snp_counts_5e8 <- sapply(retained_names, \(feat)
  if (file.exists(f <- file.path(tsv_dir, paste0(feat, ".tsv")))) sum(fread(f)$P < 5e-8) else 0L)
snp_df_5e8 <- data.frame(Feature = retained_names, snp_count_5e8 = snp_counts_5e8)

LDSC_retaind_final <- Reduce(function(x, y) merge(x, y, by = "Feature", all.x = TRUE),
                              list(LDSC_retained, snp_df, snp_df_5e8))
# firstorder_TotalEnergy has the largest h2, # of SNPs (1e-6), and # of SNPs (5e-8)

library(ggplot2)

df <- LDSC_retaind_final
df$feature_class <- sub("_.*", "", df$Feature)
df$feature_class <- ifelse(df$feature_class %in% c("firstorder", "shape"), df$feature_class, "texture")
df$contrast      <- ifelse(grepl("_inp$", df$Feature), "In-phase", "Water")
table(df$feature_class)  # 12 firstorder; 7 shape; 40 texture
table(df$contrast)       # 27 in-phase; 32 water

cols <- c(firstorder = "darkred", shape = "darkblue", texture = "darkgreen")

ggplot(df, aes(x = snp_count_5e8, y = h2, color = feature_class, shape = contrast)) +
  geom_point(size = 2.6) +
  scale_color_manual(values = cols) +
  labs(x = "Number of independent SNPs", y = "Heritability",
       color = "Feature class", shape = "MRI contrast") +
  theme_classic() +
  theme(legend.position = "right", legend.box = "vertical")


### Water vs. in-phase comparison
df_long <- pivot_longer(df, cols = c(h2, snp_count_5e8),
                         names_to = "metric", values_to = "value")
df_long$metric <- ifelse(df_long$metric == 'h2', 'Heritability (h2)', 'Count of lead SNPs')

ggplot(df_long, aes(x = contrast, y = value, fill = contrast)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1.2) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  geom_signif(comparisons = list(c("Water", "In-phase")),
              test = "t.test",
              map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
              textsize = 4, vjust = -0.2, tip_length = 0) +
  scale_fill_manual(values = c("Water" = "#377eb8", "In-phase" = "#e41a1c")) +
  theme_bw() +
  labs(x = NULL, y = NULL) +
  theme(legend.position  = "none",
        strip.text       = element_text(face = "bold", size = 12),
        axis.text.x      = element_text(size = 11, color = "black"),
        panel.grid.minor = element_blank(),
        panel.spacing    = unit(1.5, "lines"))


### Feature class comparison (firstorder vs. shape vs. texture)
df_long <- pivot_longer(df, cols = c(h2, snp_count_5e8),
                         names_to = "metric", values_to = "value")
df_long$metric <- ifelse(df_long$metric == 'h2', 'Heritability (h2)', 'Count of lead SNPs')

ggplot(df_long, aes(x = feature_class, y = value, fill = feature_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  geom_signif(comparisons     = combn(unique(df_long$feature_class), 2, simplify = FALSE),
              step_increase    = 0.1, test = "t.test",
              map_signif_level = TRUE, tip_length = 0) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(x = NULL, y = NULL) +
  theme(legend.position  = "none",
        strip.text       = element_text(face = "bold", size = 12),
        axis.text.x      = element_text(angle = 30, hjust = 1, size = 10),
        panel.grid.minor = element_blank())



### Locus summary and PheWAS annotation -----------------------------------------
### -----------------------------------------------------------------------------

# Lead SNP coordinates for all 59 retained features (combined before additional clumping)
UVMRcorrds <- read.table(paste0(pathname, "/new_UVMR/UVMRcoords_new.txt"),
                          header = FALSE, sep = " ", stringsAsFactors = FALSE)
dim(UVMRcorrds)  # 461 x 3
rsID_list <- UVMRcorrds$V3

# LD clumping at r2 < 0.2 within 1 Mb to define independent loci
dummy_sumstats <- data.table(SNP = rsID_list, P = 0.001)
dummy_file     <- file.path(pathname, 'new_UVMR', 'PLINK_files', "dummy_sumstats.txt")
fwrite(dummy_sumstats, file = dummy_file, sep = "\t", quote = FALSE)

plink_path <- file.path(pathname, "PheWAS/plink_mac_20250615/plink")
bfile_path <- file.path(pathname, "PheWAS/1000G_QC")
out_prefix <- file.path(pathname, 'new_UVMR', 'PLINK_files', "rsID_clumped_results_for_loucs")

plink_cmd <- sprintf('"%s" --bfile "%s" --clump "%s" --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 1000 --out "%s"',
                     plink_path, bfile_path, dummy_file, out_prefix)
system(plink_cmd, intern = TRUE)
clumped_results <- fread(paste0(out_prefix, ".clumped"))
final_rsID_list <- clumped_results$SNP
length(final_rsID_list)  # 327 independent loci


# Locus-level summary for the 10 most heritable MRI features + liver fat
GWAS_5e8_final <- fread("/path/to/project/GWAS_res/input.txt",
                         sep = "\t", na.strings = "NA", quote = "")
dim(GWAS_5e8_final)
length(unique(GWAS_5e8_final$LOCUS_ID))  # 57 loci
length(unique(GWAS_5e8_final$LOCUS_ID[GWAS_5e8_final$TRAIT == 'liver_fat']))  # 10 loci
View(GWAS_5e8_final)


# PheWAS via NHGRI-EBI GWAS Catalog (gwasrapidd) for the 9 lead MRI feature variants
library(gwasrapidd)

the9leadMRI <- GWAS_5e8_final[GWAS_5e8_final$TRAIT != 'liver_fat', ]
my_snps     <- unique(the9leadMRI$MARKER)
length(my_snps)  # 77 unique independent SNPs

PheWAS_list <- list()
for (my_snp in my_snps) {
  cat('-')
  NHGRI_EBI_association_res <- get_associations(variant_id = my_snp)

  if (length(NHGRI_EBI_association_res@associations$association_id) != 0) {
    ids       <- NHGRI_EBI_association_res@associations$association_id
    trait_obj <- get_traits(association_id = ids)
    PheWAS_list[[my_snp]] <- trait_obj
  } else {
    PheWAS_list[[my_snp]] <- NULL
  }
}  # ~3 hours

saveRDS(PheWAS_list, file = "/path/to/project/GWAS_res/PheWAS/PheWAS_list_9leadMRI.rds")
# PheWAS_list <- readRDS("/path/to/project/GWAS_res/PheWAS/PheWAS_list_9leadMRI.rds")

PheWAS_names_vec <- c()
num_PheWAS_vec   <- c()
for (i in 1:nrow(the9leadMRI)) {
  rsID          <- the9leadMRI$MARKER[i]
  NHGRI_EBI_res <- PheWAS_list[[rsID]]
  if (is.null(NHGRI_EBI_res)) {
    PheWAS_names_vec <- c(PheWAS_names_vec, 'NULL')
    num_PheWAS_vec   <- c(num_PheWAS_vec, 0)
  } else {
    PheWAS_names_vec <- c(PheWAS_names_vec, paste(unique(NHGRI_EBI_res@traits$trait), collapse = " | "))
    num_PheWAS_vec   <- c(num_PheWAS_vec,   length(unique(NHGRI_EBI_res@traits$trait)))
  }
}
the9leadMRI$PheWASnames   <- PheWAS_names_vec
the9leadMRI$num_of_PheWAS <- num_PheWAS_vec

saveRDS(the9leadMRI, file = "/path/to/project/GWAS_res/PheWAS/the9leadMRI.rds")
# the9leadMRI <- readRDS("/path/to/project/GWAS_res/PheWAS/the9leadMRI.rds")

# SNPs with no PheWAS hits — novel loci
the9leadMRI_nullPheWAS <- the9leadMRI[the9leadMRI$num_of_PheWAS == 0, ]
length(unique(the9leadMRI_nullPheWAS$MARKER))   # 19 unique SNPs
length(unique(the9leadMRI_nullPheWAS$LOCUS_ID))  # 15 loci
View(the9leadMRI_nullPheWAS)



### Stratified S-LDSC cell-type enrichment analysis -----------------------------
### -----------------------------------------------------------------------------

# Stratified LDSC tests enrichment of heritability within chromatin annotations.
# Cell-type S-LDSC conditions on 97 baseline annotations to isolate tissue-specific signal:
#   chi2 ~ baseline_1 + ... + baseline_97 + liver_cell_type_annotation
# Uses the 1000G_Phase3_cell_type_groups annotation (10 tissue types including Liver).
# Reference: Finucane et al. (2015) Nature Genetics.
# Baseline annotation: 1000G_Phase3_baselineLD_v2.2 (Gazal et al. 2017), 97 annotations
# covering coding regions, UTRs, enhancers, TFBS, histone marks, etc.

S_LDSC_DIR <- "/path/to/server/project/GWAS_regenie/LDSC_new/S_LDSC_results"

GenomicSEM_names <- readRDS("/path/to/project/GenomicSEM/final_GenomicSEM_names.rds")
lead_MRI_names   <- GenomicSEM_names[-length(GenomicSEM_names)]

S_LDSC_res <- do.call(rbind, lapply(lead_MRI_names, function(nm) {
  f <- file.path(S_LDSC_DIR, paste0(nm, ".results"))
  x <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  cbind(Trait = nm, x[1, , drop = FALSE])
}))

View(S_LDSC_res)

# Liver chromatin enrichment summary
min(S_LDSC_res$Enrichment); max(S_LDSC_res$Enrichment)  # 2.91 to 5.64
S_LDSC_res$Enrichment_p < 0.05 / nrow(S_LDSC_res)  # 7 of 9 pass Bonferroni
p.adjust(S_LDSC_res$Enrichment_p, method = "BH")

# Add liver fat as a comparator
f_liverfat <- "/path/to/server/project/GWAS_regenie/LDSC_new/S_LDSC_results_test/liver_fat.results"
x_liverfat <- read.table(f_liverfat, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, check.names = FALSE)
S_LDSC_res <- rbind(S_LDSC_res, cbind(Trait = "liver_fat", x_liverfat[1, , drop = FALSE]))
saveRDS(S_LDSC_res, "/path/to/project/GWAS_res/S_LDSC_res.rds")



### PDFF-lead variant effects on all 59 MRI features (scatter plot) ------------
### -----------------------------------------------------------------------------

## Compare effects of two liver fat GWAS lead variants on each MRI feature.
## Blue reference line shows the liver fat beta ratio, i.e. points on the line
## are consistent with both SNPs acting purely through the liver fat pathway.
library(data.table); library(ggplot2)

plot_two_snp_effects <- function(rsid_x, rsid_y,
                                  dir0    = '/path/to/project/targetMR/liverMRI_GWAS',
                                  lf_file = '/path/to/project/targetMR/liver_fat_leadGWAS.tsv') {
  fs <- list.files(dir0, pattern = '_over_PDFFsnps\\.tsv$', full.names = TRUE)
  dt <- rbindlist(lapply(fs, function(f) {
    x <- fread(f)[TEST == "ADD" & SNP %in% c(rsid_x, rsid_y), .(SNP, BETA, SE)]
    x[, MRI := sub('_over_PDFFsnps\\.tsv$', '', basename(f))]
    x
  }))
  plot_dt <- dcast(dt, MRI ~ SNP, value.var = c("BETA", "SE"))
  setnames(plot_dt,
           c(paste0("BETA_", rsid_x), paste0("BETA_", rsid_y),
             paste0("SE_",   rsid_x), paste0("SE_",   rsid_y)),
           c("b_x", "b_y", "se_x", "se_y"), skip_absent = TRUE)

  lf       <- fread(lf_file)[variant_id %in% c(rsid_x, rsid_y), .(variant_id, beta)]
  bx_lf    <- lf[match(rsid_x, variant_id), beta]
  by_lf    <- lf[match(rsid_y, variant_id), beta]
  slope_lf <- by_lf / bx_lf  # liver fat effect ratio used as reference slope

  p <- ggplot(plot_dt, aes(b_x, b_y)) +
    geom_abline(intercept = 0, slope = slope_lf, color = "blue", linewidth = 0.8) +
    geom_errorbar(aes(ymin = b_y - 1.96*se_y, ymax = b_y + 1.96*se_y), width = 0) +
    geom_errorbarh(aes(xmin = b_x - 1.96*se_x, xmax = b_x + 1.96*se_x), height = 0) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = paste0(rsid_x, " association with MRI"),
         y = paste0(rsid_y, " association with MRI")) +
    theme_classic()

  list(data = plot_dt, liver_fat_slope = slope_lf, plot = p)
}


res1 <- plot_two_snp_effects("rs738408", "rs58542926")  # PNPLA3 vs TM6SF2
res1$plot

res2 <- plot_two_snp_effects("rs738408", "rs429358")    # PNPLA3 vs APOE
res2$plot

library(patchwork)
res1$plot + res2$plot
