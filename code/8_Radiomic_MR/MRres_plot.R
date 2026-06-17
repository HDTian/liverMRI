# =============================================================================
# Author: Haodong Tian
# Description: Visualization of MR results — generates forest plots and summary
#              figures for UVMR/MVMR associations.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================


# Given any liver MRI feature and outcome, return the MR fitting results and scatter plot.

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

library(biomaRt)  # for mapping rsIDs to nearest gene names
library(ggrepel)  # for non-overlapping text labels in ggplot


OUTnames <- list.files(file.path(pathname, "new_UVMR"), pattern = "\\.tsv$", full.names = FALSE)
OUTnames <- OUTnames[OUTnames != "BMI.tsv"]
OUTnames <- sub("\\.tsv$", "", OUTnames)
OUTnames  # 12 outcomes of interest

retained_names <- readRDS("/path/to/project/phenotype data/features with other images/retained_names.rds")


folder <- paste0(pathname, "/new_UVMR/retained_UVMR_summary_1e6_new/")
feature_paths <- list.files(folder, pattern = "\\.tsv$", full.names = TRUE)
all_names <- c()  # all liver feature names
for (kk in 1:length(feature_paths)) {
  feature_path <- feature_paths[kk]
  feature_name <- sub("\\.tsv$", "", basename(feature_path))
  all_names    <- c(all_names, feature_name)
}
all_names  # the 59 retained features


do_Steiger <- TRUE


### Load shared data objects
## Individual-level phenotype data (rank-normal transformed and residualized)
liver_data_all <- read.table(
  "/path/to/project/phenotype data/features with other images/rint_and_residualizing.txt",
  header = TRUE, sep = "", stringsAsFactors = FALSE)
dim(liver_data_all)  # 37791 x 76 = ID + 16 covariates + 59 features

# BMI GWAS summary data (used for variant removal and MVMR adjustment)
BMI_for_UVMR <- read.table(paste0(pathname, "/new_UVMR/", 'BMI', '.tsv'),
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)


### Function 1: map rsID list to nearest gene name within ±0.1 Mb
get_closest_gene <- function(rsID_vector) {

  ## 1) rsID -> chromosomal position (GRCh38) via Ensembl BioMart
  m_snp <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp", mirror = "useast")
  v <- as.data.table(getBM(
    attributes = c("refsnp_id", "chr_name", "chrom_start"),
    filters    = "snp_filter", values = rsID_vector, mart = m_snp))
  setnames(v, c("refsnp_id", "chr_name", "chrom_start"), c("rsID", "chr", "pos"))
  v <- v[chr %in% c(as.character(1:22), "X", "Y"), unique(.SD), by = rsID]

  ## 2) For each SNP position, retrieve genes within ±0.1 Mb window
  win <- 1e5  # window size in bp
  v[, `:=`(start = pmax(1L, pos - win), end = pos + win)]
  v[, region := paste0(chr, ":", start, ":", end)]
  m_gene <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  g <- as.data.table(getBM(
    attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
    filters = "chromosomal_region", values = v$region, mart = m_gene))
  setnames(g, c("chromosome_name", "start_position", "end_position", "external_gene_name"),
           c("chr", "gstart", "gend", "gene"))

  ## 3) For each SNP, find the nearest gene
  res <- rbindlist(lapply(seq_len(nrow(v)), function(i) {
    chr <- v$chr[i]; pos <- v$pos[i]; rs <- v$rsID[i]
    g   <- as.data.frame(g)
    gg  <- g[g$chr == chr, ]
    gg  <- gg[(abs(gg$gstart - pos) < win) | (abs(gg$gend - pos) < win), ]
    gg  <- gg[!is.na(gg$gene) & gg$gene != "", ]
    if (!nrow(gg)) return(data.table(rsID = rs, nearest_gene = NA_character_, distance_bp = NA_real_))
    # Distance = 0 if SNP falls within the gene body
    d <- ifelse(pos < gg$gstart, gg$gstart - pos,
                ifelse(pos > gg$gend, pos - gg$gend, 0))
    min_dist <- min(d)
    if (min_dist == 0) {
      # SNP is within one or more genes (e.g. H2BC4/HFE); concatenate all
      target_genes <- unique(gg$gene[d == 0])
      final_gene   <- paste(target_genes, collapse = "/")
    } else {
      target_genes <- unique(gg$gene[d == min_dist])
      final_gene   <- paste(target_genes, collapse = "/")
    }
    data.table(rsID = rs, nearest_gene = final_gene, distance_bp = min_dist)
  }))
  return(res)
}
# Example usage:
rsID_vector <- c("rs1800562", "rs58542926", "rs429358", "rs738409")
gene_res <- get_closest_gene(rsID_vector)  # distance_bp = 0 means SNP is inside the gene body
gene_res


### Function 2: MR fitting for a specified liver MRI feature and outcome -> results + scatter plot
get_MR_res_plot <- function(liverMRIname,
                            outcome_name,
                            outcome_show = NA,
                            pthreshold   = 1 * 10^(-6),
                            text_style   = 'rsID'  # or 'gene'
                            ) {
  MR_RES <- list()

  if (!liverMRIname %in% retained_names) stop('No such liver MRI feature name')
  if (!outcome_name %in% OUTnames)       stop('No such outcome name')

  kk           <- match(liverMRIname, all_names)
  feature_path <- feature_paths[kk]
  feature_name <- sub("\\.tsv$", "", basename(feature_path))

  # Load outcome GWAS summary data
  OUTname         <- outcome_name
  OUT_UVMRsummary <- read.table(paste0(pathname, "/new_UVMR/", OUTname, ".tsv"),
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Load feature-specific IV summary data and apply p-value threshold
  dat <- read.table(feature_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dat <- dat[dat$P < pthreshold, ]

  dat <- dat[!is.na(match(dat$SNP, OUT_UVMRsummary$rsID)), ]  # retain SNPs in outcome GWAS
  dat <- dat[!is.na(match(dat$SNP, BMI_for_UVMR$rsID)), ]    # retain SNPs in BMI GWAS


  ### Steiger filtering
  clean_GWAS_data_path <- list.files(
    "/path/to/project/LDSC/clean_GWAS",
    pattern = paste0("^", OUTname, "\\.tsv\\.gz$"), full.names = TRUE)
  N2 <- median(fread(clean_GWAS_data_path)$effective_sample_size, na.rm = TRUE)

  # Outcome-side squared correlation per IV
  local_OUT_UVMRsummary <- OUT_UVMRsummary[match(dat$SNP, OUT_UVMRsummary$rsID), ]
  ref_effect_allele <- dat$ALLELE1
  by_vector   <- local_OUT_UVMRsummary$beta *
    (-1 + 2 * (local_OUT_UVMRsummary$trueEA == ref_effect_allele))
  byse_vector <- local_OUT_UVMRsummary$se
  t_scores    <- by_vector / byse_vector
  abs_r2      <- sqrt(t_scores^2 / (t_scores^2 + N2 - 2))

  # Exposure-side squared correlation per IV
  N1       <- nrow(liver_data_all)  # radiomics feature GWAS sample size (N = 37791)
  t_scores <- dat$BETA / dat$SE
  abs_r1   <- sqrt(t_scores^2 / (t_scores^2 + N1 - 2))

  # Fisher's z transformation for the one-sided Steiger test
  z1      <- 1/2 * log((1 + abs_r1) / (1 - abs_r1))
  z2      <- 1/2 * log((1 + abs_r2) / (1 - abs_r2))
  final_z <- (z1 - z2) / sqrt(1 / (N1 - 3) + 1 / (N2 - 3))
  # One-sided test with Bonferroni correction; retain SNPs consistent with exposure-first direction
  one_side_pvalues <- pnorm(final_z, 0, 1)
  p_threhold <- 0.05 / nrow(dat)
  dat <- dat[one_side_pvalues >= p_threhold, ]

  ###### MVMR with BMI adjustment -------------------------------------------------------
  MR_RES$MVMR_SNP_number <- nrow(dat)
  local_OUT_UVMRsummary <- OUT_UVMRsummary[match(dat$SNP, OUT_UVMRsummary$rsID), ]
  ref_effect_allele <- dat$ALLELE1
  by_vector   <- local_OUT_UVMRsummary$beta *
    (-1 + 2 * (local_OUT_UVMRsummary$trueEA == ref_effect_allele))
  byse_vector <- local_OUT_UVMRsummary$se

  BMI_data_matched <- BMI_for_UVMR[match(dat$SNP, BMI_for_UVMR$rsID), ]
  BMI_harmonized   <- BMI_data_matched$beta *
    (-1 + 2 * (BMI_data_matched$trueEA == dat$ALLELE1))

  MVMRres <- mr_mvivw(mr_mvinput(
    bx   = cbind(dat$BETA, BMI_harmonized),
    bxse = cbind(dat$SE, BMI_data_matched$se),
    by   = by_vector, byse = byse_vector))

  MR_RES$MVMR_BMI_est <- MVMRres@Estimate[1]  # estimate for the liver feature (index 1)
  MR_RES$MVMR_BMI_p   <- MVMRres@Pvalue[1]


  ####### UVMR with BMI-variant removal -------------------------------------------------
  BMI_pvalues <- (BMI_for_UVMR$pvalue)[match(dat$SNP, BMI_for_UVMR$rsID)]
  dat <- dat[BMI_pvalues > 5 * 10^(-4), ]  # key threshold for BMI-associated variant removal

  MR_RES$UVMR_SNP_number <- nrow(dat)

  # Outcome harmonization
  local_OUT_UVMRsummary <- OUT_UVMRsummary[match(dat$SNP, OUT_UVMRsummary$rsID), ]
  ref_effect_allele <- dat$ALLELE1
  by_vector   <- local_OUT_UVMRsummary$beta *
    (-1 + 2 * (local_OUT_UVMRsummary$trueEA == ref_effect_allele))
  byse_vector <- local_OUT_UVMRsummary$se

  # Primary UVMR: IVW (random-effects)
  MRres <- mr_ivw(mr_input(bx = dat$BETA, bxse = dat$SE,
                            by = by_vector, byse = byse_vector))
  MR_RES$IVW_est    <- MRres@Estimate
  MR_RES$IVW_pvalue <- MRres@Pvalue
  MR_RES$Q_pvalue   <- MRres@Heter.Stat[2]  # Cochran Q heterogeneity p-value
  MR_RES$Fvalue     <- MRres@Fstat           # F-statistic for weak instrument check

  dat$snp      <- paste0('snp_', 1:nrow(dat))
  dat$proxysnp <- local_OUT_UVMRsummary$proxyID[match(dat$SNP, local_OUT_UVMRsummary$rsID)]

  ### Annotate IVs with nearest gene name
  gene_res <- get_closest_gene(dat$SNP)
  dat$closest_gene <- gene_res$nearest_gene[match(dat$SNP, gene_res$rsID)]
  dat$closest_gene[is.na(dat$closest_gene)] <- ""  # empty string = no label shown
  MR_RES$dat <- dat

  ### MR scatter plot
  if (text_style == 'rsID') {
    plot_dat <- data.table(bx = dat$BETA, bxse = dat$SE,
                           by = by_vector, byse = byse_vector, label = dat$SNP)
  }
  if (text_style == 'gene') {
    plot_dat <- data.table(bx = dat$BETA, bxse = dat$SE,
                           by = by_vector, byse = byse_vector, label = dat$closest_gene)
  }
  # Orientate: flip SNPs with negative exposure effect so all bx values are positive
  flip_all <- plot_dat$bx < 0
  plot_dat[flip_all, `:=`(bx = -bx, by = -by)]

  xticks <- pretty(plot_dat$bx, n = 5)
  yticks <- pretty(plot_dat$by, n = 5)

  the_outcome_name <- if (is.na(outcome_show)) outcome_name else outcome_show

  # Scatter plot with IVW slope, error bars, and SNP/gene labels
  MR_scatter_plot <- ggplot(plot_dat, aes(bx, by)) +
    geom_point(size = 2, col = 'blue') +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    geom_errorbar(aes(ymin = by - 1.96 * byse, ymax = by + 1.96 * byse), width = 0, col = 'blue') +
    geom_errorbarh(aes(xmin = bx - 1.96 * bxse, xmax = bx + 1.96 * bxse), width = 0, col = 'blue') +
    geom_abline(intercept = 0, slope = MRres@Estimate, linetype = "dashed",
                linewidth = 0.8, col = 'blue') +
    geom_text_repel(aes(label = label), fontface = "italic",
                    direction = 'y', hjust = 0, vjust = 0,
                    nudge_x = 0.005 * diff(range(plot_dat$bx)),
                    nudge_y = 0.005 * diff(range(plot_dat$by)),
                    min.segment.length = Inf,
                    size = 3, box.padding = 0.3, point.padding = 0.2, max.overlaps = 2) +
    labs(x = liverMRIname, y = the_outcome_name) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
    theme_classic() +
    theme(panel.border     = element_blank(),
          axis.line.x      = element_line(color = "black", linewidth = 1),
          axis.line.y      = element_line(color = "black", linewidth = 1),
          axis.line        = element_blank(),
          axis.ticks.length = unit(0.15, "cm"))

  MR_RES$MR_scatter_plot <- MR_scatter_plot


  # Robust MR sensitivity analyses
  MRmedian <- mr_median(mr_input(bx = dat$BETA, bxse = dat$SE,
                                  by = by_vector, byse = byse_vector),
                         weighting = 'simple')
  MR_RES$weighted_median <- MRmedian@Pvalue

  MRmode <- mr_mbe(mr_input(bx = dat$BETA, bxse = dat$SE,
                             by = by_vector, byse = byse_vector))
  MR_RES$weighted_mode <- MRmode@Pvalue

  MRRAPS <- mr.raps(data.frame(beta.exposure = dat$BETA, beta.outcome = by_vector,
                                se.exposure = dat$SE, se.outcome = byse_vector),
                    diagnostics = FALSE, over.dispersion = TRUE)
  MR_RES$MR_RAPS <- 2 * pnorm(-abs(MRRAPS$beta.hat / MRRAPS$beta.se))

  return(MR_RES)
}


### Example calls — selected MRI-outcome pairs of interest

text_style_used <- 'gene'
text_style_used <- 'rsID'

### Final figure layout (3x3 panel)
library(patchwork)

# Row 1: TG | Lp(a) | T2D
MRres1$MR_scatter_plot | MRres4$MR_scatter_plot | MRres91$MR_scatter_plot

# Row 2: CAD | LDL-C | Liver fat
MRres5$MR_scatter_plot | MRres52$MR_scatter_plot | MRres7$MR_scatter_plot

# Row 3: Liver cT1 | Cirrhosis | Liver cancer
MRres100$MR_scatter_plot | MRres12$MR_scatter_plot | MRres11$MR_scatter_plot

MR_scatterplot <- (MRres1$MR_scatter_plot  | MRres4$MR_scatter_plot  | MRres91$MR_scatter_plot) /
  (MRres5$MR_scatter_plot   | MRres52$MR_scatter_plot | MRres7$MR_scatter_plot) /
  (MRres100$MR_scatter_plot | MRres12$MR_scatter_plot | MRres11$MR_scatter_plot)
MR_scatterplot & theme(plot.margin = margin(20, 20, 20, 20))
