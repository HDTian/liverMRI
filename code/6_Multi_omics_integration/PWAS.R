# =============================================================================
# Author: Haodong Tian
# Description: Proteome-wide association study mapping liver MRI features to
#              plasma protein levels using summary-level MR (two-sample MR with
#              cis-pQTLs).
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

### PWAS: two-sample summary-based MR (IVW) using cis-pQTLs from UKB-PPP ST16

### Design rationale:
# Individual-level two-sample 2SLS is used here because:
# - The MRI and protein cohorts do not fully overlap, making one-sample 2SLS underpowered
# - Using [variants + protein + covariates] for first-stage regression and
#   [variants + MRI + covariates] for second-stage regression implements two-sample 2SLS
# - Summary-based IVW is also valid; both approaches yield equivalent estimates
#   (verified by simulation: est/se/p identical to MRI ~ predicted_protein_level + covariates)
# - Only cis-pQTLs are used as instruments (as in Mendelian randomization); trans-pQTLs
#   are excluded to improve instrument validity


library(data.table)
library(readxl)
library(dplyr)
library(ggplot2)



### [ID x MRI + covariates] data: RINT-transformed MRI phenotype data ============================================================
radiomics_wide <- fread("/path/to/server/project/GWAS_regenie/my_liver_UKB_more_covar_new.txt")  # ~2 mins
dim(radiomics_wide)  # 37791 x 231: the GWAS analysis cohort

# Retain IID + 59 MRI features + covariates
retained_names <- readRDS("/path/to/project/phenotype data/features with other images/retained_names.rds")

covariates <- c("sex", "age", "age2", "sex_age", "sex_age2", paste0("PC", 1:10), "genotyping_array")  # sex=1 means female
keep_cols  <- c("IID", retained_names, covariates)
keep_cols  <- intersect(keep_cols, colnames(radiomics_wide))
radiomics_sub <- radiomics_wide[, ..keep_cols]

dim(radiomics_sub)           # 37791 x 76: 1 ID + 59 MRI features + 16 covariates
sum(!complete.cases(radiomics_sub))  # no missing values

## Restrict to unrelated EUR individuals (required for all individual-level analyses outside GWAS)
UKB_unrelated_EUR <- fread("/path/to/server/project/general/data/UKB_unrelated_EUR.txt")
dim(UKB_unrelated_EUR)  # 404585 x 1

radiomics_sub_unrelated <- radiomics_sub[IID %in% UKB_unrelated_EUR[[1]]]
dim(radiomics_sub_unrelated)  # 33122 x 76


## Apply rank-based inverse normal transformation (RINT) to the 59 MRI features
RINT <- function(x) { r <- rank(x, na.last = "keep", ties.method = "average"); qnorm((r - 0.5) / sum(!is.na(x))) }
radiomics_sub_rint <- copy(radiomics_sub_unrelated)
radiomics_sub_rint[, (retained_names) := lapply(.SD, RINT), .SDcols = retained_names]
apply(radiomics_sub_rint[, ..retained_names], 2, sd)  # all approximately 1 after RINT

dim(radiomics_sub_rint)  # 33122 x 76




### [ID x protein + covariates] data: UKB-PPP individual-level protein data ==============================
proteins_wide <- fread("/path/to/server/UK_Biobank/baskets/ukb_basket_olink/olink_data3k_pivot.tsv.gz")  # ~2 mins
dim(proteins_wide)  # 55313 x 2925: one row per participant (eid, ins_index, protein columns)
proteins_wide <- proteins_wide[which(proteins_wide$ins_index == 0), ]  # retain baseline (instance 0)
dim(proteins_wide)  # 53018 x 2925

# Decode protein column names from numeric field codes using the coding143 lookup table
linker <- fread("/path/to/server/UK_Biobank/baskets/ukb_basket_olink/coding143.tsv")
linker$meaning    <- sapply(strsplit(linker$meaning, ";"), `[`, 1)
new_colnames      <- setNames(linker$meaning, linker$coding)
num_colnames      <- names(proteins_wide)[3:ncol(proteins_wide)]
num_colnames      <- as.numeric(num_colnames)
matched_colnames  <- new_colnames[match(num_colnames, as.numeric(names(new_colnames)))]
is_na             <- is.na(matched_colnames); matched_colnames[is_na] <- as.character(num_colnames)[is_na]
names(proteins_wide)[3:ncol(proteins_wide)] <- matched_colnames

dim(proteins_wide)  # 53018 x 2925: columns are now protein gene names
# Note: protein levels are on the original (non-RINT) scale; UKB-PPP GWAS used RINT-transformed values

## Restrict to unrelated EUR individuals
proteins_wide_unrelated <- proteins_wide[eid %in% UKB_unrelated_EUR[[1]]]
dim(proteins_wide_unrelated)  # 42803 x 2925

## Append covariates to the protein data
# Covariates: sex, age, age2, sex_age, sex_age2, PC1-10, genotyping_array
complete_covariates_data <- fread("/path/to/project/PWAS_Coloc/ukbb_covars_first100col.txt",
           select = c("id","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","Sex","age","genotyping_array"))
dim(complete_covariates_data)  # 502629 x 14
complete_covariates_data$sex      <- as.numeric(complete_covariates_data$Sex == 'Female')
complete_covariates_data$sex_age  <- complete_covariates_data$sex * complete_covariates_data$age
complete_covariates_data$age2     <- complete_covariates_data$age^2
complete_covariates_data$sex_age2 <- complete_covariates_data$sex * complete_covariates_data$age2

covar_cols           <- c("sex","age","age2","sex_age","sex_age2",paste0("PC",1:10),"genotyping_array")
protein_covariates_dat <- proteins_wide_unrelated %>%
  left_join(complete_covariates_data[, c("id", covar_cols), with=FALSE], by=c("eid"="id"))
protein_covariates_dat <- protein_covariates_dat %>% dplyr::select(eid, all_of(covar_cols), everything())
dim(protein_covariates_dat)  # 42803 x 2941




### [protein x info] data: UKB-PPP fine-mapping results (ST16) ================================================
ST16 <- read_excel("/path/to/project/PWAS_Coloc/UKB_PPP_supp.xlsx",
  sheet = "ST16", skip = 4)
dim(ST16)  # 29420 x 16
names(ST16)[1] <- 'Protein'
ST16 <- ST16[ST16$`Cis/trans` == 'cis', ]  # only cis-pQTLs as instruments (Mendelian randomization principle)
dim(ST16)  # 10750 x 16

ST16$Protein <- sub(":.*", "", ST16$Protein)  # keep only readable protein name (strip isoform suffix)

# Check uniqueness of cis-region per protein
setDT(ST16)
check_region <- ST16[, .(n_rows = .N, n_unique_region = uniqueN(`Test region hg19`), regions = paste(unique(`Test region hg19`), collapse = " | ")), by = Protein]
check_region[n_unique_region > 1]
# Note: CXCL8, LMOD1, TNF have multiple overlapping cis-regions (multi-signal loci);
# EBI3_IL27 maps to two chromosomes as it encodes a heterodimeric complex — handle with care in coloc.

# Verify that all ST16 proteins are present in the server protein dataset
ST16[, in_server := Protein %in% linker$meaning]
table(ST16$in_server)  # all TRUE

length(unique(ST16$Protein))  # 1954 unique proteins; 10750 cis-pQTLs

dim(ST16); names(ST16)




### [ID x variant] dosage data: genotype dosage for all causal cis-pQTLs =====================================================================================

# Use CHR:BP:REF:ALT as variant identifiers rather than rsIDs to avoid multi-allelic mismatches
# REF:ALT is for precise variant definition, not to specify the effect allele (irrelevant for 2SLS)

tmp       <- do.call(rbind, strsplit(ST16$`Variant ID`, ":"))
ST16$CHR   <- tmp[,1]
ST16$BP    <- tmp[,2]
ST16$REF   <- tmp[,3]
ST16$ALT   <- tmp[,4]
ST16$CHR_BP <- paste0(tmp[,1], ":", tmp[,2])

length(radiomics_sub_rint$IID)       # 33122: unrelated EUR UKB-MRI individuals
length(protein_covariates_dat$eid)   # 42803: unrelated EUR UKB-PPP individuals

common_IID <- intersect(radiomics_sub_rint$IID, protein_covariates_dat$eid)  # overlap
length(common_IID)  # 4384
all_IID    <- union(radiomics_sub_rint$IID, protein_covariates_dat$eid)       # union for dosage extraction
length(all_IID)     # 71541


## STEP 1: Write sample and variant lists to server [run once]
write.table(data.frame(FID = all_IID, IID = all_IID),
            "/path/to/server/project/PWAS/target_samples.txt",
            row.names=FALSE, col.names=FALSE, quote=FALSE)
# Write per-chromosome variant ID files (UKB BGEN is split by chromosome)
for (my_chr in sort(unique(ST16$CHR))) {
  write.table(
    data.frame(VAR = paste0(ST16$CHR[ST16$CHR == my_chr], ":", ST16$BP[ST16$CHR == my_chr], ":", ST16$REF[ST16$CHR == my_chr], ":", ST16$ALT[ST16$CHR == my_chr])),
    paste0("/path/to/server/project/PWAS/target_variants_chr", my_chr, ".txt"),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

## STEP 2: On server — extract dosage using bgenix + PLINK2
# Install bgenix to build a .bgi index file; this avoids PLINK2 scanning the full chromosome BGEN (prevents OOM)
# Then submit: qsub submit_dosage.sh (~10 mins)


## STEP 3: Inspect dosage output for chromosome 1
geno_chr1 <- fread("/path/to/server/project/PWAS/geno_raw/pwas_chr1.raw")  # ~2 mins
dim(geno_chr1)  # 71278 x 1251: columns include 6 metadata fields + variant columns in rsID_Allele form

# Verify match between ST16 rsIDs (chr1) and dosage column names
length(unique(ST16$rsID[ST16$CHR == '1']))  # 1228 cis-pQTLs on chr1

rsid_ST16_chr1  <- unique(ST16$rsID[ST16$CHR == "1"])
rsid_geno_chr1  <- sub("_[^_]+$", "", colnames(geno_chr1)[-(1:6)])
sum(rsid_ST16_chr1 %in% rsid_geno_chr1)  # 1228: all ST16 chr1 variants captured
setdiff(rsid_ST16_chr1, rsid_geno_chr1)[1:20]




### First-stage 2SLS: regress RINT-protein ~ cis-pQTLs + covariates to obtain IV weights =============================

RINT      <- function(x) { r <- rank(x, na.last = "keep", ties.method = "average"); qnorm((r - 0.5) / sum(!is.na(x))) }
covariates <- c("sex","age","age2","sex_age","sex_age2",paste0("PC",1:10),"genotyping_array")

ST16_use <- copy(ST16); ST16_use[, CHR := as.character(CHR)]

prot_dt  <- copy(protein_covariates_dat); prot_dt[, IID := as.character(eid)]
base_dt  <- prot_dt[, c("IID", covariates, unique(ST16_use$Protein)), with = FALSE]
dim(base_dt)  # 42803 x 1971

## Exclude MRI individuals to enforce strict two-sample 2SLS
base_dt  <- base_dt[!base_dt$IID %in% radiomics_sub_rint$IID, ]
dim(base_dt)  # 38419 x 1971: unrelated EUR, MRI-excluded individuals


iv_stats_list <- list()  # first-stage statistics per protein
weight_list   <- list()  # per-variant weights per protein

for (chr_now in 1:22) {
  cat(paste0("==================================================================\n"))
  cat(paste0("Processing chr ", chr_now, "\n"))


  ### Load dosage data for this chromosome
  geno_file <- paste0("/path/to/server/project/PWAS/geno_raw/pwas_chr", chr_now, ".raw")
  if (!file.exists(geno_file)) next
  geno_dt <- fread(geno_file)
  geno_dt[, IID := as.character(IID)]
  old_snp_names <- names(geno_dt)[-(1:6)]  # columns 1-6 are FID/IID/PAT/MAT/SEX/PHENOTYPE
  new_snp_names <- sub("_[^_]+$", "", old_snp_names)  # strip allele suffix (e.g. rs28615823_C -> rs28615823)

  # Deduplicate rsIDs: retain the variant with the highest MAF within each duplicate group
  maf_dt <- data.table(old_name = old_snp_names, new_name = new_snp_names)
  maf_dt[, dosage_mean := sapply(old_snp_names, function(x) mean(geno_dt[[x]], na.rm=TRUE))]
  maf_dt[, MAF := ifelse(dosage_mean/2 > 0.5, 1 - dosage_mean/2, dosage_mean/2)]
  dup_snps <- maf_dt[duplicated(new_name) | duplicated(new_name, fromLast=TRUE)]
  if (nrow(dup_snps) > 0) { cat(paste0("Duplicate rsIDs found: ", length(unique(dup_snps$new_name)), " (deduplicating)\n")) }
  maf_dt  <- maf_dt[order(-MAF)][!duplicated(new_name)]
  geno_dt <- geno_dt[, c("IID", maf_dt$old_name), with=FALSE]
  setnames(geno_dt, c("IID", maf_dt$new_name))



  ### Merge: (unrelated, MRI-excluded) IID + covariates + proteins + dosage
  dat_chr <- merge(base_dt, geno_dt, by = "IID")
  cat(paste0('Current (ID+covariates+protein+dosage) data nrow: ', nrow(dat_chr), '\n'))

  proteins_chr <- unique(ST16_use[CHR == as.character(chr_now), Protein])
  cat(paste0('Current chr has ', length(proteins_chr), ' proteins\n'))

  for (prot in proteins_chr) {
    cat('-')

    # cis-pQTL rsIDs for this protein
    snps <- unique(ST16_use[CHR == as.character(chr_now) & Protein == prot, rsID])
    snps <- intersect(snps, names(dat_chr))
    if (length(snps) == 0L) next
    if (!(prot %in% names(dat_chr))) next

    dat_fit <- copy(dat_chr[, c("IID", covariates, prot, snps), with = FALSE])
    setnames(dat_fit, prot, "protein_raw")
    dat_fit <- dat_fit[!is.na(protein_raw)]
    dat_fit <- dat_fit[complete.cases(dat_fit[, c(covariates, snps), with = FALSE])]
    if (nrow(dat_fit) < 50L) next

    # RINT-transform the protein for the first-stage regression
    dat_fit[, protein_rint := RINT(protein_raw)]

    # First-stage regression: RINT-protein ~ cis-pQTLs + covariates
    bt       <- function(x) ifelse(make.names(x) == x, x, paste0("`", x, "`"))  # backtick non-standard names (e.g. 1:104093262_CA_C)
    fml_red  <- as.formula(paste("protein_rint ~", paste(bt(covariates), collapse = " + ")))
    fml_full <- as.formula(paste("protein_rint ~", paste(bt(c(covariates, snps)), collapse = " + ")))
    fit_red  <- lm(fml_red, data = dat_fit)
    fit_full <- lm(fml_full, data = dat_fit)
    sm_red   <- summary(fit_red); sm_full <- summary(fit_full); a <- anova(fit_red, fit_full)

    partial_F  <- a$F[2]
    partial_p  <- a$`Pr(>F)`[2]
    df1        <- a$Df[2]
    df2        <- a$Res.Df[2]
    R2_red     <- sm_red$r.squared
    R2_full    <- sm_full$r.squared
    partial_R2 <- (R2_full - R2_red) / (1 - R2_red)  # partial R2 (instrument strength in IV context)

    beta_hat        <- coef(fit_full)[bt(snps)]  # per-variant conditional weights for second-stage PRS
    names(beta_hat) <- snps

    iv_stats_list[[paste0(chr_now, "__", prot)]] <- data.table(CHR = chr_now, Protein = prot, n = nrow(dat_fit), n_snp = length(snps), snps = paste(snps, collapse = "|"), R2_reduced = R2_red, R2_full = R2_full, partial_R2 = partial_R2, partial_F = partial_F, partial_F_df1 = df1, partial_F_df2 = df2, partial_F_p = partial_p)
    weight_list[[paste0(chr_now, "__", prot)]]    <- data.table(CHR = chr_now, Protein = prot, rsID = names(beta_hat), weight = as.numeric(beta_hat))
  }

  cat('Current chr finished\n')

}

protein_iv_stats <- rbindlist(iv_stats_list, fill = TRUE)
dim(protein_iv_stats)  # 1916 x 12
protein_weights  <- rbindlist(weight_list, fill = TRUE)
dim(protein_weights)   # 10545 x 4

fwrite(protein_iv_stats, "/path/to/project/PWAS_Coloc/protein_iv_stats_partialF.txt", sep = "\t")
fwrite(protein_weights,  "/path/to/project/PWAS_Coloc/protein_variant_weights.txt",   sep = "\t")




### Second-stage 2SLS: compute PRS and regress each MRI feature on the predicted protein level ====================================================

covariates <- c("sex","age","age2","sex_age","sex_age2",paste0("PC",1:10),"genotyping_array")
mri_names  <- setdiff(names(radiomics_sub_rint), c("IID", covariates))  # 59 RINT-transformed MRI features
ST16_use   <- copy(ST16); ST16_use[, CHR := as.character(CHR)]

protein_weights     <- fread("/path/to/project/PWAS_Coloc/protein_variant_weights.txt")
protein_weights_use <- copy(protein_weights)
protein_weights_use[, CHR := as.character(CHR)]

base_dt <- copy(radiomics_sub_rint[, c("IID", covariates, mri_names), with = FALSE])
base_dt[, IID := as.character(IID)]
dim(base_dt)  # 33122 x 76

second_stage_list <- list()

for (chr_now in 1:22) {
  cat(paste0("==================================================================\n"))
  cat(paste0("Processing chr ", chr_now, "\n"))

  ### Load dosage data for this chromosome
  geno_file <- paste0("/path/to/server/project/PWAS/geno_raw/pwas_chr", chr_now, ".raw")
  if (!file.exists(geno_file)) next
  geno_dt <- fread(geno_file)
  geno_dt[, IID := as.character(IID)]
  old_snp_names <- names(geno_dt)[-(1:6)]
  new_snp_names <- sub("_[^_]+$", "", old_snp_names)  # strip allele suffix

  # Deduplicate rsIDs (same logic as first stage)
  maf_dt <- data.table(old_name = old_snp_names, new_name = new_snp_names)
  maf_dt[, dosage_mean := sapply(old_snp_names, function(x) mean(geno_dt[[x]], na.rm = TRUE))]
  maf_dt[, MAF := ifelse(dosage_mean/2 > 0.5, 1 - dosage_mean/2, dosage_mean/2)]
  dup_snps <- maf_dt[duplicated(new_name) | duplicated(new_name, fromLast = TRUE)]
  if (nrow(dup_snps) > 0) { cat(paste0("Duplicate rsIDs: ", length(unique(dup_snps$new_name)), " (deduplicating)\n")) }
  maf_dt  <- maf_dt[order(-MAF)][!duplicated(new_name)]
  geno_dt <- geno_dt[, c("IID", maf_dt$old_name), with = FALSE]
  setnames(geno_dt, c("IID", maf_dt$new_name))

  ### Merge: unrelated EUR MRI IDs + covariates + MRI features + dosage
  dat_chr <- merge(base_dt, geno_dt, by = "IID")
  cat(paste0("Current (ID+covariates+MRI+dosage) data nrow: ", nrow(dat_chr), "\n"))

  proteins_chr <- unique(protein_weights_use[CHR == as.character(chr_now), Protein])
  cat(paste0("Current chr has ", length(proteins_chr), " proteins\n"))

  for (prot in proteins_chr) {
    cat("-")

    w_dt <- copy(protein_weights_use[CHR == as.character(chr_now) & Protein == prot, .(rsID, weight)])
    w_dt <- w_dt[rsID %in% names(dat_chr)]
    if (nrow(w_dt) == 0L) next

    # Compute PRS = sum(SNP_j * weight_j) as the predicted protein level
    G   <- as.matrix(dat_chr[, w_dt$rsID, with = FALSE])
    storage.mode(G) <- "numeric"
    prs <- as.numeric(G %*% w_dt$weight)

    dat_fit <- copy(dat_chr[, c("IID", covariates, mri_names), with = FALSE])
    dat_fit[, PRS := prs]
    dat_fit <- dat_fit[complete.cases(dat_fit[, c(covariates, "PRS"), with = FALSE])]
    if (nrow(dat_fit) < 50L) next

    # Second-stage regression: MRI ~ PRS + covariates (for each of the 59 MRI features)
    bt <- function(x) ifelse(make.names(x) == x, x, paste0("`", x, "`"))
    for (mri in mri_names) {
      dat_fit_mri <- dat_fit[!is.na(get(mri))]
      if (nrow(dat_fit_mri) < 50L) next
      fml <- as.formula(paste(bt(mri), "~", paste(c("PRS", covariates), collapse = " + ")))
      fit <- lm(fml, data = dat_fit_mri)
      sm  <- summary(fit)
      second_stage_list[[paste0(chr_now, "__", prot, "__", mri)]] <- data.table(CHR = chr_now, Protein = prot, MRI = mri, n = nrow(dat_fit_mri), n_snp = nrow(w_dt), beta = unname(coef(fit)["PRS"]), se = unname(coef(sm)["PRS","Std. Error"]), t = unname(coef(sm)["PRS","t value"]), p = unname(coef(sm)["PRS","Pr(>|t|)"]))
    }

  }
  cat("Current chr finished\n")

}

second_stage_res <- rbindlist(second_stage_list, fill = TRUE)
dim(second_stage_res)  # 113044 = 59 MRI x 1916 proteins, 9 columns

fwrite(second_stage_res, "/path/to/project/PWAS_Coloc/second_stage_protein_MRI_results.txt", sep = "\t")




### Results visualization ================================================================================

## Protein-MRI association heatmap
# Proteins are classified into three biological axes:
# - Lipid/metabolic: APOE, NCAN, ADH4, FABP2, GSTA1/3, BCAT1, SERPINA1, SERPINA6, TREH, CELA2A
# - Fibrosis/ECM:    ACTA2, INHBB, INHBC, EFEMP1, PDE5A
# - Immune/inflammation: AIF1, FOLR2, CASP9, TNFSF10, DAPP1, ARHGAP25
# If the significant protein set changes, update the category assignments below manually.

library(data.table); library(ggplot2); library(grid); library(ggtext)

x <- as.data.table(second_stage_res)
x[, fdr := p.adjust(p, method = "fdr")]
mri_order <- cluster_df$feature[cluster_df$feature %in% unique(x$MRI)]
sig_prot  <- x[fdr < 0.05, unique(Protein)]

prot_annot <- data.table(Protein = sig_prot)
prot_annot[, cat := fcase(Protein %in% c("APOE","NCAN","ADH4","FABP2","GSTA1","GSTA3","BCAT1","SERPINA1","SERPINA6","TREH","CELA2A"), "Lipid / Metabolic", Protein %in% c("ACTA2","INHBB","INHBC","EFEMP1","PDE5A"), "Fibrosis / ECM", Protein %in% c("AIF1","FOLR2","CASP9","TNFSF10","DAPP1","ARHGAP25"), "Immune / Inflammation", default = "Other")]

prot_pos <- as.data.table(ST16)[Protein %in% sig_prot][, CHR_num := suppressWarnings(as.numeric(CHR))][order(CHR_num, BP), .SD[1], by = Protein]
prot_annot <- merge(prot_annot, prot_pos[, .(Protein, CHR_num, BP)], by = "Protein", all.x = TRUE)

cat_order  <- c("Lipid / Metabolic","Fibrosis / ECM","Immune / Inflammation","Other")
cat_colors <- c("Lipid / Metabolic"="#A65E00","Fibrosis / ECM"="#3E6652","Immune / Inflammation"="#2F4F8F","Other"="#555555")
prot_annot[, cat := factor(cat, levels = cat_order)]
prot_annot <- prot_annot[order(cat, CHR_num, BP)]
prot_order <- prot_annot$Protein

plot_dt <- merge(x[Protein %in% sig_prot & MRI %in% mri_order, .(Protein, MRI, beta, fdr)], prot_annot[, .(Protein, cat)], by = "Protein", all.x = TRUE)
plot_dt[, MRI     := factor(MRI,     levels = mri_order)]
plot_dt[, Protein := factor(Protein, levels = rev(prot_order))]
plot_dt[, cat     := factor(cat,     levels = cat_order)]

cf      <- as.data.table(cluster_df)[feature %in% mri_order][, x := match(feature, mri_order)]
clu_pos <- cf[, .(xmin = min(x)-.5, xmax = max(x)+.5, xmid = mean(range(x))), by = cluster]
vline_x <- clu_pos$xmax[-nrow(clu_pos)]

axis_colors    <- cat_colors[as.character(prot_annot$cat[match(rev(prot_order), prot_annot$Protein)])]
legend_dt      <- data.table(cat = factor(names(cat_colors), levels = cat_order), MRI = factor(mri_order[1], levels = mri_order), Protein = factor(rev(prot_order)[1], levels = rev(prot_order)))
prot_label_map <- setNames(paste0("<span style='color:", axis_colors, ";'><i>", rev(prot_order), "</i></span>"), rev(prot_order))


ggplot(plot_dt, aes(MRI, Protein, fill = beta)) +
  geom_tile(color = "grey92", linewidth = 0.15) +
  geom_point(data = plot_dt[fdr < 0.05], shape = 8, size = 0.8, color = "black") +
  geom_vline(xintercept = vline_x, color = "black", linewidth = 0.4) +
  geom_text(data = clu_pos, aes(x = xmid, y = Inf, label = paste0("C", cluster)), inherit.aes = FALSE, vjust = -0.8, size = 2.8, fontface = "bold") +
  geom_point(data = legend_dt, aes(MRI, Protein, color = cat), inherit.aes = FALSE, alpha = 0, size = 4, shape = 15, show.legend = TRUE) +
  scale_color_manual(values = cat_colors, name = "Protein category") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 15))) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, name = "Association size\n(beta)") +
  scale_x_discrete(drop = FALSE, na.translate = FALSE) +
  scale_y_discrete(drop = FALSE, na.translate = FALSE, labels = prot_label_map) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
        axis.text.y = ggtext::element_markdown(size = 6),
        panel.grid = element_blank(),
        legend.key.height = unit(0.4, "cm"),
        plot.margin = margin(30, 40, 10, 60))

dim(plot_dt)  # 1652 = 59 MRI x 28 proteins, 5 columns
saveRDS(plot_dt, file = "/path/to/project/PWAS_Coloc/PWAS_heatmap_dt.rds")



# Subset to FDR-significant protein-MRI pairs
plot_dt_sub <- plot_dt[plot_dt$fdr < 0.05, ]
nrow(plot_dt_sub)                    # 93 protein-MRI pairs
length(unique(plot_dt_sub$Protein))  # 28 unique proteins
length(unique(plot_dt_sub$MRI))      # 37 unique MRI features
names(plot_dt_sub)



### Summary bar plots: top proteins and MRI features by significant pair count
library(data.table); library(ggplot2)
prot_cols   <- c("Lipid / Metabolic"="#C97B00","Fibrosis / ECM"="#4E7F5A","Immune / Inflammation"="#4A6FA5","Other"="#7F7F7F")
simple_theme <- theme_classic() + theme(panel.border = element_blank(), axis.line = element_line(color = "black"), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(t = 10, r = 10, b = 35, l = 10))

# Top proteins by number of significant MRI associations
top_prot <- plot_dt_sub[, .N, by = .(Protein, cat)][order(-N, Protein)][1:10]
top_prot[, Protein := factor(Protein, levels = Protein)]
p1 <- ggplot(top_prot, aes(x = Protein, y = N, fill = cat)) + geom_col(width = 0.7) + scale_fill_manual(values = prot_cols) + labs(x = " ", y = "Count") + simple_theme + theme(legend.position = "none")
p1

# Top MRI features by number of significant protein associations
top_mri <- plot_dt_sub[, .N, by = MRI][order(-N, MRI)][1:10]
top_mri[, MRI := factor(MRI, levels = MRI)]
p2 <- ggplot(top_mri, aes(x = MRI, y = N)) + geom_col(width = 0.7, fill = "grey50") + labs(x = " ", y = "Count") + simple_theme +
  theme(plot.margin = margin(t = 10, r = 10, b = 35, l = 80))
p2
