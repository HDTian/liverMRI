# =============================================================================
# Author: Haodong Tian
# Description: LD score regression (LDSC) analysis via the GenomicSEM R package
#              to compute genetic correlations between liver MRI features and a
#              panel of external GWAS traits.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

### LDSC analysis via R for external traits
### => used to generate the genetic correlation heatmap

## Step 0: Standardise external trait GWAS files into a clean format for munging
## Step 1: Munge all traits/phenotypes from local PheWAS GWAS results into *.sumstats
## Step 2: Bivariate LDSC for liver_features x traits to estimate genetic correlations

pathname <- '/path/to/project'

library(readxl)
library(GenomicSEM)
library(data.table)
library(ggplot2)


### GenomicSEM input requirements -----------------------------------------
## The summary statistics files input into the munge function at a minimum need to contain five pieces of information:
# rsID of the SNP.
# A1 allele column, with A1 indicating the effect allele.
# A2 allele column, with A2 indicating the non-effect allele.
# beta: Either a logistic or continuous regression effect. (Beta) (for z-score direction)
# p-value associated with this effect estimate beta.
# (The package will automatically rename the column based on commonly observed names)


# Linear mixed models are not strictly appropriate for LDSC, heritability, and genetic correlation analyses;
# however, recent LMM-based GWAS methods are generally considered compatible.

### Sample size considerations (important for LDSC heritability and genetic correlation estimates)
# Be careful to use:
# (i)  The real GWAS sample size, not the combined discovery + validation sample size.
# (ii) The effective sample size for case/control data (logistic model), not the naive sample size.
# (iii) SNP-specific sample size when different SNPs have different sample sizes.

# If SNP-specific sample size is not provided, the SE and EAF can be used to approximate N:
#   N_approx = 1 / (2 * MAF * (1 - MAF) * SE^2)
# This yields a SNP-specific effective sample size, which is preferable when N varies across SNPs.

# It is recommended to confirm whether effective sample size columns reflect
# [sum of effective sample sizes] across cohorts,
# [total effective sample size], or alternatively the effective sample size calculated from
# total cases and controls across cohorts.
# Only the former (sum of effective sample sizes) is appropriate for GenomicSEM analyses.

### Effective sample size definition: Neff := 4 * Ncase * Ncontrol / N
# (so that Neff = N under balanced case:control ratio)
# Effective sample size is directly linked to se(beta)^2, which is why it is called "effective".
# The fixed-effect IVW meta-analyzed SE can be used to approximate Neff [with MAF]:
#   se(beta_meta)^2 = 1 / [Sum_k 1/se(beta_k)^2] ~= 4 / [Var(G) * Sum_k Neff_k]
# Here Neff = Sum_k Neff_k is the correct aggregate effective sample size.


trait_info_table <- read_excel( paste0(pathname , '/trait_info.xlsx' ) );setDT(trait_info_table)
trait_info_table <- trait_info_table[!is.na(`file name`)]

feature_paths <- paste0(pathname , '/PheWAS/' ,  trait_info_table$`file name`    )


## Step 0: Standardise GWAS data format (compute effective sample size for better munging) -----
## ---------------------------------------------------------------------------------------------

out_dir <- file.path(pathname, 'LDSC'  ,"clean_GWAS")


# HM3 SNPs — only HapMap3 SNPs are needed for LDSC
hm3 <- fread(paste0(  pathname, '/LDSC/eur_w_ld_chr/w_hm3.snplist' ))
hm3_snps<-hm3$SNP
hm3_snps <- unique(as.character(hm3_snps))
length(hm3_snps) # 1,217,311

for(kkk in 1:nrow(trait_info_table)){
  current_path <- feature_paths[kkk]
  current_trait <- sub("\\.(tsv|txt|csv)(\\.gz)?$", "", trait_info_table$`file name`[kkk], ignore.case = TRUE)
  print( '================================================' )
  print(   paste0('current phenotype:' , current_trait   )  )

  out_file <- file.path(out_dir, paste0(current_trait, ".tsv.gz"))
  if (file.exists(out_file)) {message("output file already exists; Skip this trait");next}

  # Required columns: rsID, EA, nonEA, beta, se, EAF (EAF used to compute effective sample size)
  needed <- unlist(trait_info_table[kkk, .(rsID, EA, nonEA, beta, se, EAF)], use.names = FALSE)
  needed <- trimws(needed)  # actual column names in the GWAS file for this trait

  print('reading the original GWAS...')
  minor_GWAS <- fread(current_path, select = needed )  # ~ 2 mins depending on file size
  setnames(minor_GWAS, old = needed, new = c("rsID", "EA", "nonEA", "beta", "se",'EAF'))
  minor_GWAS[, `:=`(EA = toupper(EA), nonEA = toupper(nonEA))]

  ## Recompute p-value from beta and SE to avoid inflated/adjusted p-values (e.g. from SAIGE) ----
  print('getting the honest pvalue according to beta and se (rather than the misleading pvalue)...')
  minor_GWAS[, p_value :=
               fifelse(is.finite(beta) & is.finite(se) & se > 0,
                       pchisq((beta/se)^2, df = 1, lower.tail = FALSE),
                       NA_real_) ]

  ### Filter to HM3 SNPs only ---------------------------------------------
  print('only considering HM3 SNPs ... ')

  minor_GWAS[, rsID := as.character(rsID)]
  minor_GWAS <- minor_GWAS[ !is.na(rsID) & (rsID %chin% hm3_snps) ]


  ### Add effective sample size column (continuous traits) --------------------------------
  if( trait_info_table$trait_type[kkk] == "continuous"  ){

    if(  trait_info_table$N[kkk] != "N/A"   ){
      print(  'continuous trait; sample size column detected (so directly use it) '  )
      sample_size_column <- unlist(trait_info_table[kkk, .(N)], use.names = FALSE)
      sample_size_vector <- fread(current_path, select = sample_size_column )
      sample_size_vector <- sample_size_vector[[1]]
      stopifnot(    length(sample_size_vector) == nrow(minor_GWAS)   )

      minor_GWAS[, effective_sample_size:= sample_size_vector]
    }else{
      print(  'continuous trait; no sample size column detected, so use external Excel info'  )
      minor_GWAS[, effective_sample_size:= as.numeric(trait_info_table$`sample size`[kkk])   ]
    }
  }

  ### Binary trait: use MAF + SE to derive SNP-specific effective sample size
  # This gives the sum of effective sample sizes, appropriate for GenomicSEM
  # Formula: Neff = 4 / (2 * EAF * (1 - EAF) * SE^2)
  if( trait_info_table$trait_type[kkk] == "binary"  ){
    print('binraty trait; use se and MAF to approximate the per-SNP effective sample size [sum of effective sample sizes]'  )

    # Per-SNP sum of effective sample sizes using MAF and SE
    neff<-4/(   2*minor_GWAS[,EAF] *(1-minor_GWAS[,EAF] ) * (minor_GWAS[,se])^2   )

    minor_GWAS[, effective_sample_size:= neff   ]
  }


  ### Filter to HM3 SNPs and save clean GWAS ---------------------------------------------
  print('only considering HM3 SNPs and store the partially clean GWAS ... ')

  minor_GWAS[, rsID := as.character(rsID)]
  minor_GWAS_hm3 <- minor_GWAS[ !is.na(rsID) & (rsID %chin% hm3_snps) ]

  print(paste0('clean GWAS done! dim:',  dim(minor_GWAS_hm3)[1]  ,' ', dim(minor_GWAS_hm3)[2]   ))

  ### Save clean GWAS results
  out_file <- file.path(out_dir, paste0(current_trait,".tsv.gz") )
  fwrite(minor_GWAS_hm3, out_file, sep = "\t", na = "NA", compress = "gzip")  # ~ 20 seconds
}


# Output columns: rsID  EA  nonEA   beta  se    p_value   effective_sample_size




## Step 1: Munge GWAS data -----------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------

for(kkk in 1:nrow(trait_info_table)){
  current_path <- feature_paths[kkk]
  current_trait <- sub("\\.(tsv|txt|csv)(\\.gz)?$", "", trait_info_table$`file name`[kkk], ignore.case = TRUE)
  print( '================================================' )
  print(   paste0('current phenotype:' , current_trait   )  )


  if (file.exists( file.path(pathname, "LDSC/munged_GWAS", paste0(current_trait,".sumstats.gz") ) ) ) {
    print('this phenotype already has munged GWAS results; so skip')
    next
  }

  setwd( file.path(pathname, "LDSC/munged_GWAS")   )
  out_file <- file.path(out_dir, paste0(current_trait,".tsv.gz") )
  munge( files = out_file ,  # vector of file names
         hm3 = paste0(  pathname, '/LDSC/eur_w_ld_chr/w_hm3.snplist' )   ,  # HapMap3 SNP list
         trait.names = current_trait ,
         column.names = list(SNP="rsID",A1="EA", A2="nonEA", effect="beta",P='p_value', N="effective_sample_size"),
         parallel=FALSE )  # ~ 1 min for HM3-filtered GWAS


}

# Output: munged data with columns: SNP  N  Z  A1  A2

# N: The sample sizes associated with the traits.
# Note that for binary traits reflecting a meta-analysis across multiple cohorts,
# this should reflect *the sum of effective sample sizes* across contributing cohorts.
# The sample size column should reflect the SNP-specific, sum of effective sample sizes.
# Effective sample size is defined as: 4 * v * (1 - v) * n, where v = sample prevalence (Ncase/N).
# When inputting the sum of effective sample sizes,
# the sample prevalence should be entered as 0.5 when running ldsc,
# reflecting that effective sample size already corrects for sample ascertainment.
# If SNP-specific sample sizes are provided per row, munge uses them unless the user overrides.




### Step 2: Bivariate LDSC analysis for each phenotype --------------------------------------------
## ------------------------------------------------------------------------------------------------

# Suppress verbose function output from ldsc()
quiet <- function(expr) {
  out <- file(nullfile(), "wt"); msg <- file(nullfile(), "wt")
  on.exit({ try(sink(type="message")); try(close(msg));
    try(sink());              try(close(out)) }, add = TRUE)
  sink(out); sink(msg, type="message")
  suppressWarnings(suppressMessages(force(expr)))
}


### Select liver MRI GWAS results for LDSC

### Non-BMI-adjusted GWAS (200 features) --------------------------------------
liver_features_munged_folder <- '/path/to/server/project/GWAS_regenie/LDSC_new/munged_files'
setwd( file.path(pathname, "LDSC/rg_all")   )
# Note: non-BMI-adjusted GWAS contains 200 features;
# visualization focuses on the 59 retained features.
### ------------------------------------------------------------


### BMI-adjusted GWAS (59 retained features) -----------------------------------
liver_features_munged_folder <- '/path/to/server/project/GWAS_regenie/LDSC_new_BMI/munged_files'
setwd( file.path(pathname, "LDSC/rg_all_BMI")   )
### ------------------------------------------------------------


# Liver feature munged files
liver_feature_files <- list.files(liver_features_munged_folder,pattern = "\\.sumstats\\.gz$",
  full.names = TRUE,recursive = FALSE)
liver_feature_files
liver_feature_names <- sub("\\.sumstats\\.gz$", "", basename(liver_feature_files));liver_feature_names


## Restrict to the 59 retained liver MRI features
retained_names <- readRDS("/path/to/project/phenotype data/features with other images/retained_names.rds")
keep_idx <- liver_feature_names %in% retained_names
liver_feature_files <- liver_feature_files[keep_idx]
liver_feature_names <- liver_feature_names[keep_idx]
liver_feature_files
liver_feature_names


# LD score files path — used to fit LDSC regression
LDSC_path <- file.path(pathname , 'LDSC', 'eur_w_ld_chr')


for(kkk in 1:nrow(trait_info_table)){
  current_path <- feature_paths[kkk]
  current_trait <- sub("\\.(tsv|txt|csv)(\\.gz)?$", "", trait_info_table$`file name`[kkk], ignore.case = TRUE)
  print( '================================================' )
  print(   paste0('current phenotype:' , current_trait   )  )

  if (file.exists( file.path( getwd(),
                              paste0( current_trait,'_vs_'  , liver_feature_names[length(liver_feature_names)],"_ldsc.log")
                              ) ) ) {
    print('this phenotype already has the final liver feature bivariate LDSC log results (ie all LDSC finished), so skip')
    next
  }

  out_file <- file.path(pathname, 'LDSC'  ,"munged_GWAS", paste0(current_trait,".sumstats.gz") )

  # For binary traits: use sample_prev = 0.5 because effective sample size already corrects
  # for ascertainment; see GenomicSEM wiki for rationale.
  if( trait_info_table$trait_type[kkk] == "binary"  ){
    sample_prev_used <- 0.5 }else{
      sample_prev_used <- NA
    }
   print( paste0( 'sample_prev_used: ' , sample_prev_used  )  )

    for(l in 1:length(liver_feature_names)){
      cat(  paste0( l , '-'  ) )
      LDSCout <-quiet( ldsc(
        traits          = c(out_file, liver_feature_files[l] ) ,
        trait.names     = c(current_trait , liver_feature_names[l] ),
        sample.prev     = c(sample_prev_used, NA),   # binary: Neff via 0.5; continuous: NA
        population.prev = c(NA,  NA),   # only needed when reporting liability-scale h2 for binary traits
        ld              = LDSC_path,
        wld             = LDSC_path,
        ldsc.log        = paste0(current_trait, '_vs_',liver_feature_names[l] ),
        stand           = TRUE          # also return genetic correlation matrix (rg and SE)
       )  # ~ 30s
      )
    }
   cat('\n')
}




### LDSC results: original MRI vs BMI-adjusted MRI
rg_noBMI_BMI<-fread(file.path(pathname ,'/LDSC/rg_noBMI_vs_BMI_summary_instance2.tsv'))
rg_noBMI_BMI<-fread(file.path(pathname ,'/LDSC/rg_noBMI_vs_BMI_summary_instance0.tsv'))

summary(rg_noBMI_BMI$rg  )
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.8164  0.9216  0.9694  0.9494  0.9940  1.0013
dim(rg_noBMI_BMI) # 59 7
plot(  rg_noBMI_BMI$h2 , rg_noBMI_BMI$h2_BMI,  pch = 16,
       xlab='SNP-heritability (original)', ylab='SNP-heritability (BMI-adjusted)'   )
abline(0,1,col='blue')




### Collect LDSC results  -----------------------------------------------------------------------------------
### --------------------------------------------------------------------------------------------------------
trait_names <- sub("\\.(tsv|txt|csv)(\\.gz)?$", "", trait_info_table$`file name`, ignore.case = TRUE)
logs_dir <- file.path(pathname, "LDSC/rg_all")             ## LDSC with non-BMI-adjusted liver feature GWAS
logs_dir <- file.path(pathname, "LDSC/rg_all_BMI_instance2")  ## LDSC with BMI-adjusted liver feature GWAS (instance 2)
logs_dir <- file.path(pathname, "LDSC/rg_all_BMI_instance0")  ## LDSC with BMI-adjusted liver feature GWAS (instance 0)

# Initialise result matrices
rg_mat <- matrix(NA_real_, length(liver_feature_names), length(trait_names),
                 dimnames = list(liver_feature_names, trait_names))
se_mat <- rg_mat

# List all log files
logs <- list.files(logs_dir, pattern = "_ldsc\\.log$", full.names = TRUE)
length(logs) == length( trait_names ) * length( liver_feature_names )  # must be TRUE

# Parse genetic correlation from a single log file; returns (name1, name2, rg, se)
parse_log <- function(f){
  x <- readLines(f, warn = FALSE)
  idx <- grep("^Genetic Correlation between ", x)
  if(!length(idx)) return(NULL)
  line <- x[tail(idx, 1L)]
  m <- regexec("^Genetic Correlation between\\s+(.+?)\\s+and\\s+(.+?):\\s*([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\\s*\\(([^)]+)\\)", line)
  r <- regmatches(line, m)[[1]]
  if(length(r) < 5) return(NULL)
  list(a = trimws(r[2]), b = trimws(r[3]),
       rg = as.numeric(r[4]), se = as.numeric(r[5]))
}


# Determine which of the two names in a log is the liver feature and which is the external trait;
# falls back to parsing the filename (*_vs_*_ldsc.log) if the log content is ambiguous.
who_is_who <- function(a, b, f){
  if(a %in% liver_feature_names && b %in% trait_names) return(list(feat=a, trait=b))
  if(b %in% liver_feature_names && a %in% trait_names) return(list(feat=b, trait=a))
  base <- sub("_ldsc\\.log$", "", basename(f))
  parts <- strsplit(base, "_vs_")[[1]]
  if(length(parts) == 2){
    p1 <- parts[1]; p2 <- parts[2]
    if(p1 %in% liver_feature_names) return(list(feat=p1, trait=p2))
    if(p2 %in% liver_feature_names) return(list(feat=p2, trait=p1))
  }
  NULL
}

for(f in logs){
  pr <- parse_log(f)
  if(is.null(pr) || any(!is.finite(c(pr$rg, pr$se)))) next
  id <- who_is_who(pr$a, pr$b, f)
  if(is.null(id)) next
  if(id$feat %in% rownames(rg_mat) && id$trait %in% colnames(rg_mat)){
    rg_mat[id$feat, id$trait] <- pr$rg
    se_mat[id$feat, id$trait] <- pr$se
  }
}

dim(rg_mat); dim(se_mat)  # 59 13




### Visualization ---------------------------------------------------------------------------------
### -----------------------------------------------------------------------------------------------

stopifnot(identical(rownames(rg_mat), rownames(se_mat)), identical(colnames(rg_mat), colnames(se_mat)))

## Read liver feature clustering information (from GenomicSEM.R; used for consistent ordering with MR heatmap)
cluster_df <- readRDS(  paste0( pathname, "/LDSC/Rg_cluster_df.rds") ) ; names(cluster_df)[1]<-'feature'

# 1) Sort features by cluster (preserving original order within each cluster) ----------------
cf <- data.table(feature = rownames(rg_mat), idx = seq_len(nrow(rg_mat)))[cluster_df, on = "feature"]
if (anyNA(cf$cluster)) {
  mx <- suppressWarnings(max(cf$cluster, na.rm = TRUE)); if (!is.finite(mx)) mx <- 0
  cf[is.na(cluster), cluster := mx + 1L]}
setorder(cf, cluster, idx); feature_levels <- cf$feature
# ------------------------------------------------------------------------------------------


## Set trait display order
trait_levels <- c('BMI','TG','HDL_C','TG_HDL_C','Lp_a','LDL_C','CAD',
                  'HbA1c','T2D','liver_fat','liver_iron','liver_cancer','Cirrhosis')
rg_mat <- rg_mat[, colnames(rg_mat) %in% trait_levels, drop = FALSE]
se_mat <- se_mat[, colnames(se_mat) %in% trait_levels, drop = FALSE]
rg_mat <- rg_mat[, trait_levels[trait_levels %in% colnames(rg_mat)], drop = FALSE]
se_mat <- se_mat[, trait_levels[trait_levels %in% colnames(se_mat)], drop = FALSE]
length(colnames(rg_mat)  ) == length( trait_levels )

# 2) Wide to long; compute p-values and significance flags; set tile sizes for the heatmap
dt_rg <- melt(as.data.table(rg_mat, keep.rownames = "feature"), id = "feature", variable.name = "trait", value.name = "rg")
dt_se <- melt(as.data.table(se_mat, keep.rownames = "feature"), id = "feature", variable.name = "trait", value.name = "se")
dt <- merge(dt_rg, dt_se, by = c("feature","trait"), all = TRUE)
data.table::setDT(dt)

dt[, z := rg / se]  # rg z-score
dt[, p := fifelse(is.finite(z), 2 * pnorm(abs(z), lower.tail = FALSE), NA_real_)]

dt[, rg_nom := fifelse(!is.na(p) & p < 0.05, rg, NA_real_)]   # nominally significant -> small tile
dt[, q := p.adjust(p, method = "fdr")]                        # global FDR adjustment across all pairs
dt[, rg_fdr := fifelse(!is.na(q) & q < 0.05, rg, NA_real_)]   # FDR significant -> full tile
alpha_star <- 0.05 / (dim(rg_mat)[1]*dim(rg_mat)[2])  # global Bonferroni threshold
dt[, star   := fifelse(!is.na(p) & p < alpha_star, "*", "")]   # Bonferroni significant -> star

dt[, z_abs := abs(rg / se)]  # |z| for tile sizing
z_max <- min( abs(dt$z_abs)[  !is.na(dt$rg_fdr) ] )  # minimum |z| among FDR-significant pairs (maps to full tile)
dt[, tile_size := pmin(z_abs / z_max, 1)]
z_min <- min(dt$tile_size[ (dt$tile_size!=1)&( !is.na(dt$rg_nom) ) ])
dt[, tile_size := 0.5 + 0.3 * (tile_size-z_min)/(1-z_min)]  # map to [0.5, 0.8]; 0.5 = minimum visible tile
dt$tile_size[dt$tile_size<0.5] <- 0  # non-significant tiles get zero width (invisible)


# Apply factor levels to control display order
dt[, feature := factor(feature, levels = feature_levels)]
dt[, trait   := factor(trait,   levels = trait_levels)]

# 3) Compute cluster bounding boxes and cluster labels for annotation
pos <- setNames(seq_along(feature_levels), feature_levels)
cb <- rbindlist(lapply(split(cf$feature, cf$cluster, drop = TRUE), function(v){
  r <- range(pos[v])
  data.table(xmin = 0.5, xmax = length(trait_levels) + 0.5, ymin = r[1] - 0.5, ymax = r[2] + 0.5)
}), fill = TRUE)

lab_df <- rbindlist(lapply(split(cf$feature, cf$cluster, drop = TRUE), function(v){
  r <- range(pos[v])
  data.table(cluster = unique(cf[feature %in% v, cluster]), y = mean(r))
}), fill = TRUE)
if (nrow(lab_df)) {
  lab_df[, `:=`(x = length(trait_levels) + 0.6, label = paste0("C", cluster),
                y_fac = factor(feature_levels[round(y)], levels = feature_levels))]
}

dt[, rg_show := fifelse(!is.na(p) & p <= alpha_star, rg, NA_real_)]
dt[, star_col := ifelse(star == "", NA, ifelse(abs(rg_show) >= 0.4, "black", "black"))]


# Rename traits for display
dt$trait <- as.character(dt$trait)
dt$trait[ dt$trait == "liver_iron"] <- "liver_cT1"
dt$trait[ dt$trait=='TG_HDL_C'  ] <- 'TG:HDL_C'
trait_levels <- c('BMI','TG','HDL_C','TG:HDL_C','Lp_a','LDL_C','CAD',
                  'HbA1c','T2D','liver_fat','liver_cT1','liver_cancer','Cirrhosis')
dt$trait <- factor(dt$trait,level=trait_levels)

# Save rg estimate table (rg (SE) format) for supplementary tables
dt2 <- data.table::copy(dt)[, val := sprintf("%.3f (%.3f)", rg, se)]
rg_table  <- dcast(dt2, feature ~ trait, value.var = "val")
saveRDS(rg_table, file = "/path/to/project/LDSC/rg_BMI_original.rds")
saveRDS(rg_table, file = "/path/to/project/LDSC/rg_BMI_instance0.rds")
saveRDS(rg_table, file = "/path/to/project/LDSC/rg_BMI_instance2.rds")



p_rg <- ggplot(dt, aes(trait, feature)) +
  # x = trait (columns), y = feature (rows); ordered by factor levels
  geom_tile(fill = "white", color = "grey90", linewidth = 0.2) +
  # Background grid: all tiles filled white with light grey borders
  geom_tile(aes(fill = rg_fdr), width = 0.98, height = 0.98, color = NA) +
  # FDR-significant tiles: near-full size colored tile
  geom_tile(aes(fill = rg_nom, width = tile_size, height = tile_size), color = NA) +
  # Nominally significant tiles: smaller tile sized proportionally to |z|
  geom_tile(fill = NA, color = "grey70", linewidth = 0.2) +
  # Overlay grid lines
  geom_text(
    aes(label = star, color = star_col),
    size = 3.6, vjust = 0.7, show.legend = FALSE
  ) +
  # Star annotation for Bonferroni-significant pairs
  { if (nrow(cb)) geom_rect(
    data = cb,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, color = "grey40", linewidth = 0.4
  ) } +
  # Cluster bounding boxes
  { if (nrow(lab_df)) geom_text(
    data = lab_df,
    aes(x = x, y = y_fac, label = label),
    inherit.aes = FALSE,
    hjust = 0, size = 3.3, color = "grey20"
  ) } +
  # Cluster labels on the right margin
  scale_color_identity() +
  scale_fill_gradient2(
    name = NULL,
    midpoint = 0,
    low = "blue",    # negative rg -> blue
    mid = "white",   # rg = 0 -> white
    high = "red",    # positive rg -> red
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5),
    na.value = "white"
  ) +
  scale_x_discrete(position = "top") +
  coord_fixed(
    ratio = 1,
    xlim = c(0.5, length(trait_levels) + 0.5),
    clip = "off"  # allow cluster labels to render outside the panel
  ) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colorbar(
    direction = "horizontal",
    barheight = grid::unit(0.25, "cm"),
    barwidth  = grid::unit(0.20, "npc"),
    ticks = FALSE,
    frame.colour = "black",
    title.position = "top",
    title.hjust = 0.5,
    label.position = "bottom"
  )) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0,vjust=0),
    plot.title  = element_text(hjust = 0, face = "bold"),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    plot.margin = margin(10, 10, 10, 10)
  )

p_rg

## Heatmap interpretation:
# white/null: nominal p-value > 0.05
# small tile: nominal p-value < 0.05
# full tile:  global FDR q < 0.05
# star (*):   global Bonferroni p < 0.05


# Count features with significant genetic correlation with BMI (BMI-adjusted GWAS)
dt_BMI <- dt[ dt$trait == 'BMI',  ]
dim(dt_BMI) # 59 14
sum(p.adjust(dt_BMI$p, method = "fdr") < 0.05, na.rm = TRUE) # 44




### Summary ---------
## Bivariate LDSC with original (non-BMI-adjusted) liver feature GWAS:
# 1. Most liver features are influenced by BMI; they show pervasive genetic correlation with BMI.
# 2. Because BMI is genetically correlated with many traits, liver features also show broad genetic correlations.
# 3. BMI adjustment in GWAS is motivated by the above observations.

## Bivariate LDSC with BMI-adjusted liver feature GWAS:
# 1. Adjusting for BMI substantially reduces genetic correlations with many traits, yielding a cleaner pattern.
# 2. Genetic correlations with BMI itself largely disappear, as expected.

# Features with trait-specific genetic correlations (significant with one trait but not others)
# are of particular interest for downstream MR analyses.
