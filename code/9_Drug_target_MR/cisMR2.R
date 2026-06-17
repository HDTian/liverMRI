# =============================================================================
# Author: Haodong Tian
# Description: Polygenic and network MR for modifiable exposures (BMI, LDL-C,
#              triglycerides) on liver MRI features, comparing drug-target cis-MR
#              with polygenic MR using MR-Cluster to identify pathway-specific effects.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

### cisMR2:
# For certain modifiable traits (BMI, TG, LDL-C), this script examines
# the differences in their effects on liver MRI outcomes via different pathways,
# analogous to a network/mediation MR framework.

### Output figures:
# forest plot (BMI -> MRI)   ||  forest plot (LDL-C -> MRI)  forest plot (LDL-C -> CAD)   MR-cluster scatterplot
# forest plot (TG -> MRI1)       forest plot (TG -> MRI2)    forest plot (TG -> CAD)       MR-cluster scatterplot


library(data.table)
library(biomaRt)
library(rtracklayer)
library(susieR) # fine-mapping
library(Rfast)  # fast SuSiE in susieR
library(MendelianRandomization) # MR IVW
library(mr.raps) # MR-RAPS
library(dplyr)
library(forestplot)
library(grid)

library(mrclust) # install_github("cnfoley/mrclust")
library(dplyr)
library(patchwork)
library(ggplot2)
library(scales) # manually set ggplot color order so that polygenic MR "Overall" is consistently grey


### cisMR: Gene -> BMI/LDL-C -> liver MRI features: already DONE; use previous MR-PC-GMM results

### Polygenic MR: BMI/LDL-C => liver MRI features
# NOTE: unlike cis-MR-PC, winner's curse must be considered for polygenic MR



### load some cisMR results (MR-LD-IVW)
LDaware_cisMR_res <- readRDS("/path/to/project/targetMR/cisMR_res/LDaware_cisMR_res.rds")
Postivie_LDaware_cisMR_res <- readRDS("/path/to/project/targetMR/cisMR_res/Positive_LDaware_cisMR_res.rds")
### load 1000G LD reference panel
bfile  <- "/path/to/project/PheWAS/1000G_QC"  # 1000G path
bim <- fread(paste0(bfile, ".bim"), col.names = c("CHR","SNP","CM","BP","A1","A2"))


### ========================================================================================================
### BMI ====================================================================================================
### ========================================================================================================


## Step0: Decide the SNP selection strategy ================================================================
# Use GIANT clumped SNPs (82 SNPs) - avoid the winner's curse; also keeps <100 lead SNPs.
# 3136 lead SNPs in MVP / 82 lead SNPs in GIANT — use GIANT because:
#   (1) fewer SNPs, (2) safer against winner's curse, (3) consistent with prior MR studies.
# BMI GIANT data source: https://giant-consortium.web.broadinstitute.org/GIANT_consortium_data_files -> Download BMI EUR Ancestry GZIP


## Step1: PLINK LD-clumping to independent SNPs with pvalue < 5e-8 [only once!] =============================

# PLINK configuration
plink  <- "/path/to/project/PheWAS/plink_mac_20250615/plink"
bfile  <- "/path/to/project/PheWAS/1000G_QC"
# GWAS input and PLINK output path
indir  <- "/path/to/project/targetMR/ExposureGWAS/"
outdir <- file.path(indir, "LD_clumped_results")  # PLINK results will be stored here
# GWAS data for LD-clumping (use a separate dataset to avoid winner's curse)
#DT<- fread('/path/to/project/targetMR/ProxyGWAS/BMI.tsv')  # BMI UKB+GIANT
DT<- fread('/path/to/project/targetMR/ExposureGWAS/BMI_GIANT.tsv')
#DT<- fread(file.path(indir, paste0('BMI_MVP.tsv')))  # BMI MVP
# ---- auto-detect rsID and p-value columns (with priority) ----
rs_candidates <- c("rsID","rs_id","rsid","variant_id","ID","SNP")
p_candidates  <- c("p_value","P-value","pvalue","p")
rs_col <- intersect(rs_candidates, names(DT))[1]
p_col  <- intersect(p_candidates,  names(DT))[1]
if (is.na(rs_col) || is.na(p_col)) {stop("Cannot find rsID or p-value column in ", basename(f))}
# ---- create PLINK-format dummy GWAS file (rsID + p) ----
tmp <- tempfile(fileext = ".txt")
fwrite(DT[, .(SNP = get(rs_col), P = get(p_col))], tmp, sep = "\t")
# ---- output prefix ----
#prefix <- file.path(outdir, 'BMI_UKB_GIANT')
prefix <- file.path(outdir, 'BMI_GIANT')
#prefix <- file.path(outdir, 'BMI_MVP')
# ---- PLINK clumping ----
cmd <- sprintf(
  '%s --bfile %s --clump %s --clump-p1 5e-8 --clump-p2 1e-6 --clump-r2 0.01 --clump-kb 500 --maf 0.01 --out %s',
  shQuote(plink), shQuote(bfile), shQuote(tmp), shQuote(prefix)
)
system(cmd)  # Warning like 'rs12129899' is missing from the main dataset — due to MAF filtering


## Step 2: MR estimation: match liver MRI data -> harmonization -> MR fitting ===================================

# Get the BMI summary data (use BMI_MVP — must be the same data as used for MR-PC-GMM)
## ---- find clumped file ----
clump_file <- '/path/to/project/targetMR/ExposureGWAS/LD_clumped_results/BMI_GIANT.clumped'
## ---- extract rsIDs from .clumped ----
clumped_rs <- fread(clump_file, fill=TRUE)$SNP
clumped_rs <- clumped_rs[nzchar(clumped_rs)]   # remove trailing empty strings
## ---- read exposure GWAS (used for MR estimation) ----
exp_file <- '/path/to/project/targetMR/ExposureGWAS/BMI_MVP.tsv'
BMI_DT <- data.table::fread(exp_file)  # ~ 1 min
## ---- auto-detect columns (priority order) ----
rs_col  <- intersect(c("rsID","variant_id","rsid","rs_id","ID","SNP"), names(BMI_DT))[1]
ea_col  <- intersect(c("EA","alt","ALLELE1","effect_allele","Allele1"), names(BMI_DT))[1]
est_col <- intersect(c("beta","Effect","slope","BETA"), names(BMI_DT))[1]
se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(BMI_DT))[1]
## ---- restrict to clumped SNPs & standardize output ----
GXdata <- BMI_DT[get(rs_col) %in% clumped_rs, .(rsID = get(rs_col),  EA   = get(ea_col) , est  = get(est_col), se   = get(se_col) )]
GXdata <- na.omit(GXdata)   # remove NA rows  # 81 SNPs, no NA actually
dim(GXdata) # 81 4


# Get the polygenic MR results of BMI on nominated liver MRI outcomes ------------------------
radiomicsGWAS <- '/path/to/server/project/GWAS_regenie/MR_GWAS_new'  # all QC-ed GWAS results — 200 liver MRI features
feature_paths <- list.files(radiomicsGWAS, pattern = "\\.regenie$", full.names = TRUE)
feature_paths
feat_in_dir <- sub("\\.regenie$", "", basename(feature_paths))
liverMRInames <- c( 'firstorder_Minimum_inp'  )
feature_paths_selected <- feature_paths[feat_in_dir %in% liverMRInames]


BMIres<-list()
for (kk in 1:length(feature_paths_selected)  ) {

  feature_path<- feature_paths_selected[kk]
  feature_name <- sub("\\.regenie$", "", basename(feature_path))
  feature_name

  cat(  paste0( 'current liver MRI outcome (',kk, '/', length(feature_paths_selected),'): ', feature_name , '\n'  ) )


  ## get the outcome GWAS data (rsID EA est se) -----------------------------------------
  dat <- fread(feature_path)  # ~ 1 min
  #no worry: In close.connection(con) : Problem closing connection:  Bad file descriptor
  dat <- dat[, c("SNP", "ALLELE1", "BETA", "SE")]
  GYdata <- dat[match(GXdata$rsID ,   dat$SNP)  , ]
  GYdata <- na.omit(GYdata)   # remove NA rows

  ## harmonization (according to the LD ref panel effect allele) ------------------------
  GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]  # double-insurance: only consider rsIDs present in both GX and GY, in the same order
  bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
  by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
  ref <- bim[match(GXdata_$rsID, bim$SNP), ]
  flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1
  bx[flip_e] <- -bx[flip_e];  by[flip_o] <- -by[flip_o]


  ## get the rsID bx bxse by byse table
  MRtable<- data.frame( rsID = GXdata_$rsID  , bx=bx, bxse=bxse , by=by , byse=byse  )

  ## random-effect IVW (store the F statistic and Q test results) ----------------------
  # MRres<-mr_ivw(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse ) ) # weights="delta" # really need 2rd order?
  mr_plot( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse )  , orientate=TRUE)
  # still weak heterogeneity -> may be worth running MR-cluster

  ## MR-median (more robust to heterogeneity; hence better to represent the single overall MR effect)
  ## MRmedian<-mr_median( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse) )

  ## MR-Cluster -> groups variants into multiple valid clusters
  set.seed(1123)
  MRcluster <- mr_clust_em(theta = by/bx, theta_se = byse/abs(bx),
                       bx = bx, by = by, bxse = bxse, byse = byse, obs_names =GXdata_$rsID )
  MRcluster$plots$two_stage
  best <- MRcluster$results$best  # best-allocated cluster for each SNP

  # first remove junk variants
  best_nojunk <- best[  best$cluster_class != 'Junk', ]
  Unique_Cluste_index<-unique(best_nojunk$cluster)
  cat(  paste0(   'total non-junk cluster number indicated by MR-cluster: ' , length(Unique_Cluste_index  ), '\n'  ) )

  ### MR-IVW using all variants (after removing junk/outlier variants)
  MRtable_nojunk <- MRtable[  MRtable$rsID %in% best_nojunk$observation, ]
  MRres<-mr_ivw(mr_input(bx= MRtable_nojunk$bx, bxse= MRtable_nojunk$bxse ,by=MRtable_nojunk$by, byse=MRtable_nojunk$byse ) )

  MRIVWnojunk_vector <- c(   0 , MRres@SNPs ,  MRres@Estimate, MRres@StdError, MRres@CILower, MRres@CIUpper, MRres@Pvalue    )


  cluster_MR_matrix<-c()
  for( cc in  Unique_Cluste_index ){
    current_Cluster_rsIDs <-  best_nojunk$observation[ best_nojunk$cluster == cc   ]
    ## get the current rsID MR data
    current_MRtable<-MRtable[ match( current_Cluster_rsIDs, MRtable$rsID   ),  ]
    ## random-effect IVW
    current_MRres<-mr_ivw(mr_input(bx= current_MRtable$bx, bxse= current_MRtable$bxse ,by=current_MRtable$by, byse=current_MRtable$byse ) )

    vector_res <- c( cc , nrow(current_MRtable) , current_MRres@Estimate ,current_MRres@StdError , current_MRres@CILower ,current_MRres@CIUpper, current_MRres@Pvalue  )
    cluster_MR_matrix<-rbind(   cluster_MR_matrix  , vector_res      )
  }
  colnames(  cluster_MR_matrix) <- c( 'Cluster_ID', 'num_of_SNPs','est','se','CIlow','CIup','pvalue'  ) # all numeric

  ### combined MR-IVW and cluster-specific MR
  MR_matrix<- rbind(MRIVWnojunk_vector, cluster_MR_matrix  )

  ### store the result
  BMIres[[feature_name]] <- MR_matrix


}


saveRDS(BMIres,
        "/path/to/project/targetMR/cisMR2_res/BMI_polygenicMR_res.rds")
BMIres <- readRDS("/path/to/project/targetMR/cisMR2_res/BMI_polygenicMR_res.rds")



# forest plot / ggplot ==================================================================
liverMRInames <- c( 'firstorder_Minimum_inp'  )

BMIggplots <- list()
for( outcome_name in  liverMRInames ){


    # GLP1R cis-MR (MR-LD-IVW) results
    position <- match(outcome_name, LDaware_cisMR_res$GLP1R$feature_name)
    est <-  LDaware_cisMR_res$GLP1R$LD_aware_IVWest[position]
    se <- LDaware_cisMR_res$GLP1R$LD_aware_IVWse[position]
    GLP1R_MR_LD_IVW_res<- c(est, est - 1.96*se ,est +  1.96*se )
    # polygenic MR-IVW results
    BMIres_sub <-as.data.frame(BMIres[[outcome_name]])
    MR_IVW_res<- cbind( BMIres_sub$est , BMIres_sub$CIlow , BMIres_sub$CIup )
    # ggdata
    ggdata <- as.data.frame( rbind( GLP1R_MR_LD_IVW_res , MR_IVW_res ))
    colnames(ggdata) <- c('est','CIlow','CIup'  )
    ggdata$type <- c('GLP1R', 'Overall', paste('Cluster',  BMIres_sub$Cluster_ID[-1] )  )
    rownames(ggdata)<- NULL

    # if MR-Cluster returns only one non-junk cluster, drop the redundant cluster row (identical to overall)
    if( nrow(MR_IVW_res)==2 ){   ggdata <- ggdata[- nrow(ggdata),]     }


  ### draw ggplot --------------------------
  # 1. prepare data: classify into cis-MR vs. Polygenic MR categories
  df <- transform(ggdata, type = as.character(type))
  # GLP1R belongs to cis-MR; all others belong to Polygenic MR
  lev <- unique(df$type)
  df <- df %>%mutate( category = ifelse(type %in% c("GLP1R"), "cis-MR", "Polygenic MR"),
                      display_name = factor(type, levels = lev)  )
  # 2. ordering
  lev <- unique(df$type)
  lev <- c(intersect(c("GLP1R"), lev), setdiff(lev, c("GLP1R")))
  df$display_name <- factor(df$type, levels = lev)
  df$index <- nrow(df):1
  # 3. plotting
  is_drug <- df$category == "cis-MR"
  # manually set color order so that polygenic MR "Overall" is consistently grey
  levs <- if(is.factor(df$display_name)) levels(df$display_name) else sort(unique(as.character(df$display_name)))
  cols <- setNames(scales::hue_pal()(length(levs)), levs); cols["Overall"] <- "grey50"
  ggp<- ggplot(df, aes(x = est, y = index)) +
    geom_vline(xintercept = 0, linetype = 2) +
    # --- Polygenic MR ---
    geom_errorbarh(data = subset(df, !is_drug),aes(xmin = CIlow, xmax = CIup, color = display_name),width = 0, linewidth = 0.6) +
    geom_point(data = subset(df, !is_drug),aes(color = display_name),shape = 1, size = 3, stroke = 0.9) +
    # --- cis-MR ---
    geom_errorbarh(data = subset(df, is_drug),aes(xmin = CIlow, xmax = CIup), width = 0, linewidth = 0.6, color = "black") +
    geom_point(data = subset(df, is_drug),aes(shape = display_name),size = 3, color = "black") +
    # --- style: shape legend for cis-MR, color legend for polygenic MR ---
    scale_shape_manual(values = c(GLP1R = 15)) +
    scale_color_discrete(name = "Polygenic MR") +
    guides(shape = guide_legend(title = "cis-MR", order = 1),
           color = guide_legend(title = "Polygenic MR", order = 2)
    ) +
    labs(title = sprintf("BMI → %s", outcome_name), x = NULL, y = NULL)+   # →
    scale_y_continuous(expand = expansion(mult = c(0.5, 0.5)))+
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),  # hide the default L-shaped axis lines from classic theme (prevents uneven line thickness)
      axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
      legend.title = element_text(size = 10),
      legend.key.height = unit(0.5, "cm") )+ scale_color_manual(name = "Polygenic MR", values = cols)

  # store the ggplot
  BMIggplots[[outcome_name]]<- ggp


}

BMIggplots$firstorder_Minimum_inp



saveRDS(BMIggplots,
        "/path/to/project/targetMR/cisMR2_res/BMIggplots.rds")
BMIggplots <- readRDS("/path/to/project/targetMR/cisMR2_res/BMIggplots.rds")



# final plot
# BMIggplots$firstorder_Minimum_inp|BMIggplots$glszm_GrayLevelNonUniformity_inp
#BMI_MRPCAvsMRIVW  # 617 274


### ========================================================================================================
### LDL-C ====================================================================================================
### ========================================================================================================

### Step0: Decide the SNP selection strategy ================================================================
# Previously considered using 8 canonical drug-target SNPs (HMGCR, PCSK9, NPC1L1, APOB, ABCG5/G8, SORT1, APOE, LDLR),
# but 8 SNPs is too low-power.
# An AJHG benchmark paper shows that when GWAS sample size is large (>800,000), IV selection bias is negligible.
# Therefore, use the same large exposure GWAS for both IV selection and GX data collection,
# but limit to ~100 independent lead SNPs.

# rsIDs <- c(
#   "rs12916",   # HMGCR
#   "rs2479409", # PCSK9
#   "rs2072183", # NPC1L1
#   "rs1367117", # APOB
#   "rs4299376", # ABCG5/G8
#   "rs629301",  # SORT1
#   "rs4420638", # APOE
#   "rs6511720"  # LDLR
# )


### Step1: PLINK LD-clumping to independent SNPs with pvalue < 5e-8 [only once] =========================================

# PLINK configuration
plink  <- "/path/to/project/PheWAS/plink_mac_20250615/plink"
bfile  <- "/path/to/project/PheWAS/1000G_QC"
# GWAS input and PLINK output path
indir  <- "/path/to/project/targetMR/ExposureGWAS/"
outdir <- file.path(indir, "LD_clumped_results")  # PLINK results will be stored here
# GWAS data for LD-clumping
DT<- fread(file.path(indir, paste0('LDL_C.tsv')))  # GLGC LDL-C
# ---- auto-detect rsID and p-value columns (with priority) ----
rs_candidates <- c("rsID","rs_id","rsid","variant_id","ID","SNP")
p_candidates  <- c("p_value","P-value","pvalue","p")
rs_col <- intersect(rs_candidates, names(DT))[1]
p_col  <- intersect(p_candidates,  names(DT))[1]
if (is.na(rs_col) || is.na(p_col)) {stop("Cannot find rsID or p-value column in ", basename(f))}
# ---- re-derive p-values from Z-scores to avoid truncated p-values (e.g., p=0) in PLINK input ----
# Many large GWAS report p=0 for very significant variants; PLINK cannot rank these correctly.
# Solution: compute Z = beta/se, re-derive p = 2*pnorm(-|Z|), add tiny perturbations to preserve rank order.
est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
if (is.na(est_col) || is.na(se_col)) {stop("Cannot find est or se column (for |z-score| rank) in ", basename(f))}
DT[, z := get(est_col) / get(se_col)]
DT <- DT[order(-abs(z))]  # sort by descending |Z|; critical for newP rank stability below
DT[, newP := 2 * pnorm(-abs(z))] # re-derive p-value from Z-score for PLINK clumping input
p0 <- 1e-300 # minimum p-value recognizable by PLINK
DT[!is.finite(newP) | newP < p0, newP := p0] # truncate at p0 and add small perturbations to stabilize rank
DT[newP == p0, newP := newP + seq_along(.I) * 1e-302]  # supports up to 100 rank-perturbed entries
# overwrite the original p-value column with the re-derived newP
DT[, (p_col) := newP]
# ---- create PLINK-format dummy GWAS file (rsID + p) ----
tmp <- tempfile(fileext = ".txt")
fwrite(DT[, .(SNP = get(rs_col), P = get(p_col))], tmp, sep = "\t")
# ---- output prefix ----
prefix <- file.path(outdir, 'LDL_C')
# ---- PLINK clumping ----
cmd <- sprintf(
  '%s --bfile %s --clump %s --clump-p1 5e-10 --clump-p2 5e-8 --clump-r2 0.01 --clump-kb 500 --maf 0.01 --out %s',
  shQuote(plink), shQuote(bfile), shQuote(tmp), shQuote(prefix)
)
system(cmd)  # Warning like 'rs12129899' is missing from the main dataset — due to MAF filtering
# 1037 lead SNPs in LDL_C # --clump-p1 5e-8  --clump-p2 1e-6
# 763                      # --clump-p1 5e-10 --clump-p2 5e-8




### Step 2: polygenic MR: 100 SNPs -> LDL-C -> CAD =============================================

clump_file <- '/path/to/project/targetMR/ExposureGWAS/LD_clumped_results/LDL_C.clumped'
## ---- extract rsIDs from .clumped ----
clumped_rs <- fread(clump_file, fill=TRUE)$SNP
rsIDs <- clumped_rs[nzchar(clumped_rs)]   # remove trailing empty strings


## ---- read exposure LDL-C GWAS (used for MR estimation) ----
LDL_C_file <- '/path/to/project/targetMR/ExposureGWAS/LDL_C.tsv'
LDL_C_DT <- data.table::fread(LDL_C_file)
## ---- auto-detect columns (priority order) ----
rs_col  <- intersect(c("rsID","variant_id","rsid","rs_id","ID","SNP"), names(LDL_C_DT))[1]
ea_col  <- intersect(c("EA","alt","ALLELE1","effect_allele","Allele1"), names(LDL_C_DT))[1]
est_col <- intersect(c("beta","Effect","slope","BETA"), names(LDL_C_DT))[1]
se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(LDL_C_DT))[1]
## ---- restrict to clumped SNPs & standardize output ----
GXdata <- LDL_C_DT[get(rs_col) %in% rsIDs, .(rsID = get(rs_col),  EA   = get(ea_col) , est  = get(est_col), se   = get(se_col) )]
GXdata <- na.omit(GXdata)   # remove NA rows
dim(GXdata)   # 765  4

# NOTE: must rank by |Z| = |est/se| rather than raw p-value, to correctly select the top 100 SNPs
# when some p-values are truncated to 0 in the GWAS file.
GXdata<- GXdata[  order(    -abs( GXdata$est/GXdata$se ) ) ,   ][1:min(100, nrow(GXdata)), ] # limit to first 100 most significant SNPs
GXdata$pvalue <- 2*pnorm( -abs( GXdata$est/GXdata$se  ) )  # verify p-values
# check chr:pos to confirm these are polygenic (distributed across different genes)
GXdata <- GXdata %>% left_join(bim[, .(SNP, CHR, BP)], by = c("rsID" = "SNP"))
# rs11591147    PCSK9 lead SNP
# rs12916       HMGCR lead SNP

# From JAMA paper: https://cdn.jamanetwork.com/ama/content_public/journal/jama/935764/joi160113t2.png
#                    G-LDL-C      G-CAD(OR)            G-CAD(logOR)                 MR: LDL-C on CAD(logOR)      G-T2D
# PCSK9 rs11591147:  -0.5        0.77(0.69 0.87)     -0.261 (-0.371, -0.139)        0.522 (0.278,0.742)          non-significant
# HMGCR rs12916   :  -0.07       0.97(0.95,0.98)     -0.030 (-0.051, -0.020)        0.429 (0.286,0.729)          significant

# From current data:
#                    G-LDL-C              G-CAD(logOR)                 MR: LDL-C on CAD(logOR)
# PCSK9 rs11591147:  -0.434       -0.2513 (-0.2957, -0.2069)        0.579  (0.477,0.681)
# HMGCR rs12916   :  -0.07        -0.0284 (-0.0382, -0.0186)        0.406  (0.266,0.546)



# Positive control: LDL-C on CAD ---------------------------------------------------------------------
CADdat <- fread('/path/to/project/targetMR/otherGWAS/CAD.tsv')
dat <- CADdat[, c("rsID", "effect_allele", "beta", "standard_error","effect_allele_frequency","n")]
names(dat) <- c("SNP", "ALLELE1", "BETA", "SE",'EAF','N')
GYdata <- dat[match(GXdata$rsID ,   dat$SNP)  , ]  # match retains NA for missing SNPs
GYdata <- na.omit(GYdata)   # remove NA rows
dim(GYdata)   # 99 6 [the first 100 most significant SNPs]

## harmonization (according to the LD ref panel effect allele) ------------------------
GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]
bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
ref <- bim[match(GXdata_$rsID, bim$SNP), ]
flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1
bx[flip_e] <- -bx[flip_e];by[flip_o] <- -by[flip_o]

## get the rsID bx bxse by byse table
MRtable<- data.frame( rsID = GXdata_$rsID  , bx=bx, bxse=bxse , by=by , byse=byse  )

## random-effect IVW (store the F statistic and Q test results) ----------------------
MRres<-mr_ivw(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse ) )
#  IVW    est=0.497     se=0.039 CI=(0.420, 0.574)
mr_plot( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse )  , orientate=TRUE)

## MR-Cluster (positive control: LDL-C on CAD) --------------------
set.seed(1123)
MRcluster <- mr_clust_em(theta = by/bx, theta_se = byse/abs(bx),
                         bx = bx, by = by, bxse = bxse, byse = byse, obs_names = GXdata_$rsID )
View(MRcluster$results$all)
# Cluster1 mean: 0.618  # Cluster2 mean: 0.316
# Very good results; close to our findings on LDL-C -> CAD via PCSK9 and HMGCR
MRcluster$plots$two_stage

best <- MRcluster$results$best  # best-allocated cluster for each SNP

# first remove junk variants
best_nojunk <- best[  best$cluster_class != 'Junk', ]
Unique_Cluste_index<-unique(best_nojunk$cluster)
cat(  paste0(   'total non-junk cluster number indicated by MR-cluster: ' , length(Unique_Cluste_index  ), '\n'  ) )

### MR-IVW using all variants (after removing junk/outlier variants)
MRtable_nojunk <- MRtable[  MRtable$rsID %in% best_nojunk$observation, ]
MRres<-mr_ivw(mr_input(bx= MRtable_nojunk$bx, bxse= MRtable_nojunk$bxse ,by=MRtable_nojunk$by, byse=MRtable_nojunk$byse ) )
# est=0.429     se=0.020 CI=(0.390, 0.468)
MRIVWnojunk_vector <- c(   0 , MRres@SNPs ,  MRres@Estimate, MRres@StdError, MRres@CILower, MRres@CIUpper, MRres@Pvalue    )

cluster_MR_matrix<-c()
for( cc in  Unique_Cluste_index ){
  current_Cluster_rsIDs <-  best_nojunk$observation[ best_nojunk$cluster == cc   ]
  ## get the current rsID MR data
  current_MRtable<-MRtable[ match( current_Cluster_rsIDs, MRtable$rsID   ),  ]
  ## random-effect IVW
  current_MRres<-mr_ivw(mr_input(bx= current_MRtable$bx, bxse= current_MRtable$bxse ,by=current_MRtable$by, byse=current_MRtable$byse ) )

  vector_res <- c( cc , nrow(current_MRtable) , current_MRres@Estimate ,current_MRres@StdError , current_MRres@CILower ,current_MRres@CIUpper, current_MRres@Pvalue  )
  cluster_MR_matrix<-rbind(   cluster_MR_matrix  , vector_res      )
}
colnames(  cluster_MR_matrix) <- c( 'Cluster_ID', 'num_of_SNPs','est','se','CIlow','CIup','pvalue'  ) # all numeric

### combined MR-IVW and cluster-specific MR
MR_matrix<- rbind(MRIVWnojunk_vector, cluster_MR_matrix )

LDL_C_MRcluster_matrix <- cluster_MR_matrix[-3,] # remove the null cluster (cluster_ID = 3)
# Remove the null cluster (cluster_ID = 3): LDL-C has a known causal effect on CAD,
# so mr_clust_em does not need to assign a null cluster
LDL_C_polygenicMR_matrix<- rbind(MRIVWnojunk_vector, LDL_C_MRcluster_matrix )


# save: added to the LDLres list
#LDLres <- readRDS("/path/to/project/targetMR/cisMR2_res/LDLres.rds")
#LDLres[['CAD']]<-LDL_C_polygenicMR_matrix




### appendix figure: polygenic MR scatter plot — LDL-C on CAD ------------------------------------------------------
# scatter plot + two clusters + junk & null cluster shown in grey + cluster-specific effects
# + rs11591147 (PCSK9) and rs12916 (HMGCR) highlighted
library(data.table);library(ggplot2);library(ggrepel)
plot_df<-merge(best[,.(rsID=observation,cluster)],MRtable,by="rsID")
plot_df[bx<0,`:=`(bx=-bx,by=-by)]
plot_df[,cluster_plot:=fifelse(cluster%in%c(1,2),paste0("Cluster ",cluster),"Junk/Null")]
plot_df[,highlight:=rsID%in%c("rs11591147","rs12916")]
# cluster-specific MR slopes
cluster_line<-data.table(cluster_plot=c("Cluster 2","Cluster 1"),
                         slope=c(LDL_C_MRcluster_matrix[1,3],LDL_C_MRcluster_matrix[2,3]))
cols<-c("Cluster 1"="#D55E00","Cluster 2"="#0072B2","Junk/Null"="grey70")
# final ggplot
ggplot(plot_df,aes(bx,by,color=cluster_plot))+
  geom_vline(xintercept=0,linetype=2,color="grey50")+geom_hline(yintercept=0,linetype=2,color="grey50")+
  geom_errorbar(aes(ymin=by-1.96*byse,ymax=by+1.96*byse),width=0,alpha=.35)+
  geom_errorbarh(aes(xmin=bx-1.96*bxse,xmax=bx+1.96*bxse),width=0,alpha=.35)+
  geom_point(size=2)+
  geom_abline(data=cluster_line,aes(slope=slope,intercept=0,color=cluster_plot),linetype="dashed",linewidth=1)+
  geom_point(data=plot_df[highlight==TRUE], shape=21,size=4.2,stroke=1.5, aes(fill=cluster_plot),color="red")+
  geom_text_repel(data = plot_df[rsID=="rs12916"],aes(label=rsID),color="red", fontface="bold", size=4,
                  nudge_x = 0.5*diff(range(plot_df$bx)),nudge_y = -0.10*diff(range(plot_df$by)),
                  box.padding = 0.0, point.padding = 0.5,segment.color="red", segment.size=.4,min.segment.length = 0)+
  geom_text_repel(data = plot_df[rsID=="rs11591147"],aes(label=rsID),color="red", fontface="bold", size=4,
                  nudge_x = -0.05*diff(range(plot_df$bx)),nudge_y = -0.15*diff(range(plot_df$by)),
                  box.padding = 0.0, point.padding = 0.5,segment.color="red", segment.size=.4,min.segment.length = 0)+
  scale_color_manual(values=cols)+scale_fill_manual(values=cols)+
  labs(x="Genetic association with LDL-C",y="Genetic association with CAD (logOR)",color="Cluster")+
  theme_classic()+ guides(fill = "none")
# save as PDF (due to alpha=.35)

#rs12916     Cluster 2 posterior: 0.969
#rs11591147  Cluster 1 posterior: 0.998
# --------------------------------------------------------------------------------


### STEP3: run the polygenic MR for nominated liver MRI outcomes => LDLres ==========================================
radiomicsGWAS <- '/path/to/server/project/GWAS_regenie/MR_GWAS_new'  # all QC-ed GWAS results — 200 liver MRI features
feature_paths <- list.files(radiomicsGWAS, pattern = "\\.regenie$", full.names = TRUE)
feature_paths

feat_in_dir <- sub("\\.regenie$", "", basename(feature_paths))
# liverMRInames <- highlight_tbl$Outcome[highlight_tbl$gene%in% c('PCSK9','HMGCR')  ]
liverMRInames <- 'glszm_ZoneEntropy_inp'
feature_paths_selected <- feature_paths[feat_in_dir %in% liverMRInames]
feature_paths_selected
# 9 paths/files in total that are affected by PCSK9 and HMGCR

dim(GXdata)  # 100 7  # GXdata must be the LDL-C GWAS data at this point

LDLres<-list()
for (kk in 1:length(feature_paths_selected)  ) {

  feature_path<- feature_paths_selected[kk]
  feature_name <- sub("\\.regenie$", "", basename(feature_path))
  feature_name

  cat(  paste0( 'current liver MRI outcome (',kk, '/', length(feature_paths_selected),'): ', feature_name , '\n'  ) )

  ## get the outcome GWAS data (rsID EA est se) -----------------------------------------
  dat <- fread(feature_path)  # ~ 1 min; depending on VPN connection
  #no worry: In close.connection(con) : Problem closing connection:  Bad file descriptor
  dat <- dat[, c("SNP", "ALLELE1", "BETA", "SE")]
  GYdata <- dat[match(GXdata$rsID ,   dat$SNP)  , ]
  GYdata <- na.omit(GYdata)   # remove NA rows

  ## harmonization (according to the LD ref panel effect allele) ------------------------
  GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]  # double-insurance: only consider rsIDs present in both GX and GY, in the same order
  bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
  by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
  ref <- bim[match(GXdata_$rsID, bim$SNP), ]
  flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1
  bx[flip_e] <- -bx[flip_e];  by[flip_o] <- -by[flip_o]

  ## get the rsID bx bxse by byse table
  MRtable<- data.frame( rsID = GXdata_$rsID  , bx=bx, bxse=bxse , by=by , byse=byse  )

  ## random-effect IVW (store the F statistic and Q test results) ----------------------
  #MRres<-mr_ivw(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse ) ) # weights="delta" # really need 2rd order?
  mr_plot( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse )  , orientate=TRUE)

  ## MR-median (more robust to heterogeneity; hence better to represent the single overall MR effect)
  # MRmedian<-mr_median( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse) )

  ## MR-Cluster
  set.seed(1123)
  MRcluster <- mr_clust_em(theta = by/bx, theta_se = byse/abs(bx),
                           bx = bx, by = by, bxse = bxse, byse = byse, obs_names =GXdata_$rsID ) # ~ 5 mins
  MRcluster$plots$two_stage
  best <- MRcluster$results$best  # best-allocated cluster for each SNP

  # first remove junk variants
  best_nojunk <- best[  best$cluster_class != 'Junk', ]
  Unique_Cluste_index<-unique(best_nojunk$cluster)
  cat(  paste0(   'total non-junk cluster number indicated by MR-cluster: ' , length(Unique_Cluste_index  ), '\n'  ) )
  ### MR-IVW using all variants (after removing junk/outlier variants)
  MRtable_nojunk <- MRtable[  MRtable$rsID %in% best_nojunk$observation, ]
  MRres<-mr_ivw(mr_input(bx= MRtable_nojunk$bx, bxse= MRtable_nojunk$bxse ,by=MRtable_nojunk$by, byse=MRtable_nojunk$byse ) )

  MRIVWnojunk_vector <- c(   0 , MRres@SNPs ,  MRres@Estimate, MRres@StdError, MRres@CILower, MRres@CIUpper, MRres@Pvalue    )


  cluster_MR_matrix<-c()
  for( cc in  Unique_Cluste_index ){
    current_Cluster_rsIDs <-  best_nojunk$observation[ best_nojunk$cluster == cc   ]
    ## get the current rsID MR data
    current_MRtable<-MRtable[ match( current_Cluster_rsIDs, MRtable$rsID   ),  ]
    ## random-effect IVW
    current_MRres<-mr_ivw(mr_input(bx= current_MRtable$bx, bxse= current_MRtable$bxse ,by=current_MRtable$by, byse=current_MRtable$byse ) )

    vector_res <- c( cc , nrow(current_MRtable) , current_MRres@Estimate ,current_MRres@StdError , current_MRres@CILower ,current_MRres@CIUpper, current_MRres@Pvalue  )
    cluster_MR_matrix<-rbind(   cluster_MR_matrix  , vector_res      )
  }
  colnames(  cluster_MR_matrix) <- c( 'Cluster_ID', 'num_of_SNPs','est','se','CIlow','CIup','pvalue'  ) # all numeric

  ### combined MR-IVW and cluster-specific MR
  MR_matrix<- rbind(MRIVWnojunk_vector, cluster_MR_matrix  )

  ### store the result
  LDLres[[feature_name]] <- MR_matrix
}



# glszm_GrayLevelNonUniformity_inp (kk=2): good MR results, same direction as HMGCR
# ngtdm_Busyness (kk=7): good MR results, same direction as HMGCR

# Append the polygenic MR LDL-C -> CAD results
LDLres[['CAD']]<-LDL_C_polygenicMR_matrix   # overall + cluster-specific MR (null cluster excluded)

saveRDS(LDLres,
        "/path/to/project/targetMR/cisMR2_res/LDL_polygenicMR_res.rds")
LDLres <- readRDS("/path/to/project/targetMR/cisMR2_res/LDL_polygenicMR_res.rds")




### STEP4: forest/ggplot => LDLggplots ==================================================================

# For any given outcome, draw dot-line plots comparing PCSK9/HMGCR (cis-MR) with polygenic MR (est, CIlow, CIup)


outcome_name<-'glszm_GrayLevelNonUniformity_inp'  # good MR results
outcome_name<-'ngtdm_Busyness'
outcome_name<-'ngtdm_Complexity'
outcome_name <- liverMRinames[7]


liverMRInames<- 'glszm_ZoneEntropy_inp'
all_names <- c(liverMRInames, 'CAD'   )

LDLggplots <- list()
for( outcome_name in  all_names ){


  if( outcome_name!= 'CAD'  ){
    # PCSK9 cis-MR (MR-LD-IVW) results
    position <- match(outcome_name, LDaware_cisMR_res$PCSK9$feature_name )
    est <-LDaware_cisMR_res$PCSK9$LD_aware_IVWest[position]; se <- LDaware_cisMR_res$PCSK9$LD_aware_IVWse[position]
    PCSK9_cisMR_res<- c( est ,  est- 1.96*se, est+ 1.96*se )
    # HMGCR cis-MR (MR-LD-IVW) results
    position <- match(outcome_name, LDaware_cisMR_res$HMGCR$feature_name )
    est <-LDaware_cisMR_res$HMGCR$LD_aware_IVWest[position]; se <- LDaware_cisMR_res$HMGCR$LD_aware_IVWse[position]
    HMGCR_cisMR_res<- c( est ,  est- 1.96*se, est+ 1.96*se )
    # polygenic MR results
    LDLres_sub <-as.data.frame(LDLres[[outcome_name]])
    MR_IVW_res<- cbind( LDLres_sub$est , LDLres_sub$CIlow , LDLres_sub$CIup )
    # ggdata
    ggdata <- as.data.frame( rbind( PCSK9_cisMR_res, HMGCR_cisMR_res , MR_IVW_res ))
    colnames(ggdata) <- c('est','CIlow','CIup'  )
    ggdata$type <- c('PCSK9','HMGCR', 'Overall', paste('Cluster',  LDLres_sub$Cluster_ID[-1] )  )
    rownames(ggdata)<- NULL
  }else{
    # [if outcome is CAD] -------------
    # PCSK9/HMGCR cis-MR (MR-LD-IVW) results (positive control)
    positions <- match(c('PCSK9','HMGCR'  )  , Postivie_LDaware_cisMR_res$gene)
    PCSK9_HMGCR_cisMR_res<-cbind( Postivie_LDaware_cisMR_res$Estimate[positions],
                                    Postivie_LDaware_cisMR_res$CI_low[positions],
                                    Postivie_LDaware_cisMR_res$CI_up[positions])
    # polygenic MR-IVW results
    outcome_name<-'CAD'
    LDLres_sub <-as.data.frame(LDLres[[outcome_name]])
    MR_IVW_res<- cbind( LDLres_sub$est , LDLres_sub$CIlow , LDLres_sub$CIup )
    # ggdata
    ggdata <- as.data.frame( rbind( PCSK9_HMGCR_cisMR_res , MR_IVW_res ))
    colnames(ggdata) <- c('est','CIlow','CIup'  )
    ggdata$type <- c('PCSK9','HMGCR', 'Overall', paste('Cluster',  LDLres_sub$Cluster_ID[-1] )  )
    rownames(ggdata)<- NULL
  }


  ### draw ggplot --------------------------
  # 1. prepare data: classify into cis-MR vs. Polygenic MR categories
  df <- transform(ggdata, type = as.character(type))
  # PCSK9/HMGCR belong to cis-MR; all others belong to Polygenic MR
  df <- df %>%mutate( category = ifelse(type %in% c("PCSK9", "HMGCR"), "cis-MR", "Polygenic MR"),
                      display_name = factor(type, levels = c("PCSK9", "HMGCR") )  )
  # 2. ordering
  lev <- unique(df$type)
  lev <- c(intersect(c("PCSK9","HMGCR"), lev), setdiff(lev, c("PCSK9","HMGCR")))
  df$display_name <- factor(df$type, levels = lev)
  df$index <- nrow(df):1
  # 3. plotting
  is_drug <- df$category == "cis-MR"
  # manually set color order so that polygenic MR "Overall" is consistently grey
  levs <- if(is.factor(df$display_name)) levels(df$display_name) else sort(unique(as.character(df$display_name)))
  cols <- setNames(scales::hue_pal()(length(levs)), levs); cols["Overall"] <- "grey50"
  ggp<- ggplot(df, aes(x = est, y = index)) +
    geom_vline(xintercept = 0, linetype = 2) +
    # --- Polygenic MR ---
    geom_errorbarh(data = subset(df, !is_drug),aes(xmin = CIlow, xmax = CIup, color = display_name),width = 0, linewidth = 0.6) +
    geom_point(data = subset(df, !is_drug),aes(color = display_name),shape = 1, size = 3, stroke = 0.9) +
    # --- cis-MR ---
    geom_errorbarh(data = subset(df, is_drug),aes(xmin = CIlow, xmax = CIup), width = 0, linewidth = 0.6, color = "black") +
    geom_point(data = subset(df, is_drug),aes(shape = display_name),size = 3, color = "black") +
    # --- style: shape legend for cis-MR, color legend for polygenic MR ---
    scale_shape_manual(values = c(PCSK9 = 15, HMGCR = 17)) +
    scale_color_discrete(name = "Polygenic MR") +
    guides(shape = guide_legend(title = "cis-MR", order = 1),
           color = guide_legend(title = "Polygenic MR", order = 2)
    ) +
    labs(title = sprintf("LDL-C → %s", outcome_name), x = NULL, y = NULL)+   # →
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))+
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),  # hide the default L-shaped axis lines from classic theme (prevents uneven line thickness)
      axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
      legend.title = element_text(size = 10),
      legend.key.height = unit(0.5, "cm") ) + scale_color_manual(name = "Polygenic MR", values = cols)

  # store the ggplot
  LDLggplots[[outcome_name]]<- ggp


}


LDLggplots$glszm_ZoneEntropy_inp
LDLggplots$CAD

saveRDS(LDLggplots,
        "/path/to/project/targetMR/cisMR2_res/LDLggplots.rds")
LDLggplots <- readRDS("/path/to/project/targetMR/cisMR2_res/LDLggplots.rds")



# LDLggplots$glszm_GrayLevelNonUniformity_inp  # PCSK9 is null; but HMGCR is significant
# LDLggplots$glszm_SizeZoneNonUniformity_inp   # PCSK9 is null; but HMGCR is significant
#
# LDLggplots$ngtdm_Busyness                    # both PCSK9 and HMGCR are significant
# LDLggplots$ngtdm_Strength_inp                # HMGCR is null; but PCSK9 is significant
# LDLggplots$glcm_Correlation_inp              # PCSK9 is null; but HMGCR is significant
# LDLggplots$shape_Maximum3DDiameter
# LDLggplots$firstorder_Skewness_inp
# LDLggplots$firstorder_90Percentile

LDLggplots$CAD


### Which outcomes to present?
# Just pick representative liver MRI outcomes covering different scenarios:
# Scenario A: PCSK9 and HMGCR diverge; matched by two MR-cluster clusters   [glszm_GrayLevelNonUniformity_inp]
# Scenario B: PCSK9 and HMGCR are similar, but polygenic MR is null          [ngtdm_Complexity]
# + [CAD] as positive control


# ### final plot
# LDLggplots$glszm_GrayLevelNonUniformity_inp |
#   LDLggplots$ngtdm_Strength_inp  |
#   LDLggplots$shape_Maximum3DDiameter |
#   LDLggplots$CAD
# #LDL_MRPCAvsMRIVW  # 1234 274


### ========================================================================================================
### TG ====================================================================================================
### ========================================================================================================


## Step0: Decide the SNP selection strategy ================================================================

# An AJHG benchmark paper shows that when GWAS sample size is large (>800,000), IV selection bias is negligible.
# Therefore, use the same large exposure GWAS for both IV selection and GX data collection,
# but limit to ~100 independent lead SNPs.


# Step1: PLINK LD-clumping to independent SNPs with pvalue < 5e-8 [only once] =========================================

# PLINK configuration
plink  <- "/path/to/project/PheWAS/plink_mac_20250615/plink"
bfile  <- "/path/to/project/PheWAS/1000G_QC"
# GWAS input and PLINK output path
indir  <- "/path/to/project/targetMR/ExposureGWAS/"
outdir <- file.path(indir, "LD_clumped_results")  # PLINK results will be stored here
# GWAS data for LD-clumping
DT<- fread(file.path(indir, paste0('TG.tsv')))  # GLGC TG
# ---- auto-detect rsID and p-value columns (with priority) ----
rs_candidates <- c("rsID","rs_id","rsid","variant_id","ID","SNP")
p_candidates  <- c("p_value","P-value","pvalue","p")
rs_col <- intersect(rs_candidates, names(DT))[1]
p_col  <- intersect(p_candidates,  names(DT))[1]
if (is.na(rs_col) || is.na(p_col)) {stop("Cannot find rsID or p-value column in ", basename(f))}
# ---- re-derive p-values from Z-scores to avoid truncated p-values (e.g., p=0) in PLINK input ----
# (same strategy as LDL-C above)
est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
if (is.na(est_col) || is.na(se_col)) {stop("Cannot find est or se column (for |z-score| rank) in ", basename(f))}
DT[, z := get(est_col) / get(se_col)]
DT <- DT[order(-abs(z))]  # sort by descending |Z|; critical for newP rank stability below
DT[, newP := 2 * pnorm(-abs(z))] # re-derive p-value from Z-score for PLINK clumping input
p0 <- 1e-300 # minimum p-value recognizable by PLINK
DT[!is.finite(newP) | newP < p0, newP := p0] # truncate at p0 and add small perturbations to stabilize rank
DT[newP == p0, newP := newP + seq_along(.I) * 1e-302]  # supports up to 100 rank-perturbed entries
# overwrite the original p-value column with the re-derived newP
DT[, (p_col) := newP]
# ---- create PLINK-format dummy GWAS file (rsID + p) ----
tmp <- tempfile(fileext = ".txt")
fwrite(DT[, .(SNP = get(rs_col), P = get(p_col))], tmp, sep = "\t")
# ---- output prefix ----
prefix <- file.path(outdir, 'TG')
# ---- PLINK clumping ----
cmd <- sprintf(
  '%s --bfile %s --clump %s --clump-p1 5e-10 --clump-p2 5e-8 --clump-r2 0.01 --clump-kb 500 --maf 0.01 --out %s',
  shQuote(plink), shQuote(bfile), shQuote(tmp), shQuote(prefix)
)
system(cmd)  # Warning like 'rs12129899' is missing from the main dataset — due to MAF filtering

# 807 lead SNPs # --clump-p1 5e-10 --clump-p2 5e-8




## Step 2: polygenic MR: 100 SNPs -> TG -> CAD  ==========================================================

clump_file <- '/path/to/project/targetMR/ExposureGWAS/LD_clumped_results/TG.clumped'
## ---- extract rsIDs from .clumped ----
clumped_rs <- fread(clump_file, fill=TRUE)$SNP
rsIDs <- clumped_rs[nzchar(clumped_rs)]   # remove trailing empty strings


## ---- read exposure TG GWAS (used for MR estimation) ----
TG_file <- '/path/to/project/targetMR/ExposureGWAS/TG.tsv'
TG_DT <- data.table::fread(TG_file)  # ~ 1 min
## ---- auto-detect columns (priority order) ----
rs_col  <- intersect(c("rsID","variant_id","rsid","rs_id","ID","SNP"), names(TG_DT))[1]
ea_col  <- intersect(c("EA","alt","ALLELE1","effect_allele","Allele1"), names(TG_DT))[1]
est_col <- intersect(c("beta","Effect","slope","BETA"), names(TG_DT))[1]
se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(TG_DT))[1]
## ---- restrict to clumped SNPs & standardize output ----
GXdata0 <- TG_DT[get(rs_col) %in% rsIDs, .(rsID = get(rs_col),  EA   = get(ea_col) , est  = get(est_col), se   = get(se_col) )]
GXdata0 <- na.omit(GXdata0)   # remove NA rows
dim(GXdata0)   # 807  4

# NOTE: must rank by |Z| = |est/se| to correctly select top 100 SNPs (avoids truncated p-value issue)
GXdata<- GXdata0[  order(    -abs( GXdata0$est/GXdata0$se ) ) ,   ][1:min(100, nrow(GXdata0)), ] # limit to first 100 most significant SNPs
GXdata$pvalue <- 2*pnorm( -abs( GXdata$est/GXdata$se  ) )  # verify p-values

# rs11207980 (ANGPTL3 lead SNP) was not captured by PLINK clumping (likely overshadowed by a nearby lead SNP).
# Force-add it to ensure ANGPTL3 pathway is represented.
GXdata <- rbind( GXdata,
  TG_DT[get(rs_col)=="rs11207980",
        .(rsID=get(rs_col), EA=get(ea_col), est=get(est_col), se=get(se_col))] %>%
    na.omit() %>%  mutate(pvalue = 2*pnorm(-abs(est/se)))   )
dim(GXdata)  # 101 (= 100 + rs11207980) x 5
# check chr:pos to confirm these are polygenic
GXdata <- GXdata %>% left_join(bim[, .(SNP, CHR, BP)], by = c("rsID" = "SNP"))
# rs10889333 removed as it is likely in LD with rs11207980
GXdata <- GXdata[- match('rs10889333', GXdata$rsID  ) ,  ]
dim(GXdata) # 100 7

# Positive control: TG on CAD ---------------------------------------------------------------------
CADdat <- fread('/path/to/project/targetMR/otherGWAS/CAD.tsv')
dat <- CADdat[, c("rsID", "effect_allele", "beta", "standard_error","effect_allele_frequency","n")]
names(dat) <- c("SNP", "ALLELE1", "BETA", "SE",'EAF','N')
GYdata <- dat[match(GXdata$rsID ,   dat$SNP)  , ]  # match retains NA for missing SNPs
GYdata <- na.omit(GYdata)   # remove NA rows
dim(GYdata)  # 100 6 [the first 100 most significant SNPs]

## harmonization (according to the LD ref panel effect allele) ------------------------
GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]
bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
ref <- bim[match(GXdata_$rsID, bim$SNP), ]
flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1
bx[flip_e] <- -bx[flip_e];by[flip_o] <- -by[flip_o]

## get the rsID bx bxse by byse table
MRtable<- data.frame( rsID = GXdata_$rsID  , bx=bx, bxse=bxse , by=by , byse=byse  )

## random-effect IVW (store the F statistic and Q test results) ----------------------
MRres<-mr_ivw(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse ) )
#  IVW    est=0.291     se=0.040 CI=(0.213, 0.370)    # 100 SNPs -> TG -> CAD
mr_plot( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse )  , orientate=TRUE)

## MR-Cluster (positive control: TG on CAD) --------------------
set.seed(1123)
TG_CAD_MRcluster <- mr_clust_em(theta = by/bx, theta_se = byse/abs(bx),
                         bx = bx, by = by, bxse = bxse, byse = byse, obs_names = GXdata_$rsID ) # ~ 5 mins
MRcluster<-TG_CAD_MRcluster
MRcluster$results$all
# Cluster1 mean: 1.104
# Cluster2 mean: 0.236
# Cluster3 mean: 0.570
# Cluster4 mean: 1.680
MRcluster$plots$two_stage

best <- MRcluster$results$best  # best-allocated cluster for each SNP

# first remove junk variants
best_nojunk <- best[  best$cluster_class != 'Junk', ]
Unique_Cluste_index<-unique(best_nojunk$cluster)
cat(  paste0(   'total non-junk cluster number indicated by MR-cluster: ' , length(Unique_Cluste_index  ), '\n'  ) )

### MR-IVW using all variants (after removing junk/outlier variants)
MRtable_nojunk <- MRtable[  MRtable$rsID %in% best_nojunk$observation, ]
MRres<-mr_ivw(mr_input(bx= MRtable_nojunk$bx, bxse= MRtable_nojunk$bxse ,by=MRtable_nojunk$by, byse=MRtable_nojunk$byse ) )
# est=0.306     se=0.039 CI=(0.230, 0.382)
MRIVWnojunk_vector <- c(   0 , MRres@SNPs ,  MRres@Estimate, MRres@StdError, MRres@CILower, MRres@CIUpper, MRres@Pvalue    )

cluster_MR_matrix<-c()
for( cc in  Unique_Cluste_index ){
  current_Cluster_rsIDs <-  best_nojunk$observation[ best_nojunk$cluster == cc   ]
  ## get the current rsID MR data
  current_MRtable<-MRtable[ match( current_Cluster_rsIDs, MRtable$rsID   ),  ]
  ## random-effect IVW
  current_MRres<-mr_ivw(mr_input(bx= current_MRtable$bx, bxse= current_MRtable$bxse ,by=current_MRtable$by, byse=current_MRtable$byse ) )

  vector_res <- c( cc , nrow(current_MRtable) , current_MRres@Estimate ,current_MRres@StdError , current_MRres@CILower ,current_MRres@CIUpper, current_MRres@Pvalue  )
  cluster_MR_matrix<-rbind(   cluster_MR_matrix  , vector_res      )
}
colnames(  cluster_MR_matrix) <- c( 'Cluster_ID', 'num_of_SNPs','est','se','CIlow','CIup','pvalue'  ) # all numeric

### combined MR-IVW and cluster-specific MR
MR_matrix<- rbind(MRIVWnojunk_vector, cluster_MR_matrix )

TGres <- readRDS("/path/to/project/targetMR/cisMR2_res/TG_polygenicMR_res.rds")
TGres$CAD<- MR_matrix



### appendix figure: polygenic MR scatter plot — TG on CAD ----------------------------------------------
# scatter plot + four clusters + junk & null cluster shown in grey + cluster-specific effects
# + rs964184 (APOC3 lead SNP) and rs11207980 (ANGPTL3 lead SNP) highlighted
## merge cluster assignment with MR table
plot_df <- merge(  best[, .(rsID = observation, cluster)], MRtable, by = "rsID" )
## orient effects to positive bx
plot_df[bx < 0, `:=`(bx = -bx, by = -by)]
## cluster labels
plot_df[, cluster_plot := fifelse(cluster %in% 1:4, paste0("Cluster ", cluster), "Junk/Null")]
## highlight APOC3 and ANGPTL3 lead SNPs
plot_df[, highlight := rsID %in% c("rs964184","rs11207980")]
## cluster-specific MR slopes (clusters 1-4)
cluster_line <- data.table( cluster_plot = paste0("Cluster ", CAD_MRcluster_matrix[,1]), slope = CAD_MRcluster_matrix[, 3] )
## colors
cols <- c(
  "Cluster 1" = "#D55E00","Cluster 2" = "#0072B2","Cluster 3" = "#009E73","Cluster 4" = "#CC79A7","Junk/Null" = "grey70"
)

## final plot
ggplot(plot_df, aes(bx, by, color = cluster_plot)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbar(aes(ymin = by - 1.96 * byse, ymax = by + 1.96 * byse),width = 0, alpha = 0.35) +
  geom_errorbarh(aes(xmin = bx - 1.96 * bxse,xmax = bx + 1.96 * bxse),width = 0, alpha = 0.35) +
  geom_point(size = 2) +
  geom_abline( data = cluster_line,
    aes(slope = slope, intercept = 0, color = cluster_plot),
    linetype = "dashed", linewidth = 1) +
  geom_point(
    data = plot_df[highlight == TRUE],
    shape = 21, size = 4.2, stroke = 1.5,
    aes(fill = cluster_plot), color = "red") +
  geom_text_repel(
    data = plot_df[rsID == "rs964184"], aes(label = rsID),
    color = "red", fontface = "bold", size = 4,
    nudge_x = 0.5 * diff(range(plot_df$bx)),nudge_y = -0.10 * diff(range(plot_df$by)),
    box.padding = 0.0, point.padding = 0.5, segment.color = "red", segment.size = 0.4,min.segment.length = 0) +
  geom_text_repel(
    data = plot_df[rsID == "rs11207980"], aes(label = rsID),
    color = "red", fontface = "bold", size = 4,
    nudge_x = 0.5 * diff(range(plot_df$bx)),nudge_y = -0.10 * diff(range(plot_df$by)),
    box.padding = 0.0, point.padding = 0.5, segment.color = "red", segment.size = 0.4,min.segment.length = 0) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs( x = "Genetic association with TG", y = "Genetic association with CAD (log OR)",color = "Cluster") +
  theme_classic() +
  guides(fill = "none")
# save as PDF (due to alpha=.35)
# MRIVW_TG_CAD  # 5.24  4.65
###  appendix figure: polygenic MR: TG on CAD [DONE!] -------------------------------------------




### STEP3: run the polygenic MR for nominated liver MRI outcomes => TGres ==========================================

dim(GXdata)  # 100 7 # still based on TG GXdata
# Broad Server VPN required
radiomicsGWAS <- '/path/to/server/project/GWAS_regenie/MR_GWAS_new'  # all QC-ed GWAS results — 200 liver MRI features
feature_paths <- list.files(radiomicsGWAS, pattern = "\\.regenie$", full.names = TRUE)
feature_paths


liverMRInames <- c(  'glszm_ZoneEntropy_inp', 'glszm_GrayLevelVariance_inp'  )
feat_in_dir <- sub("\\.regenie$", "", basename(feature_paths))
feature_paths_selected <- feature_paths[feat_in_dir %in% liverMRInames]
feature_paths_selected
# 2 paths/files in total (nominated MRI outcomes affected by APOC3/ANGPTL3 via TG)

dim(GXdata)  # verify this is still TG GXdata (100 SNPs)

TGres<-list()
for (kk in 1:length(feature_paths_selected)  ) {

  feature_path<- feature_paths_selected[kk]
  feature_name <- sub("\\.regenie$", "", basename(feature_path))
  feature_name

  cat(  paste0( 'current liver MRI outcome (',kk, '/', length(feature_paths_selected),'): ', feature_name , '\n'  ) )

  ## get the outcome GWAS data (rsID EA est se) -----------------------------------------
  dat <- fread(feature_path)  # ~ 1 min; depending on VPN connection
  #no worry: In close.connection(con) : Problem closing connection:  Bad file descriptor
  dat <- dat[, c("SNP", "ALLELE1", "BETA", "SE")]
  GYdata <- dat[match(GXdata$rsID ,   dat$SNP)  , ]
  GYdata <- na.omit(GYdata)   # remove NA rows

  ## harmonization (according to the LD ref panel effect allele) ------------------------
  GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]  # double-insurance: only consider rsIDs present in both GX and GY, in the same order
  bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
  by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
  ref <- bim[match(GXdata_$rsID, bim$SNP), ]
  flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1
  bx[flip_e] <- -bx[flip_e];  by[flip_o] <- -by[flip_o]

  ## get the MRtable: rsID bx bxse by byse table
  MRtable<- data.frame( rsID = GXdata_$rsID  , bx=bx, bxse=bxse , by=by , byse=byse  )

  ## random-effect IVW (store the F statistic and Q test results) ----------------------
  #MRres<-mr_ivw(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse ) ) # weights="delta" # really need 2rd order?
  #mr_plot( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse )  , orientate=TRUE)

  ## MR-median (more robust to heterogeneity; hence better to represent the single overall MR effect)
  # MRmedian<-mr_median( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse) )

  ## MR-Cluster
  set.seed(1123)
  cat( 'MR clustering (roughly 5 mins) ... \n'  )
  MRcluster <- mr_clust_em(theta = by/bx, theta_se = byse/abs(bx),
                           bx = bx, by = by, bxse = bxse, byse = byse, obs_names =GXdata_$rsID ) # ~ 5 mins
  #MRcluster$plots$two_stage
  best <- MRcluster$results$best  # best-allocated cluster for each SNP

  # first remove junk variants
  best_nojunk <- best[  best$cluster_class != 'Junk', ]
  Unique_Cluste_index<-unique(best_nojunk$cluster)
  cat(  paste0(   'total non-junk cluster number indicated by MR-cluster: ' , length(Unique_Cluste_index  ), '\n'  ) )

  ### MR-IVW using all variants (after removing junk/outlier variants)
  MRtable_nojunk <- MRtable[  MRtable$rsID %in% best_nojunk$observation, ]
  MRres<-mr_ivw(mr_input(bx= MRtable_nojunk$bx, bxse= MRtable_nojunk$bxse ,by=MRtable_nojunk$by, byse=MRtable_nojunk$byse ) )

  MRIVWnojunk_vector <- c(   0 , MRres@SNPs ,  MRres@Estimate, MRres@StdError, MRres@CILower, MRres@CIUpper, MRres@Pvalue    )


  cluster_MR_matrix<-c()
  for( cc in  Unique_Cluste_index ){
    current_Cluster_rsIDs <-  best_nojunk$observation[ best_nojunk$cluster == cc   ]
    ## get the current rsID MR data
    current_MRtable<-MRtable[ match( current_Cluster_rsIDs, MRtable$rsID   ),  ]
    ## random-effect IVW
    current_MRres<-mr_ivw(mr_input(bx= current_MRtable$bx, bxse= current_MRtable$bxse ,by=current_MRtable$by, byse=current_MRtable$byse ) )

    vector_res <- c( cc , nrow(current_MRtable) , current_MRres@Estimate ,current_MRres@StdError , current_MRres@CILower ,current_MRres@CIUpper, current_MRres@Pvalue  )
    cluster_MR_matrix<-rbind(   cluster_MR_matrix  , vector_res      )
  }
  colnames(  cluster_MR_matrix) <- c( 'Cluster_ID', 'num_of_SNPs','est','se','CIlow','CIup','pvalue'  ) # all numeric

  ### combined MR-IVW and cluster-specific MR
  MR_matrix<- rbind(MRIVWnojunk_vector, cluster_MR_matrix  )


  ### store the result
  TGres[[feature_name]] <- MR_matrix
}




saveRDS(TGres,
        "/path/to/project/targetMR/cisMR2_res/TG_polygenicMR_res.rds")
TGres <- readRDS("/path/to/project/targetMR/cisMR2_res/TG_polygenicMR_res.rds")
# polygenic MR: all_SNPs -> TG -> nominated MRI outcomes + CAD




### STEP4: forest plot / ggplot: [TG] on [CAD / nominated MRI] => TGggplots ============================================
# using cis-MR (MR-LD-IVW) + polygenic MR results

outcome_name<-'glszm_ZoneEntropy_inp'

liverMRInames <- c( 'glszm_ZoneEntropy_inp', 'glszm_GrayLevelVariance_inp' ) # the MRI outcomes
all_names <- c(liverMRInames, 'CAD'   )

TGggplots <- list()
for( outcome_name in  all_names ){

  if( outcome_name!= 'CAD'  ){
    # APOC3 cis-MR (MR-LD-IVW) results
    position <- match(outcome_name, LDaware_cisMR_res$APOC3$feature_name )
    est <-LDaware_cisMR_res$APOC3$LD_aware_IVWest[position]; se <- LDaware_cisMR_res$APOC3$LD_aware_IVWse[position]
    APOC3_cisMR_res<- c( est ,  est- 1.96*se, est+ 1.96*se )
    # ANGPTL3 cis-MR (MR-LD-IVW) results
    position <- match(outcome_name, LDaware_cisMR_res$ANGPTL3$feature_name )
    est <-LDaware_cisMR_res$ANGPTL3$LD_aware_IVWest[position]; se <- LDaware_cisMR_res$ANGPTL3$LD_aware_IVWse[position]
    ANGPTL3_cisMR_res<- c( est ,  est- 1.96*se, est+ 1.96*se )
    # polygenic MR results
    TGres_sub <-as.data.frame(TGres[[outcome_name]])
    MR_IVW_res<- cbind( TGres_sub$est , TGres_sub$CIlow , TGres_sub$CIup )
    # ggdata
    ggdata <- as.data.frame( rbind( APOC3_cisMR_res, ANGPTL3_cisMR_res , MR_IVW_res ))
    colnames(ggdata) <- c('est','CIlow','CIup'  )
    ggdata$type <- c('APOC3','ANGPTL3', 'Overall', paste('Cluster',  TGres_sub$Cluster_ID[-1] )  )
    rownames(ggdata)<- NULL
  }else{
    # [if outcome is CAD] -------------
    # APOC3/ANGPTL3 cis-MR (MR-LD-IVW) results (positive control)
    positions <- match(c('APOC3','ANGPTL3'  )  , Postivie_LDaware_cisMR_res$gene)
    APOC3_ANGPTL3_cisMR_res<-cbind( Postivie_LDaware_cisMR_res$Estimate[positions],
                      Postivie_LDaware_cisMR_res$CI_low[positions],
                      Postivie_LDaware_cisMR_res$CI_up[positions])
    # polygenic MR-IVW results
    outcome_name<-'CAD'
    TGres_sub <-as.data.frame(TGres[[outcome_name]])
    MR_IVW_res<- cbind( TGres_sub$est , TGres_sub$CIlow , TGres_sub$CIup )
    # ggdata
    ggdata <- as.data.frame( rbind( APOC3_ANGPTL3_cisMR_res , MR_IVW_res ))
    colnames(ggdata) <- c('est','CIlow','CIup'  )
    ggdata$type <- c('APOC3','ANGPTL3', 'Overall', paste('Cluster',  TGres_sub$Cluster_ID[-1] )  )
    rownames(ggdata)<- NULL
  }


  ### draw ggplot --------------------------
  # 1. prepare data: classify into cis-MR vs. Polygenic MR categories
  df <- transform(ggdata, type = as.character(type))
  # APOC3/ANGPTL3 belong to cis-MR; all others belong to Polygenic MR
  df <- df %>%mutate( category = ifelse(type %in% c("APOC3", "ANGPTL3"), "cis-MR", "Polygenic MR"),
                      display_name = factor(type, levels = c("APOC3", "ANGPTL3") )  )
  # 2. ordering
  lev <- unique(df$type)
  lev <- c(intersect(c("APOC3","ANGPTL3"), lev), setdiff(lev, c("APOC3","ANGPTL3")))
  df$display_name <- factor(df$type, levels = lev)
  df$index <- nrow(df):1
  # 3. plotting
  is_drug <- df$category == "cis-MR"
  # manually set color order so that polygenic MR "Overall" is consistently grey
  levs <- if(is.factor(df$display_name)) levels(df$display_name) else sort(unique(as.character(df$display_name)))
  cols <- setNames(scales::hue_pal()(length(levs)), levs); cols["Overall"] <- "grey50"
  ggp<- ggplot(df, aes(x = est, y = index)) +
    geom_vline(xintercept = 0, linetype = 2) +
    # --- Polygenic MR ---
    geom_errorbarh(data = subset(df, !is_drug),aes(xmin = CIlow, xmax = CIup, color = display_name),width = 0, linewidth = 0.6) +
    geom_point(data = subset(df, !is_drug),aes(color = display_name),shape = 1, size = 3, stroke = 0.9) +
    # --- cis-MR ---
    geom_errorbarh(data = subset(df, is_drug),aes(xmin = CIlow, xmax = CIup), width = 0, linewidth = 0.6, color = "black") +
    geom_point(data = subset(df, is_drug),aes(shape = display_name),size = 3, color = "black") +
    # --- style: shape legend for cis-MR, color legend for polygenic MR ---
    scale_shape_manual(values = c(APOC3 = 15, ANGPTL3 = 17)) +
    scale_color_discrete(name = "Polygenic MR") +
    guides(shape = guide_legend(title = "cis-MR", order = 1),
           color = guide_legend(title = "Polygenic MR", order = 2)
    ) +
    labs(title = sprintf("TG → %s", outcome_name), x = NULL, y = NULL)+   # →
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))+
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),  # hide the default L-shaped axis lines from classic theme (prevents uneven line thickness)
      axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
      legend.title = element_text(size = 10),
      legend.key.height = unit(0.5, "cm") )+ scale_color_manual(name = "Polygenic MR", values = cols)

  # store the ggplot
  TGggplots[[outcome_name]]<- ggp


}

TGggplots$glszm_ZoneEntropy_inp
TGggplots$glszm_GrayLevelVariance_inp  # interesting finding: cis-MR (APOC3 -> TG -> this MRI) is negative; while polygenic TG -> this MRI is positive
TGggplots$CAD

saveRDS(TGggplots,
        "/path/to/project/targetMR/cisMR2_res/TGggplots.rds")
TGggplots <- readRDS("/path/to/project/targetMR/cisMR2_res/TGggplots.rds")


### ===========================================================================================
### Final plots ===============================================================================
### ===========================================================================================

BMIggplots <- readRDS("/path/to/project/targetMR/cisMR2_res/BMIggplots.rds")
LDLggplots <- readRDS("/path/to/project/targetMR/cisMR2_res/LDLggplots.rds")
TGggplots <- readRDS("/path/to/project/targetMR/cisMR2_res/TGggplots.rds")

###  Output figures:
# forest plot (BMI -> MRI)   ||  forest plot (LDL-C -> MRI)  forest plot (LDL-C -> CAD)   MR-cluster scatterplot
# forest plot (TG -> MRI1)       forest plot (TG -> MRI2)    forest plot (TG -> CAD)       MR-cluster scatterplot

### final combined plot
BMIggplots$firstorder_Minimum_inp | LDLggplots$glszm_ZoneEntropy_inp | LDLggplots$CAD
TGggplots$glszm_GrayLevelVariance_inp | TGggplots$glszm_ZoneEntropy_inp  | TGggplots$CAD


( BMIggplots$firstorder_Minimum_inp | LDLggplots$glszm_ZoneEntropy_inp | LDLggplots$CAD ) /
  ( TGggplots$glszm_GrayLevelVariance_inp | TGggplots$glszm_ZoneEntropy_inp | TGggplots$CAD )
# cisMR_vs_polygenicMR  # 1000 500



### Heterogeneity test: cis-MR vs. polygenic MR estimates

## TG on glszm_GrayLevelVariance_inp: polygenic MR (overall) vs. cis-MR (via APOC3)
index <- LDaware_cisMR_res$APOC3$feature_name=='glszm_GrayLevelVariance_inp'
est1<- LDaware_cisMR_res$APOC3$LD_aware_IVWest[index]
se1 <- LDaware_cisMR_res$APOC3$LD_aware_IVWse[index]

est2 <- TGres$glszm_GrayLevelVariance_inp[1,3]
se2 <- TGres$glszm_GrayLevelVariance_inp[1,4]

heterogeneity_test_statistic_value <- ( est1 - est2 )/sqrt( se1^2 + se2^2   )
2*pnorm(  -abs(heterogeneity_test_statistic_value )  )     # pvalue



## BMI on firstorder_Minimum_inp: polygenic MR (overall) vs. cis-MR (via GLP1R)
index <- LDaware_cisMR_res$GLP1R$feature_name=='firstorder_Minimum_inp'
est1<- LDaware_cisMR_res$GLP1R$LD_aware_IVWest[index]        # MR-LD-IVW
se1 <- LDaware_cisMR_res$GLP1R$LD_aware_IVWse[index]

index <- cisMR_PCA_res$GLP1R$`feature name`=='firstorder_Minimum_inp'
est1<- cisMR_PCA_res$GLP1R$est_robust[index]        # MR-PC-GMM (wider CI)
se1 <- cisMR_PCA_res$GLP1R$se_robust[index]

est2 <- BMIres$firstorder_Minimum_inp[1,3]
se2 <- BMIres$firstorder_Minimum_inp[1,4]

heterogeneity_test_statistic_value <- ( est1 - est2 )/sqrt( se1^2 + se2^2   )
2*pnorm(  -abs(heterogeneity_test_statistic_value )  )     # pvalue
