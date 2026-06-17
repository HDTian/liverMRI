# =============================================================================
# Author: Haodong Tian
# Description: LD-aware cis-MR analysis for drug target genes (GLP1R, PCSK9,
#   HMGCR, APOC3, ANGPTL3, LPA) — uses PLINK LD-clumping to select independent
#   cis-IVs and runs MR against liver MRI features.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================


### LD-aware cisMR


### STEP 3: cis-IV selection and summary data collection for the gene with good proxy -------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------------------

# We use PLINK LD-clumping to select independent (LD r^2 < 0.01), significant (p < 1e-6) SNPs (MAF > 0.01).
# Other methods exist, including COJO, PCA, and fine-mapping.
# LD-clumping is preferred here because:
# (1) It produces independent SNPs, enabling robust MR methods that assume IV independence.
# (2) Methods like SuSiE require precise in-sample LD; using an external LD reference (e.g., 1000G) can be unreliable.


# About BMI:
# We use BMI GWAS without UKB samples, because GLP1R cis-SNPs with UKB+GIANT BMI suffer from weak IV (MR-PC-GMM F-statistic < 10).
# In two-sample MR, even when F < 10, the significance results remain reliable.
# MVP EUR BMI (meanBMI; RINT) GWAS: https://www.ebi.ac.uk/gwas/studies/GCST90475156  --> BMIMVP


# PLINK configuration
plink  <- "/path/to/project/PheWAS/plink_mac_20250615/plink"
bfile  <- "/path/to/project/PheWAS/1000G_QC"
# GWAS input and PLINK output path
indir  <- "/path/to/project/targetMR/ProxyGWAS"
outdir <- file.path(indir, "LD_clumped_results2")  # PLINK results will be stored here

# Gene vector and its mapping to proxy trait
genes_with_good_proxy <- c("GLP1R","PCSK9","HMGCR","APOC3","ANGPTL3","LPA")
proxy_map <- c(GLP1R ="BMI_MVP" ,   # use BMI_MVP because UKB+GIANT BMI GWAS leads to weak IV for GLP1R;
                                     # weak IV is particularly dangerous in one-sample MR
               PCSK9 = "LDL_C", HMGCR = "LDL_C",
               APOC3 = "TG",    ANGPTL3 = "TG",
               LPA="Lp_a")  # gene-to-proxy mapping


### Loop each gene: store the gene cis-variants proxy GWAS [only once]
# (convenient for future summary-data collection if the IV selection strategy changes)

# already done

### Loop each gene: PLINK LD-clumping (loose threshold r^2 = 0.4, as suggested by Ang et al. cis-MR paper)
for(g in genes_with_good_proxy){
  trait <- proxy_map[g]

  cat("clumping:", g, "cis-variants, based on GWAS of", trait, "\n")

  # ----- read the gene cis-variants proxy GWAS  ----
  DT<- fread(file.path(indir, paste0(g, "_", trait, ".cis.tsv")))
  if(!nrow(DT)) next

  # ---- auto-detect rsID and p-value columns (with priority) ----
  rs_candidates <- c("rsID","rs_id","rsid","variant_id","ID","SNP")
  p_candidates  <- c("p_value","P-value","pvalue")
  rs_col <- intersect(rs_candidates, names(DT))[1]
  p_col  <- intersect(p_candidates,  names(DT))[1]
  if (is.na(rs_col) || is.na(p_col)) {stop("Cannot find rsID or p-value column in ", basename(f))}

  # Since GWAS p-values are sometimes stored as strings (e.g., "5e-8"), and as.numeric() truncates values
  # below 1e-300 to 0, the optimal strategy is to re-derive p-values from |Z-score|.
  est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
  se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
  if (is.na(est_col) || is.na(se_col)) {stop("Cannot find est or se column (for |z-score| rank) in ", basename(f))}
  DT[, z := get(est_col) / get(se_col)]
  DT <- DT[order(-abs(z))]  # sort by significance (note: PLINK ignores row order, but ordering is critical for the newP truncation step below)
  DT[, newP := 2 * pnorm(-abs(z))]  # re-derive p-value from Z-score for PLINK clumping
  p0 <- 1e-300  # minimum p-value recognizable by PLINK
  DT[!is.finite(newP) | newP < p0, newP := p0]  # truncate + small perturbation to preserve rank stability
  DT[newP == p0, newP := newP + seq_along(.I) * 1e-302]
  # overwrite the original p-value column with the newly derived newP
  DT[, (p_col) := newP]



  # ---- create PLINK-format GWAS file (rsID + p) ----
  tmp <- tempfile(fileext = ".txt")
  fwrite(DT[, .(SNP = get(rs_col), P = get(p_col))], tmp, sep = "\t")

  # ---- output prefix: gene + trait ----
  prefix <- file.path(outdir, paste0(g, "_", trait))
  # ---- PLINK clumping ----
  cmd <- sprintf(
    '%s --bfile %s --clump %s --clump-p1 0.001 --clump-p2 0.001 --clump-r2 0.4 --clump-kb 500 --maf 0.01 --out %s',
    shQuote(plink), shQuote(bfile), shQuote(tmp), shQuote(prefix)
  )
  system(cmd)  # Warning like 'rs12129899' is missing from the main dataset -- due to MAF filtering
}



### Check the number of lead variants per gene
clump_dir <- outdir
clumped_files <- list.files(clump_dir, pattern = "\\.clumped$", full.names = TRUE)
lead_counts <- data.table(
  file   = basename(clumped_files),
  n_lead = vapply(clumped_files, function(f){
    ln <- readLines(f)
    hdr <- grep("^\\s*CHR\\b", ln)[1]      # locate the header line
    if (is.na(hdr)) return(0L)             # no header = no results
    nrow(fread(f, skip = hdr - 1L))        # each row = one lead SNP
  }, integer(1))
)
lead_counts_proxy <- lead_counts
lead_counts_proxy
# 1:    ANGPTL3_TG.clumped     32
# 2:      APOC3_TG.clumped     65
# 3: GLP1R_BMI_MVP.clumped     22
# 4:   HMGCR_LDL_C.clumped     43      # r^2 = 0.4
# 5:      LPA_Lp_a.clumped    104
# 6:   PCSK9_LDL_C.clumped     78


# 1:    ANGPTL3_TG.clumped      5
# 2:      APOC3_TG.clumped     19
# 3: GLP1R_BMI_MVP.clumped      3       # r^2 = 0.01
# 4:   HMGCR_LDL_C.clumped      7
# 5:      LPA_Lp_a.clumped     29
# 6:   PCSK9_LDL_C.clumped     14



### STEP 6: get the final exposure data and run MR analysis (LD-aware IVW) ---------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------------------



folder_map <- c(GLP1R ="ProxyGWAS" ,       # 3
                LPA = "ProxyGWAS",          # 29
                PCSK9 = "ProxyGWAS",        # 14
                HMGCR = "ProxyGWAS",        # 7
                APOC3 = "ProxyGWAS",        # 19
                ANGPTL3 = "ProxyGWAS"       # 5
)




radiomicsGWAS <- '/path/to/project/targetMR/radiomicsGWAS'
feature_paths <- list.files(radiomicsGWAS, pattern = "\\.tsv$", full.names = TRUE)
feature_paths
# feature_paths is a vector of full paths to all liver MRI feature GWAS .tsv files

LDaware_cisMR_res<-list()

for(gene_name in gene_names){
  cat(paste0( '\n=============================','\n' ))
  cat(paste0( 'Current gene: ', gene_name ,  '\n' ))


  folder_used <- paste0('/path/to/project/targetMR/'  , folder_map[gene_name])

  ## get the exposure's clumped cis-QTLs and the exposure GWAS data (rsID EA est se) ----------------------------
  ## ------------------------------------------------------------------------------------------------------------
  ## ----  find clumped file for this gene ----
  clump_dir <- file.path(folder_used, "LD_clumped_results2")   # LD_clumped_results2 is the weak-pruning (r^2=0.4) results
  clump_file <- list.files(clump_dir, pattern = gene_name, full.names = TRUE)
  clump_file <- clump_file[grepl("\\.clumped$", clump_file)][1]
  ## ---- extract rsIDs from .clumped ----
  clumped_rs <- fread(clump_file, fill=TRUE)$SNP
  clumped_rs <- clumped_rs[nzchar(clumped_rs)]   # remove empty strings
  ## ---- read exposure GWAS (cis.tsv) ----
  exp_file <- list.files(folder_used, pattern = gene_name, full.names = TRUE)
  exp_file <- exp_file[grepl("cis\\.tsv$", exp_file)][1]
  DT <- data.table::fread(exp_file)
  ## ---- auto-detect columns (priority order) ----
  rs_col  <- intersect(c("rsID","rs_id","rsid","variant_id","ID","SNP"), names(DT))[1]
  ea_col  <- intersect(c("EA","alt","ALLELE1","effect_allele","Allele1"), names(DT))[1]
  est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
  se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
  ## ---- restrict to clumped SNPs & standardize output ----
  GXdata <- DT[get(rs_col) %in% clumped_rs, .(rsID = get(rs_col),  EA   = get(ea_col) , est  = get(est_col), se   = get(se_col) )]

  ## loop all outcomes (the 59 retained liver MRI features) -----------------------------------------------------
  ## ------------------------------------------------------------------------------------------------------------
  MR_matrix_results <- c()
  for (kk in 1:length(feature_paths)  ) {
    cat(paste0(  kk, '-' )  )

    feature_path<- feature_paths[kk]
    feature_name <- sub("\\.tsv$", "", basename(feature_path))
    feature_name

    ## get the outcome GWAS data (rsID EA est se) -----------------------------------------
    dat <- read.table(feature_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    dat <- dat[, c("SNP", "ALLELE1", "BETA", "SE")]
    GYdata <- dat[match(GXdata$rsID ,   dat$SNP)  , ]
    # should be the same nrow; as liver MRI GWAS is based on 1000G, so no missing rsID

    ## harmonization (according to the LD ref panel effect allele) ------------------------
    GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]
    bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
    by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
    ref <- bim[match(GXdata_$rsID, bim$SNP), ]
    flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1
    bx[flip_e] <- -bx[flip_e];by[flip_o] <- -by[flip_o]

    ## get the corresponding LD correlation matrix (only once, as all liver MRI outcomes share the same cis-SNPs)
    if( (kk == 1)|| !exists("LD")  || (  nrow( LD )!= nrow(GXdata_)  )  ){
      cat('-(getLDmatrix)-'  )
      rs <- GXdata_$rsID
      tf  <- tempfile();out <- tempfile();writeLines(rs, tf)
      # compute LD (Pearson r)
      system(sprintf(
        "%s --bfile %s --extract %s --r square --write-snplist --out %s",
        shQuote(plink), shQuote(bfile), shQuote(tf), shQuote(out)
      ), ignore.stdout = TRUE, ignore.stderr = TRUE)
      # read the SNP order used by PLINK
      ord <- fread(paste0(out, ".snplist"), header = FALSE)[[1]]
      # read LD matrix
      LD0 <- as.matrix(fread(paste0(out, ".ld"), header = FALSE))
      # reorder LD matrix to match rs (= GXdata_$rsID) order
      idx <- match(rs, ord)
      if (any(is.na(idx))) stop("Some rsID not found in LD reference panel")
      LD <- LD0[idx, idx]
    }
    # quick check whether or not LD is invertible
    is.finite(rcond(LD)) && rcond(LD) > 0  # TRUE: numerically invertible
    all(eigen(LD, symmetric = TRUE, only.values = TRUE)$values > 0) # all eigenvalues > 0; certainly invertible

    ## LD-aware random-effect IVW (store F-statistic and Q test results) ----------------------
    MRres<-mr_ivw(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse, correlation = LD ) )
    #mr_plot( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse )  , orientate=TRUE)

    ## robust MR methods (only work with >= 3 SNPs) --------------------------------------
    #  not used for correlated SNPs

    ## Store/update the vector/matrix results ---------------------------------------------
    MR_vector_results <- c( kk,sub("\\.cis\\.tsv$", "", basename(exp_file)), feature_name ,
                            nrow(GYdata) , MRres@Fstat , MRres@Heter.Stat[2],
                            MRres@Estimate,MRres@StdError, MRres@Pvalue)
    MR_matrix_results <- rbind( MR_matrix_results , MR_vector_results )
  }
  # finish the loop

  colnames(MR_matrix_results) <- c( 'feature ID','Exposure', 'feature_name',
                                    'num_of_snps','Fstatistic','Qpvalue',
                                    'LD_aware_IVWest','LD_aware_IVWse','LD_aware_IVWp')
  MR_matrix_results <- as.data.frame(MR_matrix_results)
  MR_matrix_results[-(2:3)] <- lapply(MR_matrix_results[-(2:3)], as.numeric)

  MR_matrix_results$FDRp <- p.adjust(MR_matrix_results$LD_aware_IVWp, method = "BH")

  ### store all MR results for this gene
  LDaware_cisMR_res[[gene_name]] <- MR_matrix_results

}


View(LDaware_cisMR_res$APOC3  )
View(LDaware_cisMR_res$HMGCR  )
View(LDaware_cisMR_res$PCSK9  )
View(LDaware_cisMR_res$GLP1R  )
View(cisMR_res$GLP1R  )

saveRDS(LDaware_cisMR_res,
        "/path/to/project/targetMR/cisMR_res/LDaware_cisMR_res.rds")
LDaware_cisMR_res <- readRDS("/path/to/project/targetMR/cisMR_res/LDaware_cisMR_res.rds")


### QQ plots for [MR-PC-GMM] [LD-aware-MR] [MR-IVW]
# QQ plot of LD-aware MR
qqman::qq(LDaware_cisMR_res$GLP1R$LD_aware_IVWp, main = "Q-Q Plot")
qqman::qq(LDaware_cisMR_res$APOC3$LD_aware_IVWp, main = "Q-Q Plot")

# QQ plot of MR-PC-GMM
qqman::qq(cisMR_PCA_res$GLP1R$pvalue_robust, main = "Q-Q Plot")
qqman::qq(cisMR_PCA_res$HMGCR$pvalue_robust, main = "Q-Q Plot")
qqman::qq(cisMR_PCA_res$PCSK9$pvalue_robust, main = "Q-Q Plot")
qqman::qq(cisMR_PCA_res$ANGPTL3$pvalue_robust, main = "Q-Q Plot")
qqman::qq(cisMR_PCA_res$APOC3$pvalue_robust, main = "Q-Q Plot")
# QQ plot of MR-IVW
qqman::qq(cisMR_res$HMGCR$IVWp, main = "Q-Q Plot")

# Note: Q-Q inflation here is expected — the 59 liver MRI outcomes are highly correlated,
# so the test statistics are not i.i.d. under the global null.
# Lambda (genomic inflation) computed from correlated outcomes should not be interpreted as true inflation.


### MR-LD-IVW vs. MR-IVW scatter plot

genes <- intersect(names(LDaware_cisMR_res), names(cisMR_res))
par(mfrow=c(2,3))
for (g in genes) {
  x  <- LDaware_cisMR_res[[g]]$LD_aware_IVWest
  y  <- cisMR_res[[g]]$IVWest
  xs <- LDaware_cisMR_res[[g]]$LD_aware_IVWse
  ys <- cisMR_res[[g]]$IVWse

  plot(x, y, pch = 16, cex = 1.2, col = "black",main=g,
       xlab = "MR-LD-IVW estimates", ylab = "MR-IVW estimates")
  abline(0,1,col='blue')
}
par(mfrow=c(1,1))
# MRLDIVWvsMRIVW  # 650 500




### Step 8: Positive control: each gene -> exposure -> CAD or T2D ----------------------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------------------


positive_outcome_map<-c(GLP1R ="T2D" ,
                        LPA = "CAD",PCSK9 = "CAD",HMGCR = "CAD",APOC3 = "CAD",ANGPTL3 = "CAD" )


Postivie_LDaware_cisMR_res<-c()
for(gene_name in gene_names){
  cat(paste0( '\n=============================','\n' ))
  cat(paste0( 'Current gene: ', gene_name ,  '\n' ))

  folder_used <- paste0('/path/to/project/targetMR/'  , folder_map[gene_name])

  ## get the exposure's clumped cis-QTLs and the exposure GWAS data (rsID EA est se) ----------------------------
  ## ------------------------------------------------------------------------------------------------------------
  ## ----  find clumped file for this gene ----
  clump_dir <- file.path(folder_used, "LD_clumped_results2")   # LD_clumped_results2 is the weak-pruning (r^2=0.4) results
  clump_file <- list.files(clump_dir, pattern = gene_name, full.names = TRUE)
  clump_file <- clump_file[grepl("\\.clumped$", clump_file)][1]
  ## ---- extract rsIDs from .clumped ----
  clumped_rs <- fread(clump_file, fill=TRUE)$SNP
  clumped_rs <- clumped_rs[nzchar(clumped_rs)]   # remove empty strings
  ## ---- read exposure GWAS (cis.tsv) ----
  exp_file <- list.files(folder_used, pattern = gene_name, full.names = TRUE)
  exp_file <- exp_file[grepl("cis\\.tsv$", exp_file)][1]
  DT <- data.table::fread(exp_file)
  ## ---- auto-detect columns (priority order) ----
  rs_col  <- intersect(c("rsID","rs_id","rsid","variant_id","ID","SNP"), names(DT))[1]
  ea_col  <- intersect(c("EA","alt","ALLELE1","effect_allele","Allele1"), names(DT))[1]
  est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
  se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
  ## ---- restrict to clumped SNPs & standardize output ----
  GXdata <- DT[get(rs_col) %in% clumped_rs, .(rsID = get(rs_col),  EA   = get(ea_col) , est  = get(est_col), se   = get(se_col) )]


  ## Positive control outcome GWAS
  ## get the outcome GWAS data (rsID EA est se) -----------------------------------------
  outfile<-paste0('/path/to/project/targetMR/otherGWAS/',
                  positive_outcome_map[gene_name],'.cis.tsv'    )
  dat <- read.table(outfile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  cat(  paste0(   'positive outcome is: ' , positive_outcome_map[gene_name] , '\n'  ) )
  if( positive_outcome_map[gene_name] == 'T2D'  ){
    dat<-dat[, c("rsID", "effect_allele", "Fixed.effects_beta", "Fixed.effects_SE","effect_allele_frequency")]
    dat$N <- 933970   # total (raw) sample size from published GWAS dictionary
  }else{
    dat <- dat[, c("rsID", "effect_allele", "beta", "standard_error","effect_allele_frequency","n")]
  }
  names(dat) <- c("SNP", "ALLELE1", "BETA", "SE",'EAF','N')
  GYdata <- dat[match(GXdata$rsID ,   dat$SNP)  , ]  # NA rows are preserved by match; removed below
  GYdata <- na.omit(GYdata)   # remove NA rows

  ## harmonization (according to the LD ref panel effect allele) ------------------------
  GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]
  bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
  by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
  ref <- bim[match(GXdata_$rsID, bim$SNP), ]
  flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1
  bx[flip_e] <- -bx[flip_e];by[flip_o] <- -by[flip_o]

  ## get the corresponding LD correlation matrix (only once, as all liver MRI outcomes share the same cis-SNPs)
  if( (kk == 1)|| !exists("LD")  || (  nrow( LD )!= nrow(GXdata_)  )  ){
    cat('-(getLDmatrix)-'  )
    rs <- GXdata_$rsID
    tf  <- tempfile();out <- tempfile();writeLines(rs, tf)
    # compute LD (Pearson r)
    system(sprintf(
      "%s --bfile %s --extract %s --r square --write-snplist --out %s",
      shQuote(plink), shQuote(bfile), shQuote(tf), shQuote(out)
    ), ignore.stdout = TRUE, ignore.stderr = TRUE)
    # read the SNP order used by PLINK
    ord <- fread(paste0(out, ".snplist"), header = FALSE)[[1]]
    # read LD matrix
    LD0 <- as.matrix(fread(paste0(out, ".ld"), header = FALSE))
    # reorder LD matrix to match rs (= GXdata_$rsID) order
    idx <- match(rs, ord)
    if (any(is.na(idx))) stop("Some rsID not found in LD reference panel")
    LD <- LD0[idx, idx]
  }

  ## LD-aware random-effect IVW (store F-statistic and Q test results) ----------------------
  MRres<-mr_ivw(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse, correlation = LD ) )
  #mr_plot( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse )  , orientate=TRUE)

  ## robust MR methods (only work with >= 3 SNPs) --------------------------------------
  #  not used for correlated SNPs

  ## Store/update the vector/matrix results ---------------------------------------------
  MR_vector_results <- c( gene_name ,  folder_map[gene_name], positive_outcome_map[gene_name],
                          nrow(GYdata) , MRres@Fstat , MRres@Heter.Stat[2],
                          MRres@Estimate,MRres@StdError, MRres@CILower, MRres@CIUpper, MRres@Pvalue)


  ## Store the results
  Postivie_LDaware_cisMR_res<-rbind( Postivie_LDaware_cisMR_res , MR_vector_results )
}

colnames(Postivie_LDaware_cisMR_res) <- c( 'gene','Exposure','Outcome',
                                           'num_of_SNPs', 'Fstatistic','Qpvalue',
                                   'Estimate','se'  ,'CI_low','CI_up','pvalue')
Postivie_LDaware_cisMR_res <- as.data.frame(Postivie_LDaware_cisMR_res)
Postivie_LDaware_cisMR_res[-(1:3)] <- lapply(Postivie_LDaware_cisMR_res[-(1:3)], as.numeric)

Postivie_LDaware_cisMR_res
#    gene  Exposure Outcome num_of_SNPs Fstatistic       Qpvalue  Estimate         se      CI_low     CI_up       pvalue
#   GLP1R ProxyGWAS     T2D          45   10.98198  9.582338e-05 0.3034141 0.13643190  0.03601251 0.5708157 2.615359e-02
#     LPA ProxyGWAS     CAD         107 1449.38377 3.402038e-264 0.1626645 0.01962737  0.12419555 0.2011334 1.155222e-16
#   PCSK9 ProxyGWAS     CAD          94  132.07131  1.063004e-03 0.6102551 0.04250370  0.52694938 0.6935608 9.534600e-47
#   HMGCR ProxyGWAS     CAD          52   62.03107  5.972660e-02 0.3787397 0.07016227  0.24122417 0.5162552 6.736767e-08
#   APOC3 ProxyGWAS     CAD          69  417.36663  2.294388e-07 0.1928620 0.02806685  0.13785198 0.2478720 6.352042e-12
# ANGPTL3 ProxyGWAS     CAD          39   73.13553  8.425283e-02 0.1251874 0.07100610 -0.01398204 0.2643568 7.789194e-02

saveRDS(Postivie_LDaware_cisMR_res,
        "/path/to/project/targetMR/cisMR_res/Positive_LDaware_cisMR_res.rds")
Postivie_LDaware_cisMR_res <- readRDS("/path/to/project/targetMR/cisMR_res/Positive_LDaware_cisMR_res.rds")




### Negative control ------------------
negative_outcome_map<-c(GLP1R ="Tanning" ,
                        LPA = "Tanning",PCSK9 = "Tanning",HMGCR = "Tanning",APOC3 = "Tanning",ANGPTL3 = "Tanning" )


Negative_LDaware_cisMR_res<-c()
for(gene_name in gene_names){
  cat(paste0( '\n=============================','\n' ))
  cat(paste0( 'Current gene: ', gene_name ,  '\n' ))

  folder_used <- paste0('/path/to/project/targetMR/'  , folder_map[gene_name])

  ## get the exposure's clumped cis-QTLs and the exposure GWAS data (rsID EA est se) ----------------------------
  ## ------------------------------------------------------------------------------------------------------------
  ## ----  find clumped file for this gene ----
  clump_dir <- file.path(folder_used, "LD_clumped_results2")   # LD_clumped_results2 is the weak-pruning (r^2=0.4) results
  clump_file <- list.files(clump_dir, pattern = gene_name, full.names = TRUE)
  clump_file <- clump_file[grepl("\\.clumped$", clump_file)][1]
  ## ---- extract rsIDs from .clumped ----
  clumped_rs <- fread(clump_file, fill=TRUE)$SNP
  clumped_rs <- clumped_rs[nzchar(clumped_rs)]   # remove empty strings
  ## ---- read exposure GWAS (cis.tsv) ----
  exp_file <- list.files(folder_used, pattern = gene_name, full.names = TRUE)
  exp_file <- exp_file[grepl("cis\\.tsv$", exp_file)][1]
  DT <- data.table::fread(exp_file)
  ## ---- auto-detect columns (priority order) ----
  rs_col  <- intersect(c("rsID","rs_id","rsid","variant_id","ID","SNP"), names(DT))[1]
  ea_col  <- intersect(c("EA","alt","ALLELE1","effect_allele","Allele1"), names(DT))[1]
  est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
  se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
  ## ---- restrict to clumped SNPs & standardize output ----
  GXdata <- DT[get(rs_col) %in% clumped_rs, .(rsID = get(rs_col),  EA   = get(ea_col) , est  = get(est_col), se   = get(se_col) )]


  ## Negative control outcome GWAS
  ## get the outcome GWAS data (rsID EA est se) -----------------------------------------
  outfile<-paste0('/path/to/project/targetMR/otherGWAS/',
                  negative_outcome_map[gene_name],'.cis.tsv'    )
  dat <- read.table(outfile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  cat(  paste0(   'negative outcome is: ' , negative_outcome_map[gene_name] , '\n'  ) )

  dat <- dat[, c("rsID", "effect_allele", "beta", "se","minor_AF","n_complete_samples")]  # note: MAF != EAF but does not affect harmonization here

  names(dat) <- c("SNP", "ALLELE1", "BETA", "SE",'EAF','N')
  GYdata <- dat[match(GXdata$rsID ,   dat$SNP)  , ]
  GYdata <- na.omit(GYdata)   # remove NA rows

  ## harmonization (according to the LD ref panel effect allele) ------------------------
  GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]
  bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
  by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
  ref <- bim[match(GXdata_$rsID, bim$SNP), ]
  flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1
  bx[flip_e] <- -bx[flip_e];by[flip_o] <- -by[flip_o]

  ## get the corresponding LD correlation matrix (only once, as all liver MRI outcomes share the same cis-SNPs)
  if( (kk == 1)|| !exists("LD")  || (  nrow( LD )!= nrow(GXdata_)  )  ){
    cat('-(getLDmatrix)-'  )
    rs <- GXdata_$rsID
    tf  <- tempfile();out <- tempfile();writeLines(rs, tf)
    # compute LD (Pearson r)
    system(sprintf(
      "%s --bfile %s --extract %s --r square --write-snplist --out %s",
      shQuote(plink), shQuote(bfile), shQuote(tf), shQuote(out)
    ), ignore.stdout = TRUE, ignore.stderr = TRUE)
    # read the SNP order used by PLINK
    ord <- fread(paste0(out, ".snplist"), header = FALSE)[[1]]
    # read LD matrix
    LD0 <- as.matrix(fread(paste0(out, ".ld"), header = FALSE))
    # reorder LD matrix to match rs (= GXdata_$rsID) order
    idx <- match(rs, ord)
    if (any(is.na(idx))) stop("Some rsID not found in LD reference panel")
    LD <- LD0[idx, idx]
  }

  ## LD-aware random-effect IVW (store F-statistic and Q test results) ----------------------
  MRres<-mr_ivw(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse, correlation = LD ) )
  #mr_plot( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse )  , orientate=TRUE)

  ## robust MR methods (only work with >= 3 SNPs) --------------------------------------
  #  not used for correlated SNPs

  ## Store/update the vector/matrix results ---------------------------------------------
  MR_vector_results <- c( gene_name ,  folder_map[gene_name], negative_outcome_map[gene_name],
                          nrow(GYdata) , MRres@Fstat , MRres@Heter.Stat[2],
                          MRres@Estimate,MRres@StdError, MRres@CILower, MRres@CIUpper, MRres@Pvalue)


  ## Store the results
  Negative_LDaware_cisMR_res<-rbind( Negative_LDaware_cisMR_res , MR_vector_results )
}

colnames(Negative_LDaware_cisMR_res) <- c( 'gene','Exposure','Outcome',
                                           'num_of_SNPs', 'Fstatistic','Qpvalue',
                                           'Estimate','se'  ,'CI_low','CI_up','pvalue')
Negative_LDaware_cisMR_res <- as.data.frame(Negative_LDaware_cisMR_res)
Negative_LDaware_cisMR_res[-(1:3)] <- lapply(Negative_LDaware_cisMR_res[-(1:3)], as.numeric)

Negative_LDaware_cisMR_res
#     gene  Exposure Outcome num_of_SNPs Fstatistic      Qpvalue     Estimate          se       CI_low       CI_up     pvalue
#    GLP1R ProxyGWAS Tanning          48   11.08942 7.804742e-01 -0.025290518 0.033433896 -0.090819749 0.040238714 0.44938935
#      LPA ProxyGWAS Tanning         106 1495.07969 8.321798e-05 -0.005338358 0.003344231 -0.011892929 0.001216214 0.11042436
#    PCSK9 ProxyGWAS Tanning          90  137.26048 6.989706e-03 -0.011258520 0.017767844 -0.046082855 0.023565814 0.52631199
#    HMGCR ProxyGWAS Tanning          51   62.47248 1.396082e-01 -0.005982923 0.030899697 -0.066545216 0.054579371 0.84647031
#    APOC3 ProxyGWAS Tanning          68  423.38376 2.530049e-01  0.005593602 0.011092754 -0.016147796 0.027335000 0.61408060
#  ANGPTL3 ProxyGWAS Tanning          38   74.83694 1.007035e-01  0.067888176 0.034968233 -0.000648302 0.136424654 0.05220684



saveRDS(Negative_LDaware_cisMR_res,
        "/path/to/project/targetMR/cisMR_res/Negative_LDaware_cisMR_res.rds")
Negative_LDaware_cisMR_res <- readRDS("/path/to/project/targetMR/cisMR_res/Negative_LDaware_cisMR_res.rds")




### Supplementary Figure: scatter plot of MR-LD-IVW for GLP1R cis-SNPs on BMI and CAD
gene_name <- 'GLP1R'


folder_used <-"/path/to/project/targetMR/ProxyGWAS"  # BMI_MVP

## get the exposure's clumped cis-QTLs and the exposure GWAS data (rsID EA est se) ----------------------------
## ------------------------------------------------------------------------------------------------------------
## ----  find clumped file for this gene ----
clump_dir <- file.path(folder_used, "LD_clumped_results2")   # LD_clumped_results2 is the weak-pruning (r^2=0.4) results
clump_file <- list.files(clump_dir, pattern = gene_name, full.names = TRUE)
clump_file <- clump_file[grepl("\\.clumped$", clump_file)][1]
## ---- extract rsIDs from .clumped ----
clumped_rs <- fread(clump_file, fill=TRUE)$SNP
clumped_rs <- clumped_rs[nzchar(clumped_rs)]   # remove empty strings
## ---- read exposure GWAS (cis.tsv) ----
exp_file <- list.files(folder_used, pattern = gene_name, full.names = TRUE)
exp_file <- exp_file[grepl("cis\\.tsv$", exp_file)][1]
DT <- data.table::fread(exp_file)
## ---- auto-detect columns (priority order) ----
rs_col  <- intersect(c("rsID","rs_id","rsid","variant_id","ID","SNP"), names(DT))[1]
ea_col  <- intersect(c("EA","alt","ALLELE1","effect_allele","Allele1"), names(DT))[1]
est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
## ---- restrict to clumped SNPs & standardize output ----
GXdata <- DT[get(rs_col) %in% clumped_rs, .(rsID = get(rs_col),  EA   = get(ea_col) , est  = get(est_col), se   = get(se_col) )]


## Positive control outcome GWAS (CAD)
## get the outcome GWAS data (rsID EA est se) -----------------------------------------
outfile<-paste0('/path/to/project/targetMR/otherGWAS/',
                'CAD','.cis.tsv'    )
dat <- read.table(outfile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
cat(  paste0(   'positive outcome is: ' , positive_outcome_map[gene_name] , '\n'  ) )

dat <- dat[, c("rsID", "effect_allele", "beta", "standard_error","effect_allele_frequency","n")]

names(dat) <- c("SNP", "ALLELE1", "BETA", "SE",'EAF','N')
GYdata <- dat[match(GXdata$rsID ,   dat$SNP)  , ]
GYdata <- na.omit(GYdata)   # remove NA rows

## harmonization (according to the LD ref panel effect allele) ------------------------
GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]
bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
ref <- bim[match(GXdata_$rsID, bim$SNP), ]
flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1
bx[flip_e] <- -bx[flip_e];by[flip_o] <- -by[flip_o]

## get the corresponding LD correlation matrix
if( (kk == 1)|| !exists("LD")  || (  nrow( LD )!= nrow(GXdata_)  )  ){
  cat('-(getLDmatrix)-'  )
  rs <- GXdata_$rsID
  tf  <- tempfile();out <- tempfile();writeLines(rs, tf)
  # compute LD (Pearson r)
  system(sprintf(
    "%s --bfile %s --extract %s --r square --write-snplist --out %s",
    shQuote(plink), shQuote(bfile), shQuote(tf), shQuote(out)
  ), ignore.stdout = TRUE, ignore.stderr = TRUE)
  # read the SNP order used by PLINK
  ord <- fread(paste0(out, ".snplist"), header = FALSE)[[1]]
  # read LD matrix
  LD0 <- as.matrix(fread(paste0(out, ".ld"), header = FALSE))
  # reorder LD matrix to match rs (= GXdata_$rsID) order
  idx <- match(rs, ord)
  if (any(is.na(idx))) stop("Some rsID not found in LD reference panel")
  LD <- LD0[idx, idx]
}

## LD-aware random-effect IVW (store F-statistic and Q test results) ----------------------
MRres<-mr_ivw(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse, correlation = LD ) )
MRres # est = -0.114   se = 0.084  CI = (-0.279, 0.052)  pvalue = 0.178

mr_plot( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse )  , orientate=TRUE)


## draw the final ggplot
library(ggplot2)

d <- data.frame(bx, bxse, by, byse)
s <- sign(d$bx); s[s==0] <- 1
d$bx <- d$bx*s; d$by <- d$by*s

ggplot(d, aes(bx, by)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = by - 1.96*byse, ymax = by + 1.96*byse),
                width = 0, colour = "grey60") +
  geom_errorbarh(aes(xmin = bx - 1.96*bxse, xmax = bx + 1.96*bxse),
                 height = 0, colour = "grey60") +
  geom_point(colour = "black") +
  labs(x = "Genetic association with BMI [SD]", y = "Genetic association with CAD [logOR]") +
  theme_classic()
