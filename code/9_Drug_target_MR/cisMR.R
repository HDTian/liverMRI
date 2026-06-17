# =============================================================================
# Author: Haodong Tian
# Description: Cis-Mendelian Randomization analysis using cis-pQTLs as
#   instruments for liver-relevant protein targets, with multiple MR methods
#   (IVW, weighted median, MR-Egger).
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================


### cis-MR analysis

# STEP 0: prepare the gene information
# STEP 1: decide and prepare the exposure data used
# STEP 2: get the hg19 cis-region rsID for each gene
# STEP 3: cis-IV selection and summary data collection for the gene with good proxy
# STEP 4: cis-IV selection and summary data collection for the gene with good protein
# STEP 5: cis-IV selection and summary data collection for the gene with good expression
# STEP 6: get the final exposure data and run MR analysis
# STEP 7: cisMR with correlated SNPs
# STEP 8: combine the results and visualization




library(data.table)
library(biomaRt)
library(rtracklayer)
library(susieR)   # fine-mapping
library(Rfast)    # fast SuSiE in susieR
library(MendelianRandomization)  # MR IVW
library(mr.raps)  # MR-RAPS
library(dplyr)
library(forestplot)
library(grid)




### STEP 0: prepare the gene information ----------------------------------------------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------------------
gene_names <- c( 'GLP1R',  'LPA','PCSK9','HMGCR', 'APOC3', 'ANGPTL3'  )

# Connect to Ensembl human gene dataset
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "asia")
# Query: gene symbol -> Ensembl ID + chromosome + position
gene_result <- getBM(attributes = c('external_gene_name',
                                    'ensembl_gene_id',
                                    'chromosome_name',
                                    'start_position',  'end_position', 'gene_biotype'),
                     filters = 'external_gene_name', values = gene_names, mart = ensembl)
gene_result <- gene_result[match(gene_names, gene_result$external_gene_name), ]

gene_result  # note: start_position and end_position refer to hg38 transcript (including UTRs) coordinates




### STEP 1: decide and prepare the exposure data used ---------------------------------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------------------

# What is the purpose?
# Test the sharp null hypothesis that the drug target has no causal effect on the outcome.
# We are able to use proxy exposure (i.e., traits in the downstream of the drug target).

# Which exposure to choose: tissue-specific gene expression (eQTL; e.g., GTEx) vs.
# plasma protein (pQTL; e.g., UKB-PPP) vs. proxy trait?
# Choose the one yielding the most independent SNPs -> this increases power and enables robust MR methods.

# GWAS resources:
# Broad Server pQTL: /path/to/server/UK_Biobank/baskets/ukb_basket_pqtl
# UKB-PPP: https://www.synapse.org/Synapse:syn51364943/files/
# deCODE pQTL: https://www.decode.com/summarydata/  (Ferkingstad et al.; note: deCODE provides genome-wide, not cis-only, regions)
# GTEx v10 eQTL: /path/to/server/project/GTEx/v10/gtex_v10_eqtl

# Minimum GWAS requirements:
# ideally: [rsID CHR BP(hg19) BP(hg38)] [EA nonEA EAF] [est se pvalue] [N]
# minimally: [rsID or CHR+BP] [EA] [est se]
# If only CHR+BP is available, we prefer to annotate rsID via liftover and the 1000G reference panel.




### STEP 2: get the hg19 cis-region rsID for each gene --------------------------------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------------------

## Strategy
# MR analysis uses only rsID, not the original GWAS chr:pos.
# Since we use 1000G as the LD reference panel (which defines all post-QC SNPs),
# we only need to identify rsIDs in the cis-region (+-100 kb, hg19) for each gene.
# This is also preferred when the proxy GWAS lacks chr:pos or when positions are unreliable (e.g., hg19/hg38 confusion).

# Since 1000G is hg19-based, we need hg19 gene positions.
ensembl37 <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 37)
gene_result_hg19 <- getBM(
  attributes = c("external_gene_name","ensembl_gene_id", "chromosome_name","start_position","end_position", "gene_biotype"),
  filters = "external_gene_name", values = gene_names, mart = ensembl37 )
gene_result_hg19 <- gene_result_hg19[match(gene_names, gene_result_hg19$external_gene_name), ]
gene_result_hg19 # hg19
gene_result; gene_result_hg38 <- gene_result     # hg38

## About the cis-window:
# +-100 kb is most common; +-250 kb is also acceptable.

## get the hg19 cis-region (+-100 kb) rsIDs for each gene
bfile  <- "/path/to/project/PheWAS/1000G_QC"  # 1000G reference panel path
bim <- fread(paste0(bfile, ".bim"),col.names = c("CHR","SNP","CM","BP","A1","A2"))

# ensure chromosome types are consistent
setDT(gene_result_hg19)
setDT(bim)
gene_result_hg19[, chromosome_name := as.character(chromosome_name)]
bim[, CHR := as.character(CHR)]

# +-100 kb (1e5) cis-SNPs (list-column)
cis_snps <- gene_result_hg19[ , .(rsID =
  list(bim[ CHR == chromosome_name & BP >= pmax(start_position-1e5, 1L) & BP <= end_position+1e5, SNP ] )  ),
  by = external_gene_name]  # ~ 1 second
cis_snps[, .(n_snps = lengths(rsID)), by = external_gene_name]  # quick check number of rsIDs included
cis_snps




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
outdir <- file.path(indir, "LD_clumped_results")  # PLINK results will be stored here

# Gene vector and its mapping to proxy trait
genes_with_good_proxy <- c("GLP1R","PCSK9","HMGCR","APOC3","ANGPTL3","LPA")
proxy_map <- c(GLP1R ="BMI_MVP" ,   # use BMI_MVP because UKB+GIANT BMI GWAS leads to weak IV for GLP1R;
                                     # weak IV is particularly dangerous in one-sample MR
               PCSK9 = "LDL_C", HMGCR = "LDL_C",
               APOC3 = "TG",    ANGPTL3 = "TG",
               LPA="Lp_a")  # gene-to-proxy mapping


### Loop each gene: store the gene cis-variants proxy GWAS [only once]
# (convenient for future summary-data collection if the IV selection strategy changes)
for(g in genes_with_good_proxy){
  trait <- proxy_map[g]
  f <- file.path(indir, paste0(trait, ".tsv"))
  if(!file.exists(f)) {print('no proxy trait GWAS exists for this gene');next}

  cat("cis_variants GWAS:", g, "cis-variants, based on GWAS of", trait, "\n")

  out_f <- file.path(indir, paste0(g, "_", trait, ".cis.tsv"))
  if(file.exists(out_f)) {print('already have cis_variants GWAS; skip'); next}

  DT <- fread(f)
  if(!nrow(DT)) next

  # ---- auto-detect rsID and p-value columns (with priority) ----
  rs_candidates <- c("rsID","rs_id","rsid","variant_id","ID","SNP")
  p_candidates  <- c("p_value","P-value","pvalue")
  rs_col <- intersect(rs_candidates, names(DT))[1]
  p_col  <- intersect(p_candidates,  names(DT))[1]
  if (is.na(rs_col) || is.na(p_col)) {stop("Cannot find rsID or p-value column in ", basename(f))}

  # ---- cis rsIDs for this gene ----
  cis_rs <- cis_snps[external_gene_name == g, rsID[[1]]]
  if(length(cis_rs) == 0) next
  # ---- restrict GWAS to cis-variants only ----
  DT <- DT[get(rs_col) %in% cis_rs]
  if(!nrow(DT)) next
  # ----- store the gene cis-variants proxy GWAS ----
  fwrite(DT,file.path(indir, paste0(g, "_", trait, ".cis.tsv")),sep = "\t")
}



### Loop each gene: PLINK LD-clumping
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
  # ---- PLINK clumping (r^2 < 0.01, p < 1e-6, MAF > 0.01) ----
  cmd <- sprintf(
    '%s --bfile %s --clump %s --clump-p1 1e-6 --clump-p2 1e-4 --clump-r2 0.01 --clump-kb 500 --maf 0.01 --out %s',
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
# 1:    ANGPTL3_TG.clumped      6
# 2:      APOC3_TG.clumped     17
# 3: GLP1R_BMI_MVP.clumped      3
# 4:   HMGCR_LDL_C.clumped      5
# 5:      LPA_Lp_a.clumped     20
# 6:   PCSK9_LDL_C.clumped     13

# 1:    ANGPTL3_TG.clumped      5
# 2:      APOC3_TG.clumped     19
# 3: GLP1R_BMI_MVP.clumped      3
# 4:   HMGCR_LDL_C.clumped      7
# 5:      LPA_Lp_a.clumped     29
# 6:   PCSK9_LDL_C.clumped     14
# Note: the number of lead SNPs can change when the Z-score-based non-truncated p-values are used,
# because the lead SNP identity changes, which in turn alters the LD clumping logic.




### STEP 4: cis-IV selection and summary data collection for the gene with good protein ----------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------------------


### About the UKB-PPP plasma protein GWAS:
# UKB-PPP is fitted by REGENIE, with hg38 positions; we use liftover to convert to hg19 and annotate rsIDs using our 1000G panel.
# UKB-PPP does not include p-values directly; we compute them from est and se.
# UKB-PPP GWAS is split by chromosome.
# Workflow: know the chr -> prepare GENE.gz -> run bash liftover_rsID_pvalue.sh -> GENE.gz.tsv (with rsID and p-value)


# PLINK configuration
plink  <- "/path/to/project/PheWAS/plink_mac_20250615/plink"
bfile  <- "/path/to/project/PheWAS/1000G_QC"
# GWAS input and PLINK output path
indir  <- "/path/to/project/targetMR/ProteinGWAS"
outdir <- file.path(indir, "LD_clumped_results")  # PLINK results will be stored here

genes_with_good_protein <- c("GLP1R","LPA","PCSK9","ANGPTL3")
# No mapping needed here as each gene maps to its uniquely named protein GWAS file (GENE.gz.tsv)



### Loop each gene: store the gene cis-variants protein GWAS [only once]
# (convenient for future summary-data collection if the IV selection strategy changes)
for(g in genes_with_good_protein){
  cat("UKB-PPP protein GWAS for :", g, "\n")
  f <- file.path(indir, paste0(g, ".gz.tsv"))
  if(!file.exists(f)) {print('no protein GWAS exists for this gene');next}

  out_f <- file.path(indir, paste0(g, ".", 'cis', ".tsv"))
  if(file.exists(out_f)) {print('already have cis_variants GWAS; skip'); next}

  DT <- fread(f)
  if(!nrow(DT)) next

  # ---- auto-detect rsID and p-value columns (with priority) ----
  rs_candidates <- c("rsID","variant_id","rsid","rs_id","ID","SNP")
  p_candidates  <- c("p_value","P-value","pvalue")
  rs_col <- intersect(rs_candidates, names(DT))[1]
  p_col  <- intersect(p_candidates,  names(DT))[1]
  if (is.na(rs_col) || is.na(p_col)) {stop("Cannot find rsID or p-value column in ", basename(f))}

  # ---- cis rsIDs for this gene ----
  cis_rs <- cis_snps[external_gene_name == g, rsID[[1]]]
  if(length(cis_rs) == 0) next
  # ---- restrict GWAS to cis-variants only ----
  DT <- DT[get(rs_col) %in% cis_rs]
  if(!nrow(DT)) next
  # ----- store the gene cis-variants protein GWAS ----
  fwrite(DT,file.path(indir, paste0(g, ".", 'cis', ".tsv")),sep = "\t")
}


### Loop each gene: PLINK LD-clumping
for(g in genes_with_good_protein){

  cat("clumping:", g, "cis-variants, based on its UKB-PPP plasma protein GWAS",  "\n")

  # ----- read the gene cis-variants protein GWAS  ----
  DT<- fread(file.path(indir, paste0(g, ".", 'cis', ".tsv")))
  if(!nrow(DT)) next

  # ---- auto-detect rsID and p-value columns (with priority) ----
  rs_candidates <- c("rsID","variant_id","rsid","rs_id","ID","SNP")
  p_candidates  <- c("p_value","P-value","pvalue")
  rs_col <- intersect(rs_candidates, names(DT))[1]
  p_col  <- intersect(p_candidates,  names(DT))[1]
  if (is.na(rs_col) || is.na(p_col)) {stop("Cannot find rsID or p-value column in ", basename(f))}

  # Re-derive p-values from |Z-score| to avoid truncation to 0 for very small p-values
  est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
  se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
  if (is.na(est_col) || is.na(se_col)) {stop("Cannot find est or se column (for |z-score| rank) in ", basename(f))}
  DT[, z := get(est_col) / get(se_col)]
  DT <- DT[order(-abs(z))]
  DT[, newP := 2 * pnorm(-abs(z))]
  p0 <- 1e-300
  DT[!is.finite(newP) | newP < p0, newP := p0]
  DT[newP == p0, newP := newP + seq_along(.I) * 1e-302]
  DT[, (p_col) := newP]


  # ---- create PLINK-format GWAS file (rsID + p) ----
  tmp <- tempfile(fileext = ".txt")
  fwrite(DT[, .(SNP = get(rs_col), P = get(p_col))], tmp, sep = "\t")
  # ---- output prefix: gene ----
  prefix <- file.path(outdir, paste0(g, ".", 'cis'))
  # ---- PLINK clumping (r^2 < 0.01, p < 1e-6, MAF > 0.01) ----
  cmd <- sprintf(
    '%s --bfile %s --clump %s --clump-p1 1e-6 --clump-p2 1e-4 --clump-r2 0.01 --clump-kb 500 --maf 0.01 --out %s',
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
lead_counts_protein <- lead_counts
lead_counts_protein
# ANGPTL3.cis.clumped      2
#     LPA.cis.clumped     21
#   PCSK9.cis.clumped      6
lead_counts_proxy




### STEP 5: cis-IV selection and summary data collection for the gene with good expression --------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------------------

### About the GTEx expression data:
# GTEx GWAS is based on hg38 without rsIDs.
# GTEx eQTL effect allele is the ALT allele (not the minor allele) -- important for harmonization in MR.
# Workflow: know the gene ID and suitable tissue -> prepare the expression GWAS ->
#   run bash liftover_rsID_n.sh -> GENE_Tissue.tsv (with rsID and sample size n)

# GTEx v10 eQTL: /path/to/server/project/GTEx/v10/gtex_v10_eqtl
# Columns: chr pos ref alt build gene_id variant_id tss_distance af ma_samples ma_count pval_nominal slope slope_se idx rsID
# Effective N: ma_count / af

# For GLP1R, tissue choices:
#   Pancreas: HbA1c/T2D pathway
#   Brain_Hypothalamus: BMI/obesity pathway



# PLINK configuration
plink  <- "/path/to/project/PheWAS/plink_mac_20250615/plink"
bfile  <- "/path/to/project/PheWAS/1000G_QC"
# GWAS input and PLINK output path
indir  <- "/path/to/project/targetMR/ExpressionGWAS"
outdir <- file.path(indir, "LD_clumped_results")  # PLINK results will be stored here

genes_with_good_expression <- c("GLP1R","LPA","APOC3","PCSK9","HMGCR","ANGPTL3")
tissue_map <- c(GLP1R ="Brain_Hypothalamus" ,  # GLP1R ="Pancreas" (alternative tissue)
                LPA = "Liver", PCSK9 = "Liver", HMGCR = "Liver", APOC3 = "Liver", ANGPTL3 = "Liver")



### Loop each gene: store the gene cis-variants expression GWAS [only once]
# (convenient for future summary-data collection if the IV selection strategy changes)
for(g in genes_with_good_expression){
  tissue <- tissue_map[g]
  f <- file.path(indir, paste0(g,'_',tissue, ".tsv"))
  if(!file.exists(f)) {print('no tissue GWAS exists for this gene');next}

  cat("cis_variants GWAS:", g, "cis-variants, based on GWAS of GTEx with tissue:", tissue, "\n")

  out_f <- file.path(indir, paste0(g, "_", tissue, ".cis.tsv"))
  if(file.exists(out_f)) {print('already have cis_variants GWAS; skip'); next}

  DT <- fread(f)
  if(!nrow(DT)) next

  # ---- auto-detect rsID and p-value columns (with priority) ----
  rs_candidates <- c("rsID","variant_id","rsid","rs_id","ID","SNP")
  p_candidates  <- c("p_value","P-value","pvalue","pval_nominal")
  rs_col <- intersect(rs_candidates, names(DT))[1]
  p_col  <- intersect(p_candidates,  names(DT))[1]
  if (is.na(rs_col) || is.na(p_col)) {stop("Cannot find rsID or p-value column in ", basename(f))}

  # ---- cis rsIDs for this gene ----
  cis_rs <- cis_snps[external_gene_name == g, rsID[[1]]]
  if(length(cis_rs) == 0) next
  # ---- restrict GWAS to cis-variants only ----
  DT <- DT[get(rs_col) %in% cis_rs]
  if(!nrow(DT)) next
  # ----- store the gene cis-variants expression GWAS ----
  fwrite(DT,file.path(indir, paste0(g, "_", tissue, ".cis.tsv")),sep = "\t")
} # ~ 1 s


### Loop each gene: PLINK LD-clumping
for(g in genes_with_good_expression){
  tissue <- tissue_map[g]

  cat("clumping:", g, "cis-variants, based on GWAS of GTEx with tissue: ", tissue, "\n")

  # ----- read the gene cis-variants expression GWAS  ----
  DT<- fread(file.path(indir, paste0(g, "_", tissue, ".cis.tsv")))
  if(!nrow(DT)) next

  # ---- auto-detect rsID and p-value columns (with priority) ----
  rs_candidates <- c("rsID","variant_id","rsid","rs_id","ID","SNP")
  p_candidates  <- c("p_value","P-value","pvalue","pval_nominal")
  rs_col <- intersect(rs_candidates, names(DT))[1]
  p_col  <- intersect(p_candidates,  names(DT))[1]
  if (is.na(rs_col) || is.na(p_col)) {stop("Cannot find rsID or p-value column in ", basename(f))}

  # Re-derive p-values from |Z-score| to avoid truncation to 0 for very small p-values
  est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
  se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
  if (is.na(est_col) || is.na(se_col)) {stop("Cannot find est or se column (for |z-score| rank) in ", basename(f))}
  DT[, z := get(est_col) / get(se_col)]
  DT <- DT[order(-abs(z))]
  DT[, newP := 2 * pnorm(-abs(z))]
  p0 <- 1e-300
  DT[!is.finite(newP) | newP < p0, newP := p0]
  DT[newP == p0, newP := newP + seq_along(.I) * 1e-302]
  DT[, (p_col) := newP]

  # ---- create PLINK-format GWAS file (rsID + p) ----
  tmp <- tempfile(fileext = ".txt")
  fwrite(DT[, .(SNP = get(rs_col), P = get(p_col))], tmp, sep = "\t")
  # ---- output prefix: gene + tissue ----
  prefix <- file.path(outdir, paste0(g, "_", tissue))
  # ---- PLINK clumping (r^2 < 0.01, p < 1e-6, MAF > 0.01) ----
  cmd <- sprintf(
    '%s --bfile %s --clump %s --clump-p1 1e-6 --clump-p2 1e-4 --clump-r2 0.01 --clump-kb 500 --maf 0.01 --out %s',
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
lead_counts_expression <- lead_counts
lead_counts_expression
# ANGPTL3_Liver.clumped      1
#     LPA_Liver.clumped      1




### STEP 6: get the final exposure data and run MR analysis (independent IVs) ---------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------------------

# Current IV selection strategy: PLINK --clump-p1 1e-6 --clump-p2 1e-4 --clump-r2 0.01 --maf 0.01
lead_counts_proxy
lead_counts_protein
lead_counts_expression
gene_names

# Build the folder mapping for each gene: use the dataset yielding the largest number of independent SNPs
folder_map <- c(GLP1R ="ProxyGWAS" ,       # 3
                LPA = "ProxyGWAS",          # 29
                PCSK9 = "ProxyGWAS",        # 14
                HMGCR = "ProxyGWAS",        # 7
                APOC3 = "ProxyGWAS",        # 19
                ANGPTL3 = "ProxyGWAS"       # 5
                )

## For each gene, run MR (PCA + independent IV version) for all outcomes:
# all SNPs:          MR-PCA-GMM + F-statistic + Q p-value
# independent SNPs:  IVW + F-statistic + Q p-value + robust MR (Median, Mode, RAPS)

# Note: cisMR can still benefit from robust methods -- although cis-IVs are typically all-valid-or-all-invalid,
# some SNPs may be in LD with pleiotropic variants (analogous to SuSiE-style TWAS)


### Prepare the liver MRI GWAS outcomes, restricted to all cis-SNPs [only once; skip next time] -------------------
# (this accelerates the MR pipeline by avoiding repeated reads of large GWAS files)
retained_names <- readRDS("/path/to/project/phenotype data/features with other images/retained_names.rds")
liver_gwas_dir <- "/path/to/server/project/GWAS_regenie/GWAS_results_new"   # original GWAS data on server
over_cisSNPs_dir <- '/path/to/project/targetMR/radiomicsGWAS'

### Option A: extract cis-SNPs on the Broad Server (fast) ---------------
all_rsID <- unique(unlist(cis_snps$rsID)) # collect all cis rsIDs across genes (4625 SNPs)
rs_file_local <- "/path/to/server/project/quick_running2/cis_rsID.txt"
writeLines(all_rsID, rs_file_local) # write one rsID per line
# Then on the Broad Server:
# ssh Broad Server -> cd /path/to/server/project/quick_running2 -> R-4.3 -> bash extract_cis_snps.sh  (~30s per feature)
# Then copy all *.tsv to local:
# mv /path/to/server/... /path/to/project/targetMR/radiomicsGWAS/

### Option B: extract cis-SNPs locally (slower) ------------
all_rsID <- unique(unlist(cis_snps$rsID))    # 4625 SNPs
rs_file <- tempfile(fileext = ".txt")
writeLines(all_rsID, rs_file)
for(i in seq_along(retained_names)){
  cat(i, "- ")
  tr <- retained_names[i]
  outfile <- file.path(over_cisSNPs_dir, paste0(tr, ".tsv"))
  if (file.exists(outfile) && file.size(outfile) > 0) next  # skip if already extracted
  f <- file.path(liver_gwas_dir, paste0(tr, ".regenie"))
  if (!file.exists(f)){print('no GWAS file on server!');next}
  cmd <- sprintf( "{ head -n 1 %s; fgrep -w -f %s %s; }", shQuote(f), shQuote(rs_file), shQuote(f))
  DT <- try(fread(cmd = cmd), silent = TRUE)
  if (inherits(DT, "try-error") || !nrow(DT)) next
  fwrite(DT, outfile, sep = "\t")
}




### MR fitting --------------------------------------------------------------------------------------------------
radiomicsGWAS <- '/path/to/project/targetMR/radiomicsGWAS'
feature_paths <- list.files(radiomicsGWAS, pattern = "\\.tsv$", full.names = TRUE)
feature_paths
cisMR_res<-list()

for(gene_name in gene_names){
  cat(paste0( '\n=============================','\n' ))
  cat(paste0( 'Current gene: ', gene_name ,  '\n' ))


  folder_used <- paste0('/path/to/project/targetMR/'  , folder_map[gene_name])

  ## get the exposure's clumped cis-QTLs and the exposure GWAS data (rsID EA est se) ----------------------------
  ## ------------------------------------------------------------------------------------------------------------
  ## ----  find clumped file for this gene ----
  clump_dir <- file.path(folder_used, "LD_clumped_results")
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
    GYdata <- dat[match(GXdata$rsID ,   dat$SNP)  , ]  # should be the same nrow; liver MRI GWAS is 1000G-based so no missing rsID

    ## harmonization (according to the LD ref panel effect allele) ------------------------
    GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]
    bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
    by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
    ref <- bim[match(GXdata_$rsID, bim$SNP), ]
    flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1
    bx[flip_e] <- -bx[flip_e];by[flip_o] <- -by[flip_o]

    ## random-effect IVW (store F-statistic and Q test results) ----------------------
    MRres<-mr_ivw(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse ) )
    #mr_plot( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse )  , orientate=TRUE)

    ## robust MR methods (only work with >= 3 SNPs) --------------------------------------
    if( nrow(GYdata) >=3  ){
      MRmedian<-mr_median( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse) )
      Median_est<-MRmedian@Estimate  ; Median_p<- MRmedian@Pvalue
      MRmode<-mr_mbe(  mr_input(bx= bx, bxse= bxse ,by=by, byse=byse ) )
      Mode_est<-MRmode@Estimate  ; Mode_p<- MRmode@Pvalue
      set.seed(1123) # set a seed to ensure MR-RAPS convergence (avoids occasional subscript-out-of-bounds error for LPA)
      MRRAPS <- suppressWarnings(
        mr.raps( data.frame(beta.exposure=bx, beta.outcome=by ,se.exposure=bxse , se.outcome =byse ) ,diagnostics=FALSE, over.dispersion = TRUE )
      )
      MRRAP_est <- MRRAPS$beta.hat; MRRAP_p <-  2*pnorm(-abs(MRRAPS$beta.hat/MRRAPS$beta.se) )
    }else{
      Median_est <- Median_p <- Mode_est <- Mode_p <- MRRAP_est <- MRRAP_p <- NA
    }

    ## Store/update the vector/matrix results ---------------------------------------------
    MR_vector_results <- c( kk,sub("\\.cis\\.tsv$", "", basename(exp_file)), feature_name , nrow(GYdata) ,  MRres@Fstat , MRres@Heter.Stat[2],
                            MRres@Estimate,MRres@StdError, MRres@Pvalue, # MR-IVW is the primary method
                            Median_est  ,Median_p ,
                            Mode_est  ,Mode_p ,
                            MRRAP_est , MRRAP_p)
    MR_matrix_results <- rbind( MR_matrix_results , MR_vector_results )
  }
  # finish the loop

  colnames(MR_matrix_results) <- c( 'feature ID','Exposure', 'feature_name', 'num_of_snps','Fstatistic','Qpvalue',
                                'IVWest','IVWse','IVWp',
                                'Median_est','Median_p', 'Mode_est','Mode_p', 'MRRAP_est','MRRAP_p')
  MR_matrix_results <- as.data.frame(MR_matrix_results)
  MR_matrix_results[-(2:3)] <- lapply(MR_matrix_results[-(2:3)], as.numeric)

  MR_matrix_results$FDRp <- p.adjust(MR_matrix_results$IVWp, method = "BH")

  ### store all MR results for this gene
  cisMR_res[[gene_name]] <- MR_matrix_results

}

saveRDS(cisMR_res,
        "/path/to/project/targetMR/cisMR_res/cisMR_res.rds")
cisMR_res <- readRDS("/path/to/project/targetMR/cisMR_res/cisMR_res.rds")



View(cisMR_res$GLP1R)
View(cisMR_res$PCSK9)
View(cisMR_res$HMGCR)




### STEP 7: cisMR with correlated SNPs (MR-PC-GMM) ------------------------------------------------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------------------


# This step uses all cis-SNPs (with pairwise r^2 < 0.95 pruning) rather than only independent lead SNPs.
# Ash's PCA method (MR-PC-GMM) is used, which has the highest statistical power among LD-SNP MR methods.

# LD-SNP MR methods available: LD-SNP IVW (requires SNP selection) and MR-PC-GMM (no SNP selection needed).
# SNP selection with LD is sensitive to the reference panel -- 1000G European matches UKB-PPP and European proxy GWAS well.

# Note: GTEx eQTLs are not used here, as GTEx in-sample LD can differ substantially from the 1000G reference panel,
# making LD-aware MR estimation unreliable. For GLP1R, BMI proxy is preferred.

# Limitation of LD-SNP MR: classical robust MR methods (MR-Egger, MR-Median) are not applicable.
# Recall LD-aware robust alternatives: MR-Egger, MRmedian, cisMR-cMR.


### About MR-PCA ---------------------------------------------------------------
# Input SNPs: no pairwise r^2 > 0.95, and nominally associated with the exposure (p < 0.05)
# Inputs: GX_est, GX_se, GY_est, GY_se, LD matrix (all harmonized to LD ref effect allele), nx, ny
# Outputs: MR-PCA-GMM estimate, SE, p-value; F-statistic; Q p-value

## Technical note on MR-PC-GMM:
# MR-PC-GMM uses GMM (generalized method of moments) for inference, based on multivariate regression.
# This is why nx and ny (sample sizes) are needed: to transform univariate summary statistics
# to multivariate regression summary statistics.
# robust = TRUE: introduces overdispersion, allowing residual variance > theoretical value
# (analogous to multiplicative random effects).

# Sensitivity analysis on ny input: the MR-PC-GMM results are not sensitive to the ny input,
# as ny is largely offset in the key formula when genetic effect sizes are small.


### Prepare the sample size for each exposure GWAS
# Use UKB-PPP or proxy GWAS only (not GTEx, due to LD mismatch concerns)
folder_map <- c(GLP1R ="ProxyGWAS" ,       # 3 (MVP GWAS)
                LPA = "ProxyGWAS",          # 22
                PCSK9 = "ProxyGWAS",        # 13
                HMGCR = "ProxyGWAS",        # 5
                APOC3 = "ProxyGWAS",        # 17
                ANGPTL3 = "ProxyGWAS"       # 6
)

## Method 1: get sample size from the N column in the GWAS data
exposure_samplesize_map<-c()
for(gene_name in gene_names){
  cat(paste0( '\n=============================','\n' ))
  cat(paste0( 'Current gene: ', gene_name ,  '\n' ))
  folder_used <- paste0('/path/to/project/targetMR/'  , folder_map[gene_name])
  exp_file <- list.files(folder_used, pattern = gene_name, full.names = TRUE)
  exp_file <- exp_file[grepl("cis\\.tsv$", exp_file)][1]  # if multiple matches (e.g., BMI), the first is BMI_MVP
  DT <- data.table::fread(exp_file)
  # get the median sample size
  n_col  <- intersect(c("n","N","sample_size"), names(DT))[1]
  if(is.na(n_col)){      exposure_samplesize_map <- c( exposure_samplesize_map, NA )     }else{
    GXdata <- DT[, .(N   = get(n_col) )]
    exposure_samplesize_map <- c( exposure_samplesize_map, median( GXdata$N ) )
  }
}
names(exposure_samplesize_map ) <- gene_names

## Method 2: approximate N from se and MAF using nx ~= 1 / (GXse^2 * 2 * MAF * (1 - MAF))
# This approximation requires Var(phenotype) = 1 (e.g., RINT-normalized traits) and small genetic effects.
exposure_samplesize_map2<-c()
for(gene_name in gene_names){
  cat(paste0( '\n=============================','\n' ))
  cat(paste0( 'Current gene: ', gene_name ,  '\n' ))
  folder_used <- paste0('/path/to/project/targetMR/'  , folder_map[gene_name])
  exp_file <- list.files(folder_used, pattern = gene_name, full.names = TRUE)
  exp_file <- exp_file[grepl("cis\\.tsv$", exp_file)][1]
  DT <- data.table::fread(exp_file)
  rs_col  <- intersect(c("rsID","rs_id","rsid","variant_id","ID","SNP"), names(DT))[1]
  se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
  eaf_col  <- intersect(c("effect_allele_frequency","A1FREQ", "MinFreq"), names(DT))[1]
  GXdata <- DT[ , .(rsID = get(rs_col), se   = get(se_col) , eaf  = get(eaf_col))]
  putative_n <- median( 1/(  GXdata$se^2 * 2 * GXdata$eaf * (1- GXdata$eaf) )  )
  exposure_samplesize_map2 <- c( exposure_samplesize_map2, median( putative_n ) )
}
names(exposure_samplesize_map2 ) <- gene_names

exposure_samplesize_map;exposure_samplesize_map2
# The two methods give similar results.
# Strategy: use N column directly if available; otherwise use the approximation.
exposure_samplesize_map


### Some useful functions
prune_ld95 <- function(LD, thr = 0.95){  # returns rsID vector of SNPs with no pairwise r^2 > thr
  keep <- colnames(LD)
  R2   <- LD^2
  diag(R2) <- 0
  repeat{
    if (max(R2) <= thr) break
    # drop the SNP with the most r^2 > thr connections
    drop <- names(which.max(rowSums(R2 > thr)))
    keep <- setdiff(keep, drop)
    R2   <- R2[keep, keep, drop = FALSE]
  }
  return(keep)
}


### MR-PCA fitting --------------------------------------------------------------------------------------------------
radiomicsGWAS <- '/path/to/project/targetMR/radiomicsGWAS'
feature_paths <- list.files(radiomicsGWAS, pattern = "\\.tsv$", full.names = TRUE)

cisMR_PCA_res<-list()

for(gene_name in gene_names){
  cat(paste0( '\n=============================','\n' ))
  cat(paste0( 'Current gene: ', gene_name ,  '\n' ))


  folder_used <- paste0('/path/to/project/targetMR/'  , folder_map[gene_name])

  ## QC cis-SNPs: prune to pairwise r^2 < 0.95 using the 1000G LD reference ----------------------------
  cat('getting LD matrix for pruning pairwise r^2 < 0.95 ... \n'  )
  rs <- unlist(cis_snps[external_gene_name == gene_name, rsID][[1]])
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
  # reorder LD matrix to match rs order
  idx <- match(rs, ord)
  if (any(is.na(idx))) stop("Some rsID not found in LD reference panel")
  LD <- LD0[idx, idx]
  colnames(LD)<-rownames(LD)<- rs
  # prune to no pairwise r^2 > 0.95
  rs_ld95 <- prune_ld95(LD, thr = 0.95)
  cat(paste0(   'finished pruning: from ' , length(rs) , ' to ', length(rs_ld95) , ' SNPs (r^2 < 0.95) \n' )     )

  ## get the exposure GWAS data for QC-ed SNPs (rsID EA est se) ----------------------------
  exp_file <- list.files(folder_used, pattern = gene_name, full.names = TRUE)
  exp_file <- exp_file[grepl("cis\\.tsv$", exp_file)][1]  # if multiple matches, the first is BMI_MVP
  DT <- data.table::fread(exp_file)
  # Note: duplicated rsIDs can occur for non-biallelic SNPs; resolved by keeping the one with higher MAF
  ## ---- auto-detect columns (priority order) ----
  rs_col  <- intersect(c("rsID","rs_id","rsid","variant_id","ID","SNP"), names(DT))[1]
  ea_col  <- intersect(c("EA","alt","ALLELE1","effect_allele","Allele1"), names(DT))[1]
  est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
  se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
  eaf_col  <- intersect(c("effect_allele_frequency","A1FREQ", "MinFreq"), names(DT))[1]
  ## ---- restrict to QC-ed SNPs (r^2 < 0.95) & standardize output ----
  GXdata <- DT[get(rs_col) %in% rs_ld95, .(rsID = get(rs_col),  EA   = get(ea_col) , est  = get(est_col), se   = get(se_col), eaf =get(eaf_col)  )]
  ## for duplicated rsIDs, retain the one with higher MAF
  GXdata <- GXdata[order(rsID, -pmin(eaf, 1 - eaf))]
  GXdata <- GXdata[!duplicated(rsID)]
  ## ensure MAF > 0.01 (avoids bad GMM approximation)
  GXdata <- GXdata[pmin(GXdata$eaf, 1 - GXdata$eaf) > 0.01, ]
  ## only retain nominally significant exposure SNPs (p < 0.05)
  GXdata <- GXdata[   2*pnorm( -abs(GXdata$est/GXdata$se) ) <0.05  ,   ]
  cat(  paste0( 'final num of SNPs for MR-PCA (p < 0.05; MAF > 0.01; LD r^2 < 0.95): ', nrow(GXdata) ,'\n' ) )



  ## loop all outcomes (the 59 retained liver MRI features) -----------------------------------------------------
  ## ------------------------------------------------------------------------------------------------------------
  nx_used<-exposure_samplesize_map[gene_name] ; ny_used <- 37725
  cat(  paste0(   'nx and ny used in MR-PCA-GMM: '  ,  nx_used, ' ',ny_used, '\n'  ) )
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
    GYdata <- na.omit(GYdata)   # remove NA rows
    # liver MRI GWAS is 1000G-based so missing rsIDs are rare (only occasionally removed by REGENIE QC)

    ## harmonization (according to the LD ref panel effect allele) ------------------------
    GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]  # ensure same order for GX and GY
    bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
    by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
    ref <- bim[match(GXdata_$rsID, bim$SNP), ]
    flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1
    bx[flip_e] <- -bx[flip_e];  by[flip_o] <- -by[flip_o]

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


    ## MR-PCA-GMM fitting (store F-statistic and Q test results) ----------------------
    MRPCAres <- mr_pcgmm(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse, correlation = LD ) ,  robust = FALSE, thres=0.99,
                      nx =nx_used , ny= ny_used )
    MRPCAres_robust <- mr_pcgmm(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse, correlation = LD ) , robust= TRUE,thres=0.99,
                         nx =nx_used , ny= ny_used )

    ## Store/update the vector/matrix results ---------------------------------------------
    MR_vector_results <- c( kk, sub("\\.cis\\.tsv$", "", basename(exp_file)) ,feature_name ,
                            nrow(GYdata) , MRPCAres@PCs,
                            MRPCAres@Fstat , MRPCAres@Heter.Stat[1] , MRPCAres@Heter.Stat[2],MRPCAres_robust@Heter.Stat,
                            MRPCAres@Estimate, MRPCAres@StdError , MRPCAres@CILower, MRPCAres@CIUpper,  MRPCAres@Pvalue,
                            MRPCAres_robust@Estimate, MRPCAres_robust@StdError ,MRPCAres_robust@CILower , MRPCAres_robust@CIUpper  ,MRPCAres_robust@Pvalue )
    MR_matrix_results <- rbind( MR_matrix_results , MR_vector_results )
  }
  ### finish the loop!

  colnames(MR_matrix_results) <- c( 'feature ID', 'Exposure','feature name',
                                    'num_of_snps', 'num_of_PCs',
                                    'Fstatistic', 'Qvalue','Qpvalue', 'Qvalue_robust',
                                    'est','se','CIlow','CIup','pvalue',
                                    'est_robust','se_robust','CIlow_robust','CIup_robust','pvalue_robust')

  MR_matrix_results <- as.data.frame(MR_matrix_results)
  MR_matrix_results[-(2:3)] <- lapply(MR_matrix_results[-(2:3)], as.numeric)

  MR_matrix_results$FDRp <- p.adjust(MR_matrix_results$pvalue_robust, method = "BH")

  ### store all MR results for this gene
  cisMR_PCA_res[[gene_name]] <- MR_matrix_results

}


saveRDS(cisMR_PCA_res,
        "/path/to/project/targetMR/cisMR_res/cisMR_PCA_res.rds")
cisMR_PCA_res <- readRDS("/path/to/project/targetMR/cisMR_res/cisMR_PCA_res.rds")


View(cisMR_PCA_res$GLP1R)  # note: GLP1R uses MVP BMI to avoid weak IV (F < 10 with UKB+GIANT BMI)
View(cisMR_PCA_res$LPA)
View(cisMR_PCA_res$PCSK9)
View(cisMR_PCA_res$HMGCR)
View(cisMR_PCA_res$APOC3)   # 1 MRI outcome passing FDR p-value < 0.05
View(cisMR_PCA_res$ANGPTL3)

# Conclusion: only one MRI feature passes FDR correction: APOC3 -> TG -> glszm_ZoneEntropy_inp


### MR-PC-GMM vs. MR-IVW scatter plot (sanity check)

plot(cisMR_PCA_res$HMGCR$est_robust, cisMR_res$HMGCR$IVWest,
     xlab='MR-PC-GMM estimates',ylab='MR-IVW estimates')
abline(0,1,col='red')

genes <- intersect(names(cisMR_PCA_res), names(cisMR_res))
par(mfrow=c(2,3))
for (g in genes) {
  x  <- cisMR_PCA_res[[g]]$est_robust
  y  <- cisMR_res[[g]]$IVWest
  xs <- cisMR_PCA_res[[g]]$se_robust
  ys <- cisMR_res[[g]]$IVWse

  plot(x, y, pch = 16, cex = 1.2, col = "black",main=g,
       xlab = "MR-PC-GMM estimates", ylab = "MR-IVW estimates")
  abline(0,1,col='blue')
}
par(mfrow=c(1,1))




### Step 8: Positive control using MR-PC-GMM: each gene -> exposure -> CAD or T2D ----------------------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------------------------------

### Original CAD GWAS: meta-analysis GWAS, EUR participants
### Original T2D GWAS: meta-analysis GWAS, EUR participants


### STEP 1: prepare the cis CAD or T2D GWAS data [only once; skip next time] ----------------------------------------------
all_rsID <- unique(unlist(cis_snps$rsID))    # 4625 SNPs
originalfile <- '/path/to/project/targetMR/otherGWAS/CAD.tsv'
DT <- fread(originalfile)
# ---- auto-detect rsID columns (with priority) ----
rs_candidates <- c("rsID","rs_id","rsid","variant_id","ID","SNP")
rs_col <- intersect(rs_candidates, names(DT))[1]
# ---- restrict GWAS to cis-variants only ----
DT <- DT[get(rs_col) %in% all_rsID]
if(!nrow(DT)) { print('no rsID column detected!') }
# ----- store the cis-restricted GWAS ----
outfile <- '/path/to/project/targetMR/otherGWAS/CAD.cis.tsv'
fwrite(DT,outfile,sep = "\t")

### STEP 2: note on ny input for MR-PC-GMM --------------------------------
# MR-PC-GMM uses ny in a rigorous formula; however, the results are not sensitive to the ny input
# because ny is effectively offset when genetic effect sizes are small.
# Sensitivity checks confirm stability across different ny values.


### STEP 3: loop all SNPs, and obtain the MR-PC-GMM results on CAD/T2D ------------------------

positive_outcome_map<-c(GLP1R ="T2D" ,
                        LPA = "CAD",PCSK9 = "CAD",HMGCR = "CAD",APOC3 = "CAD",ANGPTL3 = "CAD" )


Postivie_cisMR_res<-c()
for(gene_name in gene_names){
  cat(paste0( '\n=============================','\n' ))
  cat(paste0( 'Current gene: ', gene_name ,  '\n' ))


  folder_used <- paste0('/path/to/project/targetMR/'  , folder_map[gene_name])

  ## QC cis-SNPs: prune to pairwise r^2 < 0.95 ----------------------------
  cat('getting LD matrix for pruning pairwise r^2 < 0.95 ... \n'  )
  rs <- unlist(cis_snps[external_gene_name == gene_name, rsID][[1]])
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
  # reorder LD matrix to match rs order
  idx <- match(rs, ord)
  if (any(is.na(idx))) stop("Some rsID not found in LD reference panel")
  LD <- LD0[idx, idx]
  colnames(LD)<-rownames(LD)<- rs
  # prune to no pairwise r^2 > 0.95
  rs_ld95 <- prune_ld95(LD, thr = 0.95)
  cat(paste0(   'finished pruning: from ' , length(rs) , ' to ', length(rs_ld95) , ' SNPs \n' )    )

  ## get the exposure GWAS data for QC-ed SNPs ----------------------------
  exp_file <- list.files(folder_used, pattern = gene_name, full.names = TRUE)
  exp_file <- exp_file[grepl("cis\\.tsv$", exp_file)][1]  # if multiple matches, the first is BMI_MVP
  DT <- data.table::fread(exp_file)
  # Note: duplicated rsIDs can occur for non-biallelic SNPs; resolved by keeping the one with higher MAF
  ## ---- auto-detect columns (priority order) ----
  rs_col  <- intersect(c("rsID","rs_id","rsid","variant_id","ID","SNP"), names(DT))[1]
  ea_col  <- intersect(c("EA","alt","ALLELE1","effect_allele","Allele1"), names(DT))[1]
  est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
  se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
  eaf_col  <- intersect(c("effect_allele_frequency","A1FREQ", "MinFreq","af"), names(DT))[1]
  ## ---- restrict to QC-ed SNPs & standardize output ----
  GXdata <- DT[get(rs_col) %in% rs_ld95, .(rsID = get(rs_col),  EA   = get(ea_col) , est  = get(est_col), se   = get(se_col), eaf =get(eaf_col)  )]
  ## for duplicated rsIDs, retain the one with higher MAF
  GXdata <- GXdata[order(rsID, -pmin(eaf, 1 - eaf))]
  GXdata <- GXdata[!duplicated(rsID)]
  ## ensure MAF > 0.01
  GXdata <- GXdata[pmin(GXdata$eaf, 1 - GXdata$eaf) > 0.01, ]
  ## only retain nominally significant exposure SNPs
  GXdata <- GXdata[   2*pnorm( -abs(GXdata$est/GXdata$se) ) <0.05  ,   ]
  cat(  paste0( 'final num of SNPs for MR-PCA (p < 0.05): ', nrow(GXdata) ,'\n' ) )


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
  GYdata <- dat[match(GXdata$rsID ,   dat$SNP)  , ]
  GYdata <- na.omit(GYdata)   # remove NA rows

  ## harmonization (according to the LD ref panel effect allele) ------------------------
  GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]  # ensure same order for GX and GY
  bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
  by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
  ref <- bim[match(GXdata_$rsID, bim$SNP), ]
  flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1; flip alleles if = A2
  bx[flip_e] <- -bx[flip_e];  by[flip_o] <- -by[flip_o]

  ## get the corresponding LD correlation matrix
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


  ## get ny for MR-PC-GMM (not sensitive to the choice of ny)
  ny_used <- median( GYdata$N  )  # total sample size


  ## MR-PCA-GMM fitting (store F-statistic and Q test results) ----------------------
  nx_used<-exposure_samplesize_map[gene_name]
  MRPCAres <- mr_pcgmm(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse, correlation = LD ) ,  robust = FALSE, thres=0.99,
                       nx =nx_used , ny= ny_used )
  MRPCAres_robust <- tryCatch(
    mr_pcgmm(
      mr_input(bx = bx, bxse = bxse, by = by, byse = byse, correlation = LD),
      robust = TRUE, thres = 0.99,
      nx = nx_used, ny = ny_used
    ),
    error = function(e){
      warning(sprintf(
        "mr_pcgmm robust failed for gene %s: %s; fallback to non-robust result",
        gene_name, e$message
      ))
      MRPCAres
    }
  )


  ## Store the results
  res <-  c( gene_name ,sub("\\.cis\\.tsv$", "", basename(exp_file)) , positive_outcome_map[gene_name] ,
             length(bx), MRPCAres@PCs,  MRPCAres@Fstat,
             MRPCAres@Estimate, MRPCAres@CILower, MRPCAres@CIUpper,  MRPCAres@Pvalue ,  MRPCAres@Heter.Stat[1], MRPCAres@Heter.Stat[2],
             MRPCAres_robust@Estimate ,MRPCAres_robust@CILower, MRPCAres_robust@CIUpper, MRPCAres_robust@Pvalue, MRPCAres_robust@Heter.Stat[1])
  Postivie_cisMR_res<-rbind( Postivie_cisMR_res , res )
}

colnames(Postivie_cisMR_res) <- c( 'gene','Exposure','Outcome', 'num_of_SNPs','num_of_PCs', 'Fstatistic',
                              'Estimate', 'CI_low','CI_up','pvalue','Qvalue','Qpvalue',
                              'Estimate_robust','CI_low_robust','CI_upper_robust','pvalue_robust','Qvalue_robust')
Postivie_cisMR_res <- as.data.frame(Postivie_cisMR_res)
Postivie_cisMR_res[-(1:3)] <- lapply(Postivie_cisMR_res[-(1:3)], as.numeric)

saveRDS(Postivie_cisMR_res,
        "/path/to/project/targetMR/cisMR_res/Positive_cisMR_PCA_res.rds")
Postivie_cisMR_res <- readRDS("/path/to/project/targetMR/cisMR_res/Positive_cisMR_PCA_res.rds")


Postivie_cisMR_res
### MR-PC-GMM
#    gene      Exposure Outcome num_of_SNPs num_of_PCs Fstatistic   Estimate      CI_low     CI_up        pvalue
#   GLP1R GLP1R_BMI_MVP     T2D         229         13   25.79972 0.65841456  0.41384795 0.9029812  1.316360e-07
#     LPA      LPA_Lp_a     CAD         274         14 7137.21842 0.25398818  0.23787097 0.2701054 1.801024e-209
#   PCSK9   PCSK9_LDL_C     CAD         253         16  510.72125 0.59143599  0.50722286 0.6756491  4.139175e-43
#   HMGCR   HMGCR_LDL_C     CAD         157          4  608.02874 0.33995333  0.19698162 0.4829250  3.156930e-06
#   APOC3      APOC3_TG     CAD         265          8 1918.77352 0.24667789  0.19064414 0.3027116  6.223369e-18
# ANGPTL3    ANGPTL3_TG     CAD          71          3  756.06032 0.05155968 -0.09902637 0.2021457  5.021703e-01




### Positive control: MR-IVW (independent SNPs)
MR_matrix_results<-c()
for(gene_name in gene_names){
  cat(paste0( '\n=============================','\n' ))
  cat(paste0( 'Current gene: ', gene_name ,  '\n' ))


  folder_used <- paste0('/path/to/project/targetMR/'  , folder_map[gene_name])

  ## get the exposure's clumped cis-QTLs and the exposure GWAS data (rsID EA est se) ----------------------------
  ## ------------------------------------------------------------------------------------------------------------
  ## ----  find clumped file for this gene ----
  clump_dir <- file.path(folder_used, "LD_clumped_results")
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
  rs_col  <- intersect(c("rsID","rs_id","variant_id","rsid","ID","SNP"), names(DT))[1]
  ea_col  <- intersect(c("EA","alt","ALLELE1","effect_allele","Allele1"), names(DT))[1]
  est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
  se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
  ## ---- restrict to clumped SNPs & standardize output ----
  GXdata <- DT[get(rs_col) %in% clumped_rs, .(rsID = get(rs_col),  EA   = get(ea_col) , est  = get(est_col), se   = get(se_col) )]


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
    GYdata <- dat[match(GXdata$rsID ,   dat$SNP)  , ]
    GYdata <- na.omit(GYdata)   # remove NA rows

    ## harmonization (according to the LD ref panel effect allele) ------------------------
    GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]
    bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
    by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
    ref <- bim[match(GXdata_$rsID, bim$SNP), ]
    flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1
    bx[flip_e] <- -bx[flip_e];by[flip_o] <- -by[flip_o]

    ## random-effect IVW (store F-statistic and Q test results) ----------------------
    MRres<-mr_ivw(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse ) )
    mr_plot( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse )  , orientate=TRUE)

    ## robust MR methods (only work with >= 3 SNPs) --------------------------------------
    if( nrow(GYdata) >=3  ){
      MRmedian<-mr_median( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse) )
      Median_est<-MRmedian@Estimate  ; Median_p<- MRmedian@Pvalue
      MRmode<-mr_mbe(  mr_input(bx= bx, bxse= bxse ,by=by, byse=byse ) )
      Mode_est<-MRmode@Estimate  ; Mode_p<- MRmode@Pvalue
      set.seed(1123)
      MRRAPS <- suppressWarnings(
        mr.raps( data.frame(beta.exposure=bx, beta.outcome=by ,se.exposure=bxse , se.outcome =byse ) ,diagnostics=FALSE, over.dispersion = TRUE )
      )
      MRRAP_est <- MRRAPS$beta.hat; MRRAP_p <-  2*pnorm(-abs(MRRAPS$beta.hat/MRRAPS$beta.se) )
    }else{
      Median_est <- Median_p <- Mode_est <- Mode_p <- MRRAP_est <- MRRAP_p <- NA
    }

    ## Store/update the vector/matrix results ---------------------------------------------
    MR_vector_results <- c( gene_name , folder_map[gene_name], positive_outcome_map[gene_name] ,
                            nrow(GYdata) ,  MRres@Fstat , MRres@Heter.Stat[2],
                            MRres@Estimate,MRres@StdError ,MRres@Pvalue,
                            Median_est  ,Median_p ,
                            Mode_est  ,Mode_p ,
                            MRRAP_est , MRRAP_p)
    MR_matrix_results <- rbind( MR_matrix_results , MR_vector_results )
}
# finish the loop
colnames(MR_matrix_results) <- c( 'gene',  'Exposure', 'Outcome',
                                  'num_of_snps','Fstatistic','Qpvalue',
                                    'IVWest', 'IVWse', 'IVWp',
                                    'Median_est','Median_p', 'Mode_est','Mode_p', 'MRRAP_est','MRRAP_p')

MR_matrix_results <- as.data.frame(MR_matrix_results)
MR_matrix_results[-(1:3)] <- lapply(MR_matrix_results[-(1:3)], as.numeric)
MR_matrix_results
saveRDS(MR_matrix_results,
        "/path/to/project/targetMR/cisMR_res/Positive_MRIVW_res.rds")
MR_matrix_results <- readRDS("/path/to/project/targetMR/cisMR_res/Positive_MRIVW_res.rds")

### MR-IVW
#    gene  Exposure Outcome num_of_snps Fstatistic      Qpvalue     IVWest      IVWse         IVWp
#   GLP1R ProxyGWAS     T2D           2   78.91369 5.212829e-01 0.34325735 0.17223976 4.627193e-02
#     LPA ProxyGWAS     CAD          29 3708.48224 2.017663e-74 0.27170069 0.02992091 1.079501e-19
#   PCSK9 ProxyGWAS     CAD          14  692.43476 5.863010e-01 0.60966986 0.03986237 8.336361e-53
#   HMGCR ProxyGWAS     CAD           7  387.86833 9.915915e-02 0.35133152 0.09069150 1.070990e-04
#   APOC3 ProxyGWAS     CAD          19 1068.16057 1.602229e-01 0.23093946 0.02884898 1.193536e-15
# ANGPTL3 ProxyGWAS     CAD           4  618.73220 9.892162e-01 0.05062726 0.07103690 4.760379e-01

## PCSK9: 0.610 (0.532, 0.688)
## HMGCR: 0.351 (0.173, 0.530)
# Note: HMGCR effect on CAD is underestimated relative to PCSK9. Both target LDL-C reduction;
# however, HMGCR has fewer independent cis-IVs, which reduces the precision of the estimate.




### Negative control ------------------

### STEP 1: prepare the cis negative control GWAS data [only once; skip next time] ----------------------------------------------
all_rsID <- unique(unlist(cis_snps$rsID))    # 4625 SNPs
originalfile <- '/path/to/project/targetMR/otherGWAS/Hair_color.tsv'
DT <- fread(originalfile)
# ---- auto-detect rsID columns (with priority) ----
rs_candidates <- c("rsID","rs_id","rsid","variant_id","ID","SNP")
rs_col <- intersect(rs_candidates, names(DT))[1]
# ---- restrict GWAS to cis-variants only ----
DT <- DT[get(rs_col) %in% all_rsID]
if(!nrow(DT)) { print('no rsID column detected!') }
# ----- store the cis-restricted GWAS ----
outfile <- '/path/to/project/targetMR/otherGWAS/Hair_color.cis.tsv'
fwrite(DT,outfile,sep = "\t")

### STEP 2: MR-PC-GMM is not sensitive to ny input

### STEP 3: loop all SNPs, and obtain the MR-PC-GMM results for tanning / hair color ------------------------
negative_outcome_map<-c(GLP1R ="Hair_color" ,
                        LPA = "Hair_color",PCSK9 = "Hair_color",HMGCR = "Hair_color",APOC3 = "Hair_color",ANGPTL3 = "Hair_color" )

negative_outcome_map<-c(GLP1R ="Tanning" ,
                        LPA = "Tanning",PCSK9 = "Tanning",HMGCR = "Tanning",APOC3 = "Tanning",ANGPTL3 = "Tanning" )

Negative_cisMR_res<-c()
for(gene_name in gene_names){
  cat(paste0( '\n=============================','\n' ))
  cat(paste0( 'Current gene: ', gene_name ,  '\n' ))


  folder_used <- paste0('/path/to/project/targetMR/'  , folder_map[gene_name])

  ## QC cis-SNPs: prune to pairwise r^2 < 0.95 ----------------------------
  cat('getting LD matrix for pruning pairwise r^2 < 0.95 ... \n'  )
  rs <- unlist(cis_snps[external_gene_name == gene_name, rsID][[1]])
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
  # reorder LD matrix to match rs order
  idx <- match(rs, ord)
  if (any(is.na(idx))) stop("Some rsID not found in LD reference panel")
  LD <- LD0[idx, idx]
  colnames(LD)<-rownames(LD)<- rs
  # prune to no pairwise r^2 > 0.95
  rs_ld95 <- prune_ld95(LD, thr = 0.95)
  cat(paste0(   'finished pruning: from ' , length(rs) , ' to ', length(rs_ld95) , ' SNPs \n' )     )

  ## get the exposure GWAS data for QC-ed SNPs ----------------------------
  exp_file <- list.files(folder_used, pattern = gene_name, full.names = TRUE)
  exp_file <- exp_file[grepl("cis\\.tsv$", exp_file)][1]  # if multiple matches, the first is BMI_MVP
  DT <- data.table::fread(exp_file)
  # Note: duplicated rsIDs can occur for non-biallelic SNPs; resolved by keeping the one with higher MAF
  ## ---- auto-detect columns (priority order) ----
  rs_col  <- intersect(c("rsID","rs_id","rsid","variant_id","ID","SNP"), names(DT))[1]
  ea_col  <- intersect(c("EA","alt","ALLELE1","effect_allele","Allele1"), names(DT))[1]
  est_col <- intersect(c("beta","Effect","slope","BETA"), names(DT))[1]
  se_col  <- intersect(c("standard_error","StdErr","slope_se","SE"), names(DT))[1]
  eaf_col  <- intersect(c("effect_allele_frequency","A1FREQ", "MinFreq","af"), names(DT))[1]
  ## ---- restrict to QC-ed SNPs & standardize output ----
  GXdata <- DT[get(rs_col) %in% rs_ld95, .(rsID = get(rs_col),  EA   = get(ea_col) , est  = get(est_col), se   = get(se_col), eaf =get(eaf_col)  )]
  ## for duplicated rsIDs, retain the one with higher MAF
  GXdata <- GXdata[order(rsID, -pmin(eaf, 1 - eaf))]
  GXdata <- GXdata[!duplicated(rsID)]
  ## ensure MAF > 0.01
  GXdata <- GXdata[pmin(GXdata$eaf, 1 - GXdata$eaf) > 0.01, ]
  ## only retain nominally significant exposure SNPs
  GXdata <- GXdata[   2*pnorm( -abs(GXdata$est/GXdata$se) ) <0.05  ,   ]
  cat(  paste0( 'final num of SNPs for MR-PCA (p < 0.05): ', nrow(GXdata) ,'\n' ) )


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
  GXdata_ <- GXdata[ match(GYdata$SNP ,GXdata$rsID ),   ]  # ensure same order for GX and GY
  bx   <- GXdata_$est;bxse <- GXdata_$se;EA_e <- toupper(GXdata_$EA) # exposure
  by   <- GYdata$BETA;byse <- GYdata$SE;EA_o <- toupper(GYdata$ALLELE1)  # outcome
  ref <- bim[match(GXdata_$rsID, bim$SNP), ]
  flip_e <- EA_e == ref$A2  ; flip_o <- EA_o == ref$A2 # harmonize both to bim$A1; flip alleles if = A2
  bx[flip_e] <- -bx[flip_e];  by[flip_o] <- -by[flip_o]

  ## get the corresponding LD correlation matrix
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

  ## get ny for MR-PC-GMM (not sensitive to the choice of ny)
  ny_used <- median( GYdata$N  )  # total sample size



  ## MR-PCA-GMM fitting (store F-statistic and Q test results) ----------------------
  nx_used<-exposure_samplesize_map[gene_name]
  MRPCAres <- mr_pcgmm(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse, correlation = LD ) ,  robust = FALSE, thres=0.99,
                       nx =nx_used , ny= ny_used )
  MRPCAres_robust <- mr_pcgmm(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse, correlation = LD ) , robust= TRUE,thres=0.99,
                              nx =nx_used , ny= ny_used )


  ## Store the results
  res <-  c( gene_name ,sub("\\.cis\\.tsv$", "", basename(exp_file)) , negative_outcome_map[gene_name] ,
             length(bx), MRPCAres@PCs,  MRPCAres@Fstat,
             MRPCAres@Estimate, MRPCAres@CILower, MRPCAres@CIUpper,  MRPCAres@Pvalue ,  MRPCAres@Heter.Stat[1], MRPCAres@Heter.Stat[2],
             MRPCAres_robust@Estimate ,MRPCAres_robust@CILower, MRPCAres_robust@CIUpper, MRPCAres_robust@Pvalue, MRPCAres_robust@Heter.Stat)
  Negative_cisMR_res<-rbind( Negative_cisMR_res , res )
}
colnames(Negative_cisMR_res) <- c( 'gene','Exposure','Outcome', 'num_of_SNPs','num_of_PCs', 'Fstatistic',
                                   'Estimate', 'CI_low','CI_up','pvalue','Qvalue','Qpvalue',
                                   'Estimate_robust','CI_low_robust','CI_upper_robust','pvalue_robust','Qvalue_robust')
Negative_cisMR_res <- as.data.frame(Negative_cisMR_res)
Negative_cisMR_res[-(1:3)] <- lapply(Negative_cisMR_res[-(1:3)], as.numeric)

saveRDS(Negative_cisMR_res,
        "/path/to/project/targetMR/cisMR_res/Negative_cisMR_PCA_res.rds")
Negative_cisMR_res <- readRDS("/path/to/project/targetMR/cisMR_res/Negative_cisMR_PCA_res.rds")

Negative_cisMR_res
### MR-PC-GMM
#   gene      Exposure Outcome num_of_SNPs num_of_PCs Fstatistic      Estimate      CI_low        CI_up     pvalue
#  GLP1R GLP1R_BMI_MVP Tanning         236         13   26.31195 -0.0095027863 -0.09419315  0.075187576 0.82593333
#    LPA      LPA_Lp_a Tanning         274         14 7276.82403 -0.0071915612 -0.01411618 -0.000266946 0.04179846
#  PCSK9   PCSK9_LDL_C Tanning         246         15  561.45963 -0.0120507656 -0.04690998  0.022808450 0.49805333
#  HMGCR   HMGCR_LDL_C Tanning         156          4  603.58592  0.0006691044 -0.06482188  0.066160091 0.98402387
#  APOC3      APOC3_TG Tanning         265          8 1915.07275  0.0015219144 -0.02592136  0.028965193 0.91344585
#ANGPTL3    ANGPTL3_TG Tanning          70          3  751.98123  0.0861366535  0.01912247  0.153150837 0.01176088

# Both tanning ability and hair color negative controls show no consistent signal, supporting specificity.


### Negative control for MR-IVW (independent SNPs)
MR_matrix_results<-c()
for(gene_name in gene_names){
  cat(paste0( '\n=============================','\n' ))
  cat(paste0( 'Current gene: ', gene_name ,  '\n' ))


  folder_used <- paste0('/path/to/project/targetMR/'  , folder_map[gene_name])

  ## get the exposure's clumped cis-QTLs and the exposure GWAS data (rsID EA est se) ----------------------------
  ## ------------------------------------------------------------------------------------------------------------
  ## ----  find clumped file for this gene ----
  clump_dir <- file.path(folder_used, "LD_clumped_results")
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
  rs_col  <- intersect(c("rsID","rs_id","variant_id","rsid","ID","SNP"), names(DT))[1]
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

  ## random-effect IVW (store F-statistic and Q test results) ----------------------
  MRres<-mr_ivw(mr_input(bx= bx, bxse= bxse ,by=by, byse=byse ) )
  mr_plot( mr_input(bx= bx, bxse= bxse ,by=by, byse=byse )  , orientate=TRUE)

  ## robust MR methods not used for positive/negative controls

  ## Store/update the vector/matrix results ---------------------------------------------
  MR_vector_results <- c( gene_name , folder_map[gene_name], negative_outcome_map[gene_name] ,
                          nrow(GYdata) ,  MRres@Fstat , MRres@Heter.Stat[2],
                          MRres@Estimate,MRres@StdError , MRres@CILower, MRres@CIUpper , MRres@Pvalue )
  MR_matrix_results <- rbind( MR_matrix_results , MR_vector_results )
}
# finish the loop
colnames(MR_matrix_results) <- c( 'gene',  'Exposure', 'Outcome',
                                  'num_of_snps','Fstatistic','Qpvalue',
                                  'IVWest', 'IVWse',  'IVW_CI_low', 'IVW_CI_up' , 'IVWp')

MR_matrix_results <- as.data.frame(MR_matrix_results)
MR_matrix_results[-(1:3)] <- lapply(MR_matrix_results[-(1:3)], as.numeric)
MR_matrix_results
### MR-IVW
#    gene  Exposure Outcome num_of_snps Fstatistic    Qpvalue       IVWest       IVWse   IVW_CI_low   IVW_CI_up      IVWp
#   GLP1R ProxyGWAS Tanning           3   109.7102 0.79268286 -0.016530643 0.043457500 -0.101705777 0.068644491 0.7036586
#     LPA ProxyGWAS Tanning          28  3768.0927 0.85532251 -0.002675395 0.003232692 -0.009011356 0.003660565 0.4078937
#   PCSK9 ProxyGWAS Tanning          13   737.8408 0.26489669 -0.014071566 0.018509277 -0.050349083 0.022205951 0.4471088
#   HMGCR ProxyGWAS Tanning           6   447.5651 0.67427524 -0.015120242 0.031864300 -0.077573122 0.047332638 0.6351293
#   APOC3 ProxyGWAS Tanning          19  1068.1606 0.06778408  0.006430661 0.014867016 -0.022708154 0.035569477 0.6653450
# ANGPTL3 ProxyGWAS Tanning           4   618.7322 0.08683953  0.061526039 0.049417712 -0.035330897 0.158382975 0.2131244


saveRDS(MR_matrix_results,
        "/path/to/project/targetMR/cisMR_res/Negative_MRIVW_res.rds")
MR_matrix_results <- readRDS("/path/to/project/targetMR/cisMR_res/Negative_MRIVW_res.rds")


### [visualization code has been moved to a separate script: forestplot.R]
# The following commented-out block constructs a combined highlight table from MR-PC-GMM results
# (FDR p-value < 0.05, Q p-value > 0.01, IVW p-value < 0.05, consistent direction across robust methods)
# and generates a forest plot combining MRI outcomes, positive controls (CAD/T2D), and negative controls (Tanning).
# See forestplot.R for the complete visualization code.
