# =============================================================================
# Author: Haodong Tian
# Description: Colocalization analysis (coloc.susie) between liver MRI features
#              and plasma proteins, using GWAS summary statistics and in-sample
#              LD matrices.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

### Colocalization analysis for highlighted MRI-protein pairs identified by PWAS

library(coloc)
library(stringr)
library(data.table)

## Chris Wallace's coloc.susie requirements:
# - The same SNPs must be present in both trait 1 (T1) and trait 2 (T2)
# - Include ALL SNPs in the region: do not trim by significance, MAF, or any other criterion
# - The standard deviation of both traits is required (if absent, coloc will estimate it from MAF)

## Additional requirement for coloc.susie (SuSiE-based fine-mapping):
# - An LD matrix is required (one per trait; they can differ since fine-mapping runs separately)
# - The LD matrix must be a square numeric matrix with dimnames corresponding to SNP IDs

## LD matrix strategy:
# - Must use UKB in-sample LD for fine-mapping; reference panel LD is not appropriate here
# - Use LD computed from unrelated EUR UKB-MRI participants
# - PLINK2 is used with ref-first so the REF allele is the dosage reference;
#   REF/ALT are fixed in the .pvar header: #CHROM, POS, ID, REF, ALT
# - Post-hoc validation of dosage allele direction via dbSNP and LDlink is recommended

## Analysis strategy:
# 1. Obtain each protein's full UKB-PPP GWAS
# 2. Define the cis-region (±500 kb window around lead cis-pQTL) and extract cis-UKB-PPP GWAS
# 3. Extract the MRI cis-GWAS for the same region
# 4. QC variants and match UKB-PPP and MRI cis-GWAS to obtain the shared variant list
# 5. Compute in-sample LD matrix (UKB unrelated EUR MRI IDs) using PLINK2
# 6. Prepare list inputs and run coloc / coloc.susie



protein_vector <- unique(plot_dt$Protein)  # 28 highlighted proteins



### 1. Obtain each protein's UKB-PPP GWAS ======================================================================================
outdir <- "/path/to/project/PWAS_Coloc/UKB_PPP_GWAS"


### only run once --------
writeLines(paste0('
import synapseclient, json
syn = synapseclient.Synapse()
syn.login(authToken="<SYNAPSE_TOKEN>")
children = list(syn.getChildren("syn51365303"))  # syn51365303 is the UKB-PPP EUR (discovery) batch ID
with open("', outdir, '/ukbppp_manifest.json", "w") as f:
    json.dump(children, f, indent=2)
print("Done, total files:", len(children))
'), con = paste0(outdir, "/get_manifest.py"))
system(paste0("python3 '", outdir, "/get_manifest.py'"))  # Run Python script to obtain manifest of all UKB-PPP files
### only run once --------


### Download target protein UKB-PPP GWAS files
writeLines(paste0('
import synapseclient, json, re, os
syn = synapseclient.Synapse()
syn.login(authToken="<SYNAPSE_TOKEN>")
with open("', outdir, '/ukbppp_manifest.json") as f:
    children = json.load(f)
targets = ', paste0('["', paste(protein_vector, collapse='","'), '"]'), '
pattern = re.compile(r"^(" + "|".join(targets) + r")_")
matched = [c for c in children if pattern.match(c["name"])]
print("Matched:", len(matched), "files")
for c in matched:
    print(c["name"], c["id"])
    syn.get(c["id"], downloadLocation="', outdir, '")
print("All done.")
'), con = paste0(outdir, "/download_targets.py"))
system(paste0("python3 '", outdir, "/download_targets.py'"))   # Downloads UKB-PPP GWAS for each protein in protein_vector (~1 min per protein)


length(list.files(outdir, pattern = "\\.tar$"))  # 28

# UKB-PPP GWAS format (REGENIE output):
# ALLELE1 is the effect allele (see: https://github.com/rgcgithub/regenie/issues/50)
# CHR BP coordinates are hg38




### 2. Define cis-region and extract cis-UKB-PPP GWAS ===========================================================
# UKB-PPP uses dynamic test regions based on significant SNP positions and LD structure.
# We instead use a fixed ±500 kb window around each lead cis-pQTL for simplicity.

ST16_selected_protein <- ST16[ST16$Protein %in% protein_vector, ]
dim(ST16_selected_protein)  # 125 22
names(ST16_selected_protein)
ST16_selected_protein[, my_window := {tmp <- as.integer(tstrsplit(`Variant ID`, ":")[[2]]); paste(tstrsplit(`Variant ID`, ":")[[1]][1L], floor(min(tmp)-5e5), ceiling(max(tmp)+5e5), sep="_")}, by=Protein]
ST16_selected_protein[, window_length := as.integer(tstrsplit(my_window, "_")[[3]]) - as.integer(tstrsplit(my_window, "_")[[2]]), by=Protein]



library(data.table)
library(stringr)

in_dir  <- "/path/to/project/PWAS_Coloc/UKB_PPP_GWAS"   # directory containing 28 UKB-PPP GWAS tar files
out_dir <- "/path/to/project/PWAS_Coloc/UKB_PPP_cisGWAS" # output directory for cis-GWAS files
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

prot_region_dt <- unique(ST16_selected_protein[, c("Protein","my_window")])  # per-protein cis-window
if (nrow(prot_region_dt) != length(protein_vector)) { print('WARNING: cis-region is not unique for some proteins!') }

setDT(prot_region_dt)
tar_files <- list.files(in_dir, pattern = "\\.tar$", full.names = TRUE)

for (i in 1:nrow(prot_region_dt)) {
  prot   <- prot_region_dt$Protein[i]
  region <- prot_region_dt$my_window[i]
  cat("Processing:", prot, "| region:", region, "\n")

  ss      <- strsplit(region, "_", fixed = TRUE)[[1]]
  chr_now <- ss[1]; bp_start <- as.numeric(ss[2]); bp_end <- as.numeric(ss[3])

  tar_hit <- tar_files[grepl(paste0("^", prot, "_.*\\.tar$"), basename(tar_files))]
  if (length(tar_hit) == 0) { cat("  No tar found for", prot, "\n"); next }
  if (length(tar_hit) > 1)  { cat("  Multiple tar found for", prot, "-> use first\n") }
  tar_hit <- tar_hit[1]

  members   <- system2("tar", c("-tf", shQuote(tar_hit)), stdout = TRUE)
  member_hit <- members[grepl(paste0("/discovery_chr", chr_now, "_"), members)]
  if (length(member_hit) == 0) { cat("  No chr", chr_now, "file in", basename(tar_hit), "\n"); next }
  member_hit <- member_hit[1]

  cmd    <- paste0("tar -xOf ", shQuote(tar_hit), " ", shQuote(member_hit), " | gzip -dc")
  gwas_dt <- fread(cmd = cmd)

  # UKB-PPP GENPOS is hg38; reconstruct hg19 BP from the ID field
  gwas_dt[, GENPOS := as.integer(tstrsplit(ID, ":", fixed = TRUE)[[2]])]
  cis_dt <- gwas_dt[GENPOS >= bp_start & GENPOS <= bp_end]

  out_file <- file.path(out_dir, paste0(prot, "_", region, "_cisGWAS.txt"))
  fwrite(cis_dt, out_file, sep = "\t")
  cat("  Saved:", basename(out_file), "| nrow =", nrow(cis_dt), "\n")
}



### 3. Extract MRI cis-GWAS ==============================================================================================

# QC-filtered (in-sample MAF>0.01; INFO>0.8) MRI GWAS stored on server
mri_gwas_dir <- "/path/to/server/project/GWAS_regenie/Step2_new"  # chr-specific GWAS files for faster I/O
files    <- list.files(mri_gwas_dir, pattern = "\\.regenie$")
mri_names <- gsub("^liver_chr[0-9]+_", "", gsub("\\.regenie$", "", files))
length(unique(mri_names))  # should be 200 unique MRI features

out_dir <- "/path/to/project/PWAS_Coloc/MRI_cisGWAS"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sig_dt <- unique(merge(plot_dt[fdr < 0.05, .(Protein, MRI)], unique(ST16_selected_protein[, c("Protein","my_window")]), by = "Protein", all.x = TRUE))
setDT(sig_dt)

for (i in 1:nrow(sig_dt)) {
  prot <- sig_dt$Protein[i]; mri <- sig_dt$MRI[i]; region <- sig_dt$my_window[i]
  ss      <- strsplit(region, "_", fixed = TRUE)[[1]]
  chr_now <- as.integer(ss[1]); bp_start <- as.numeric(ss[2]); bp_end <- as.numeric(ss[3])

  hit <- list.files(mri_gwas_dir, full.names = TRUE, pattern = paste0("liver_chr", chr_now, "_", mri, "\\.regenie$"))
  if (length(hit) == 0) { cat("MISSING:", prot, as.character(mri), "\n"); next }

  prot_dir <- file.path(out_dir, prot); dir.create(prot_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(prot_dir, paste0(mri, "_", region, "_cisGWAS.txt"))
  if (file.exists(out_file)) { cat("SKIP (exists):", prot, as.character(mri), "\n"); next }

  x     <- fread(hit[1])  # ~20 seconds per file
  cis_x <- x[CHROM == chr_now & GENPOS >= bp_start & GENPOS <= bp_end]

  fwrite(cis_x, out_file, sep = "\t")
  cat("Saved:", prot, "|", as.character(mri), "| nrow =", nrow(cis_x), "\n")
  rm(x); gc()
}  # ~20 mins total




### 4. QC variants and match UKB-PPP and MRI cis-GWAS to obtain shared variant lists ===================================================

## QC filters applied to MRI GWAS: in-sample MAF>0.01, INFO>0.8 (same as main MRI GWAS workflow)
## Broad UKB imputed data uses hg19 coordinates

## For each protein:
# => shared variants table: [CHR | BP | rsID (from MRI GWAS) | Broad_ID (e.g. 11:5248232_T_A)]
# => shared variants LD R matrix (must have rownames and colnames)


protein_cis_dir <- "/path/to/project/PWAS_Coloc/UKB_PPP_cisGWAS"
mri_cis_dir     <- "/path/to/project/PWAS_Coloc/MRI_cisGWAS"
out_dir         <- "/path/to/project/PWAS_Coloc/shared_variant_lists"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

prot_vec <- intersect(basename(list.dirs(mri_cis_dir, recursive = FALSE, full.names = FALSE)), unique(ST16_selected_protein$Protein))

for (prot in prot_vec) {
  cat("Processing:", prot, "\n")
  out_file <- file.path(out_dir, paste0(prot, "_shared_variants.txt"))
  if (file.exists(out_file)) { cat("Skipping:", prot, "| already exists\n"); next }

  prot_file <- list.files(protein_cis_dir, pattern = paste0("^", prot, "_.*cisGWAS\\.txt(\\.gz)?$"), full.names = TRUE)
  if (length(prot_file) == 0) { cat("  No protein cisGWAS file\n"); next }
  prot_gwas <- fread(prot_file[1])

  mri_file <- list.files(file.path(mri_cis_dir, prot), pattern = "_cisGWAS\\.txt(\\.gz)?$", full.names = TRUE)
  if (length(mri_file) == 0) { cat("  No MRI cisGWAS file\n"); next }
  mri_gwas <- fread(mri_file[1])  # one MRI file is sufficient since all share the same variant set
  mri_gwas <- mri_gwas[pmin(A1FREQ, 1 - A1FREQ) > 0.01 & INFO > 0.8]  # MRI GWAS QC: in-sample MAF>0.01, INFO>0.8

  # Match by CHR+BP; rsID and EAF taken from MRI GWAS
  shared_dt <- merge(mri_gwas[, .(CHROM, GENPOS, rsID = ID, EAF_MRI_GWAS = A1FREQ, ALLELE0_MRI = ALLELE0, ALLELE1_MRI = ALLELE1)], unique(prot_gwas[, .(CHROM, GENPOS)]), by = c("CHROM","GENPOS"))
  if (nrow(shared_dt) == 0) { cat("  No shared variants\n"); next }

  ### Match Broad variant IDs from the UKB imputation manifest (mfi) file -------------
  chr_now      <- shared_dt$CHROM[1]
  bp_min       <- min(shared_dt$GENPOS)
  bp_max       <- max(shared_dt$GENPOS)
  bp_file      <- tempfile(fileext = ".txt")
  mfi_sub_file <- tempfile(fileext = ".txt")
  fwrite(data.table(GENPOS = unique(shared_dt$GENPOS)), bp_file, col.names = FALSE, sep = "\t")

  mfi_file <- paste0("/path/to/server/ukb/imputed_v3/ukb_mfi_chr", chr_now, "_v3.txt")
  cmd <- sprintf("awk 'NR==FNR{bp[$1]=1; next} ($3<%d){next} ($3>%d){exit} ($3 in bp){print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8}' %s %s > %s",
                 bp_min, bp_max, shQuote(bp_file), shQuote(mfi_file), shQuote(mfi_sub_file))
  cat('Scanning Broad mfi file...\n')
  system(cmd)  # awk scans only the [min, max] BP range for efficiency (~1 min)

  mfi_dt <- fread(mfi_sub_file, header = FALSE)
  setnames(mfi_dt, c("Broad_ID","Broad_rsID","Broad_BP","Broad_REF","Broad_ALT","Broad_MAF","Broad_MINOR","Broad_INFO"))
  mfi_dt[, `:=`(CHROM = chr_now, GENPOS = Broad_BP)]
  mfi_dt[, allele_key    := paste0(pmin(Broad_REF, Broad_ALT), "_", pmax(Broad_REF, Broad_ALT))]
  shared_dt[, allele_key := paste0(pmin(ALLELE0_MRI, ALLELE1_MRI), "_", pmax(ALLELE0_MRI, ALLELE1_MRI))]
  # Match by CHR+BP+unordered allele pair to obtain Broad variant IDs
  shared_dt <- merge(shared_dt, mfi_dt[, .(CHROM, GENPOS, allele_key, Broad_ID, Broad_rsID, Broad_BP, Broad_REF, Broad_ALT, Broad_MAF, Broad_MINOR, Broad_INFO)], by = c("CHROM","GENPOS","allele_key"), all.x = TRUE)
  shared_dt[, allele_key := NULL]
  ### ---------------------------------------------------------------------------

  out_file <- file.path(out_dir, paste0(prot, "_shared_variants.txt"))
  fwrite(shared_dt, out_file, sep = "\t")
  cat("  Saved:", basename(out_file), "| nrow =", nrow(shared_dt), "\n")
}




### 5. Compute in-sample LD matrix ====================================================================================================================

## Use unrelated UKB MRI individual IDs for LD computation
## (protein cohort may include non-EUR ancestry; MRI cohort is restricted to unrelated EUR)

## Post-hoc check: LD correlations should reflect the REF allele as defined in the .pvar file.
## Validated against LDlink for a subset of variant pairs.

# Write UKB MRI unrelated EUR IID list to server
write.table(data.frame(FID = radiomics_sub_rint$IID, IID = radiomics_sub_rint$IID),
            "/path/to/server/project/PWAS_Coloc/UKB_MRI_unrelated_EUR.txt",
            row.names=FALSE, col.names=FALSE, quote=FALSE)
# Submit LD matrix computation job: qsub get_LDmatrix.sh (~2 hrs)
# Output: *.unphased.vcor1  (LD matrix), *.unphased.vcor1.vars (variant order), *.pvar (REF/ALT info)


# Example: read LD matrix for a single protein
protein_dosage  <- fread("/path/to/server/project/PWAS_Coloc/dosage/ABO.raw")        # ~2 mins
dim(protein_dosage)  # 33059 (UKB-MRI unrelated IDs) x 17583 columns (6 metadata + 17577 variants in rsID_Allele form)
protein_variant <- fread("/path/to/server/project/PWAS_Coloc/dosage/ABO.pvar")
protein_variant$dosage_matrix_colname <- colnames(protein_dosage)[-(1:6)]  # the allele appended to rsID in the column name matches the REF allele in pvar
dim(protein_variant)  # 17577 x 5

# Verify dosage allele direction: the allele appended to rsID in column names corresponds to the REF allele
# Validated via dbSNP: sum(protein_dosage$rs11669668_G) / (2*nrow(protein_dosage)) gives AF ~0.966,
# consistent with G being the reference (dbSNP reports AF of A = 0.024).
# Harmonization downstream is performed relative to the LD REF allele.


prefix <- "/path/to/server/project/PWAS_Coloc/LDmatrix/ABO"

# Read pvar for REF/ALT information
pvar <- fread(paste0(prefix, ".pvar"))  # columns: #CHROM POS ID REF ALT
setnames(pvar, 1, "CHROM")
# Read variant order for LD matrix rows/columns
vars <- fread(paste0(prefix, ".unphased.vcor1.vars"), header = FALSE)
setnames(vars, "V1", "ID")
# Read LD matrix
LD <- as.matrix(fread(paste0(prefix, ".unphased.vcor1"), header = FALSE))  # ~1 min
dim(LD)  # e.g. 4500 x 4500
# Sanity checks
stopifnot(nrow(vars) == nrow(LD), nrow(LD) == ncol(LD))
stopifnot(all(vars$ID == pvar$ID))
# Assign row/column names
rownames(LD) <- vars$ID
colnames(LD) <- vars$ID




### 6. Run coloc, SuSiE fine-mapping, and coloc.susie =================================================================================================================================================

## Harmonization strategy:
# - Restrict to variants present in the LD matrix => final variant list
# - Harmonize UKB-PPP and UKB-MRI cis-GWAS effect alleles to the LD REF allele (from pvar)
# - coloc.susie input lists:
#   T1 / T2: $snp, $position, $beta, $varbeta, $type, $MAF, $N, $LD
#   All betas and LD correlations must reference the same effect allele definition

## Strategy for fine-mapping and colocalization:
# - Always report coloc (coloc.abf) regardless of fine-mapping outcome
# - If T2 (MRI) fine-mapping converges with a credible set, also report fine-mapping results
# - If both T1 and T2 SuSiE fine-mapping converge, run coloc.susie


current_protein <- "ABO";   current_MRI <- "firstorder_90Percentile"
current_protein <- 'APOE';  current_MRI <- "glszm_SizeZoneNonUniformity"

current_protein <- "ACTA2";    current_MRI <- "shape_LeastAxisLength"
current_protein <- "DAPP1";    current_MRI <- "gldm_DependenceEntropy"
current_protein <- "NCAN";     current_MRI <- "firstorder_90Percentile"
current_protein <- "NCAN";     current_MRI <- "ngtdm_Busyness"
current_protein <- "SERPINA1"; current_MRI <- "gldm_LargeDependenceLowGrayLevelEmphasis"


all_proteins <- unique(ST16_selected_protein$Protein)

all_results <- list()
for (current_protein in all_proteins) {
  cat(paste0('================================================================\n'))
  cat(paste0(current_protein, '\n'))

  ### Load LD matrix for this protein's cis-region
  prefix <- paste0("/path/to/server/project/PWAS_Coloc/LDmatrix/", current_protein)
  pvar <- fread(paste0(prefix, ".pvar"))
  setnames(pvar, 1, "CHROM")
  vars <- fread(paste0(prefix, ".unphased.vcor1.vars"), header = FALSE)
  setnames(vars, "V1", "ID")
  cat('Reading LD matrix from server...\n')
  LD <- as.matrix(fread(paste0(prefix, ".unphased.vcor1"), header = FALSE))  # ~1 min
  dim(LD)
  stopifnot(nrow(vars) == nrow(LD), nrow(LD) == ncol(LD))
  stopifnot(all(vars$ID == pvar$ID))
  rownames(LD) <- vars$ID; colnames(LD) <- vars$ID

  ### Load UKB-PPP and MRI cis-GWAS
  PPP_MRI_shared_variant <- fread(paste0("/path/to/project/PWAS_Coloc/shared_variant_lists/", current_protein, "_shared_variants.txt"))
  if (!all(PPP_MRI_shared_variant$rsID == pvar$ID)) { cat("Shared variants must align with LD pvar (same IDs and order).\n"); stop() }
  PPP_file <- list.files("/path/to/project/PWAS_Coloc/UKB_PPP_cisGWAS", pattern = paste0("^", current_protein, "_.*_cisGWAS\\.txt$"), full.names = TRUE)[1]


  ### Prepare UKB-PPP dataset
  PPP_raw <- fread(PPP_file)[TEST == "ADD", .(CHROM, GENPOS, ALLELE0_PPP = ALLELE0, ALLELE1_PPP = ALLELE1, BETA_PPP = BETA, SE_PPP = SE, N_PPP = N)]
  PPP_raw[, allele_key := paste0(pmin(ALLELE0_PPP, ALLELE1_PPP), "_", pmax(ALLELE0_PPP, ALLELE1_PPP))]

  # Anchor table: shared variants with LD REF/ALT appended
  anchor <- copy(PPP_MRI_shared_variant)[, .(CHROM, GENPOS, snp = rsID, MAF = Broad_MAF, ALLELE0_MRI, ALLELE1_MRI)]
  anchor[, allele_key := paste0(pmin(ALLELE0_MRI, ALLELE1_MRI), "_", pmax(ALLELE0_MRI, ALLELE1_MRI))]
  anchor[, `:=`(REF = pvar$REF, ALT = pvar$ALT)]

  PPP_dt <- merge(anchor, PPP_raw, by = c("CHROM", "GENPOS", "allele_key"), all.x = TRUE)
  if (anyDuplicated(PPP_dt[, .(CHROM, GENPOS, snp)])) stop("Duplicate variants in UKB-PPP GWAS at same CHR:BP; please deduplicate.")

  ### Harmonize: flip BETA to REF allele direction
  PPP_dt[, beta_ref := fifelse(ALLELE1_PPP == REF, BETA_PPP, fifelse(ALLELE1_PPP == ALT, -BETA_PPP, NA_real_))]
  original_num_variants <- nrow(PPP_dt)

  ## If variant count exceeds 5000, retain the 5000 most significant variants in UKB-PPP GWAS
  # (coloc.susie computational limit)
  if (original_num_variants > 5000) {
    cat('This protein has >5000 variants; retaining the 5000 most significant variants in UKB-PPP GWAS.\n')
    PPP_dt[, p_PPP := 2 * pnorm(-abs(BETA_PPP / SE_PPP))]
    keep_idx <- order(PPP_dt$p_PPP)[1:5000]
    PPP_dt <- PPP_dt[keep_idx]
    LD     <- LD[keep_idx, keep_idx]
    pvar   <- pvar[keep_idx, ]
  }
  T1 <- list(snp = PPP_dt$snp, position = PPP_dt$GENPOS,
             beta = PPP_dt$beta_ref, varbeta = PPP_dt$SE_PPP^2, type = "quant",
             MAF = PPP_dt$MAF, N = median(PPP_dt$N_PPP, na.rm = TRUE), LD = LD)

  ### SuSiE fine-mapping for UKB-PPP (run once per protein, shared across all MRI features)
  susie_T1 <- susie_rss(z = T1$beta/sqrt(T1$varbeta), n = T1$N, R = T1$LD, L = 10)
  cat(paste0('UKB-PPP SuSiE converged: ', susie_T1$converged, '\n'))

  ## Prepare MRI file list for the inner loop
  mri_files <- list.files(file.path("/path/to/project/PWAS_Coloc/MRI_cisGWAS", current_protein), pattern = "_cisGWAS\\.txt$", full.names = FALSE)
  all_MRIs  <- sub("_\\d+_\\d+_\\d+_cisGWAS\\.txt$", "", mri_files)

  all_results[[current_protein]] <- list()
  for (current_MRI in all_MRIs) {
    cat(paste0('---------------------------------\n'))
    cat(paste0(current_MRI, '\n'))

    MRI_file <- list.files(file.path("/path/to/project/PWAS_Coloc/MRI_cisGWAS", current_protein), pattern = paste0("^", current_MRI, "_.*_cisGWAS\\.txt$"), full.names = TRUE)[1]

    MRI_raw <- fread(MRI_file)[TEST == "ADD", .(CHROM, GENPOS, ALLELE0_MRI_GWAS = ALLELE0, ALLELE1_MRI_GWAS = ALLELE1, BETA_MRI = BETA, SE_MRI = SE, N_MRI = N)]

    MRI_dt <- merge(anchor, MRI_raw, by = c("CHROM", "GENPOS"), all.x = TRUE)

    if (original_num_variants > 5000) { MRI_dt <- MRI_dt[keep_idx] }

    if (nrow(MRI_dt) != nrow(PPP_dt)) stop("UKB MRI and PPP datasets must have the same number of rows.")
    if (!all(MRI_dt$GENPOS == PPP_dt$GENPOS)) stop("UKB MRI and PPP datasets must contain identical variants in identical order.")

    ### Harmonize: flip MRI BETA to REF allele direction
    MRI_dt[, beta_ref := fifelse(ALLELE1_MRI_GWAS == REF, BETA_MRI, fifelse(ALLELE1_MRI_GWAS == ALT, -BETA_MRI, NA_real_))]

    if (!all(PPP_dt$snp == pvar$ID) || !all(MRI_dt$snp == pvar$ID)) stop("Final PPP/MRI GWAS variant order does not match pvar.")
    if (anyNA(PPP_dt$beta_ref) || anyNA(MRI_dt$beta_ref)) stop("NA values after harmonization: some variants have alleles inconsistent with LD REF/ALT.")
    if (!all(PPP_dt$REF == pvar$REF) || !all(MRI_dt$REF == pvar$REF)) stop("REF allele mismatch between GWAS and LD matrix pvar.")
    if (!is.null(rownames(LD)) && !all(PPP_dt$snp == rownames(LD))) stop("LD matrix row names do not match final PPP GWAS variant order.")


    ### Prepare coloc.susie input list for T2 (MRI)
    T2 <- list(snp = MRI_dt$snp, position = MRI_dt$GENPOS,
               beta = MRI_dt$beta_ref, varbeta = MRI_dt$SE_MRI^2, type = "quant",
               MAF = MRI_dt$MAF, N = median(MRI_dt$N_MRI, na.rm = TRUE), LD = LD)

    ### Diagnostic Manhattan plot (harmonized region)
    PPP_dir <- "/path/to/project/PWAS_Coloc/UKB_PPP_cisGWAS"
    MRI_dir <- "/path/to/project/PWAS_Coloc/MRI_cisGWAS"
    PPP <- fread(list.files(PPP_dir, paste0("^", current_protein, "_.*_cisGWAS\\.txt$"), full.names=TRUE)[1])[TEST=="ADD", .(CHROM, GENPOS, logp=LOG10P, trait="PPP")]
    MRI <- fread(list.files(file.path(MRI_dir, current_protein), paste0("^", current_MRI, "_.*_cisGWAS\\.txt$"), full.names=TRUE)[1])[TEST=="ADD", .(CHROM, GENPOS, logp=LOG10P, trait="MRI")]
    Manhanttan_ggdata <- rbind(PPP, MRI)
    Manhanttan_plot <- ggplot(Manhanttan_ggdata, aes(GENPOS, logp)) + geom_point(size=0.8) + geom_hline(yintercept=0, linewidth=0.3) + facet_wrap(~trait, ncol=1, scales="free_y") + labs(x=paste0("Chr", unique(PPP$CHROM), " position"), y=expression(-log[10](italic(P))), title=paste(current_protein, "vs", current_MRI)) + theme_bw()
    Manhanttan_plot
    ## Harmonized z-score scatter (UKB-PPP vs UKB-MRI, restricted to the 5000 retained variants if applicable)
    z_ggdata <- merge(PPP_dt[, .(snp, z_PPP=beta_ref/SE_PPP)], MRI_dt[, .(snp, z_MRI=beta_ref/SE_MRI)], by="snp")
    z_plot <- ggplot(z_ggdata, aes(z_PPP, z_MRI)) + geom_point(size=0.8) + geom_smooth(method="lm", se=FALSE) + labs(x=paste0(current_protein, " Z"), y=paste0(current_MRI, " Z"), title=paste(current_protein, "vs", current_MRI)) + theme_bw()
    z_plot


    ### Colocalization assuming a single causal variant (coloc.abf) --------------------------------
    coloc_res <- coloc.abf(T1, T2)  # fast; does not require LD matrix
    coloc_res$summary

    ### SuSiE fine-mapping for UKB-MRI -------------------------------------------------------
    susie_T2 <- susie_rss(z = T2$beta/sqrt(T2$varbeta), n = T2$N, R = T2$LD, L = 10)
    if (is.null(susie_T2$sets$cs)) { T2_signal <- FALSE } else { T2_signal <- TRUE }

    ### Store summary results
    results_summary <- c(current_MRI, coloc_res$summary,
                         susie_T2$converged, T2_signal, susie_T1$converged)
    names(results_summary) <- c('MRI_name', 'num_of_variants', paste0('PP_H', 0:4),
                                'MRI_susie_coverge', 'MRI_susie_signal', 'PPP_susie_coverge')

    ### coloc.susie (run only when both SuSiE models converge and MRI has a credible set) -----
    if (T2_signal & susie_T2$converged & susie_T1$converged) {
      cat('Running coloc.susie...\n')
      coloc_susie_res  <- coloc.susie(T1, T2)
      Coloc_SuSiE_res  <- coloc_susie_res$summary
    } else {
      Coloc_SuSiE_res <- NULL
    }


    ### Store all results for the current MRI feature
    all_results[[current_protein]][[current_MRI]] <- list(
      Manhanttan_plot  = Manhanttan_plot,
      Manhanttan_ggdata = Manhanttan_ggdata,
      z_plot           = z_plot,
      z_ggdata         = z_ggdata,           # harmonized z-scores: UKB-PPP vs UKB-MRI over the cis-region
      susie_UKB_MRI    = susie_T2,           # UKB-MRI SuSiE results
      susie_UKB_PPP    = susie_T1,           # UKB-PPP SuSiE results
      results_summary  = results_summary,    # coloc PP_H0-H4 and SuSiE convergence indicators
      Coloc_SuSiE_res  = Coloc_SuSiE_res    # coloc.susie result (if applicable)
    )
    cat('\n')

  }
  # End MRI loop for current protein

}

saveRDS(all_results, file = "/path/to/project/PWAS_Coloc/all_results.rds")




### Supplementary table ========================================================

library(data.table)

supp_tab <- rbindlist(lapply(names(all_results), function(prot) {
  rbindlist(lapply(names(all_results[[prot]]), function(mri) {
    x <- all_results[[prot]][[mri]]$results_summary
    x <- as.list(x)
    data.table(Protein = prot, MRI = mri)[, c(names(x)) := x]
  }), fill = TRUE)
}), fill = TRUE)

setnames(supp_tab,
         old = c("MRI_name", "MRI_susie_coverge", "MRI_susie_signal", "PPP_susie_coverge"),
         new = c("MRI_name_in_vector", "MRI_converge", "MRI_signal", "PPP_converge"),
         skip_absent = TRUE)

target_cols <- c("Protein", "MRI", "num_of_variants", paste0("PP_H", 0:4), "MRI_converge", "MRI_signal", "PPP_converge")
target_cols <- intersect(target_cols, names(supp_tab))
supp_tab    <- supp_tab[, ..target_cols]

num_cols <- setdiff(names(supp_tab), c("Protein", "MRI","MRI_converge", "MRI_signal", "PPP_converge"))
supp_tab[, (num_cols) := lapply(.SD, function(z) suppressWarnings(as.numeric(as.character(z)))), .SDcols = num_cols]

View(supp_tab)
sum(supp_tab$PP_H4 > 0.8) / nrow(supp_tab)      # fraction of pairs with coloc.abf evidence (PP_H4>0.8)
unique(supp_tab$Protein[supp_tab$PP_H4 > 0.8])   # proteins with at least one colocalized MRI feature

### Append coloc.susie H4 for pairs where both SuSiE models converged and MRI has a credible set
supp_tab$coloc_susie_H4 <- NA
supp_tab[MRI_converge==TRUE & MRI_signal==TRUE & PPP_converge==TRUE,
         coloc_susie_H4 := mapply(function(prot, mri) max(all_results[[prot]][[mri]]$Coloc_SuSiE_res$PP.H4.abf, na.rm=TRUE), Protein, MRI)]

sum(supp_tab$coloc_susie_H4 > 0.8, na.rm = TRUE) / nrow(supp_tab)
# Note: when the single-causal-variant assumption is violated, coloc.abf may have both
# reduced power and inflated type I error; coloc.susie is preferred in those cases.


### Final H4: use coloc.susie H4 if available, otherwise fall back to coloc.abf H4
supp_tab$final_H4 <- ifelse(!is.na(supp_tab$coloc_susie_H4), supp_tab$coloc_susie_H4, supp_tab$PP_H4)
sum(supp_tab$final_H4 > 0.8) / nrow(supp_tab)

saveRDS(supp_tab, file = "/path/to/project/PWAS_Coloc/Coloc_supp_tab.rds")



length(unique(supp_tab$Protein[supp_tab$final_H4 > 0.8]))  # proteins with colocalization evidence
length(unique(supp_tab$MRI[supp_tab$final_H4 > 0.8]))      # MRI features with colocalization evidence




### Visualization: chord diagram of colocalized protein-MRI pairs ======================================================================

library(data.table)
library(circlize)

dt_chord <- supp_tab[final_H4 > 0.8, .(Protein, MRI)]

prot_annot <- data.table(Protein = unique(dt_chord$Protein))
prot_annot[, cat := fcase(
  Protein %in% c("APOE","NCAN","ADH4","FABP2","GSTA1","GSTA3","BCAT1","SERPINA1","SERPINA6","TREH","CELA2A"), "Lipid / Metabolic",
  Protein %in% c("ACTA2","INHBB","INHBC","EFEMP1","PDE5A"), "Fibrosis / ECM",
  Protein %in% c("AIF1","FOLR2","CASP9","TNFSF10","DAPP1","ARHGAP25"), "Immune / Inflammation",
  default = "Other"
)]
prot_annot[, cat := factor(cat, levels = c("Lipid / Metabolic","Fibrosis / ECM","Immune / Inflammation","Other"))]
prot_annot <- prot_annot[order(cat, Protein)]

proteins     <- prot_annot$Protein
mris         <- dt_chord[, .N, by = MRI][order(-N, MRI)]$MRI
sector_order <- c(proteins, mris)

# Assign distinct colors within each biological category
make_pal <- function(cols, n) {
  if (n == 1) cols[2] else grDevices::colorRampPalette(cols)(n)
}

prot_annot[cat == "Lipid / Metabolic",       col := make_pal(c("#F6D58A", "#C97B00", "#7A4A00"), .N)]
prot_annot[cat == "Fibrosis / ECM",          col := make_pal(c("#B8D8BE", "#4E7F5A", "#21442B"), .N)]
prot_annot[cat == "Immune / Inflammation",   col := make_pal(c("#B9CBEA", "#4A6FA5", "#1F3E6E"), .N)]
prot_annot[cat == "Other",                   col := make_pal(c("#B8B8B8", "#7F7F7F", "#3F3F3F"), .N)]

protein_cols <- setNames(prot_annot$col, prot_annot$Protein)
base_cols    <- c("Lipid / Metabolic"="#C97B00", "Fibrosis / ECM"="#4E7F5A", "Immune / Inflammation"="#4A6FA5", "Other"="#7F7F7F")
grid.col     <- c(protein_cols, setNames(rep("grey88", length(mris)), mris))

mat <- table(factor(dt_chord$Protein, levels = proteins), factor(dt_chord$MRI, levels = mris))

circos.clear()
circos.par(start.degree = 270, canvas.xlim = c(-1.3, 1.3), canvas.ylim = c(-1.3, 1.3))
chordDiagram(mat, order = sector_order, grid.col = grid.col, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sn      <- CELL_META$sector.index
  lab_col <- if (sn %in% proteins) protein_cols[sn] else "black"
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], sn, facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 0.6, col = lab_col)
}, bg.border = NA)
legend("topleft", legend = names(base_cols), fill = base_cols, bty = "n", cex = 0.8, inset = c(0, 0.3))
