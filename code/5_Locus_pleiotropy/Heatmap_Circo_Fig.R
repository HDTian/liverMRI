# =============================================================================
# Author: Haodong Tian
# Description: Generates genetic correlation heatmaps (within liver MRI features
#              and between MRI features and GenomicSEM factors) and Circo-style
#              GWAS figures for the paper.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================


### Genetic correlation (rg) heatmap and Circo GWAS figure


### (1) Genetic correlation heatmap: within liver MRI features and GenomicSEM factors --------
### --------------------------------------------------------------------------------------------

# Load the GenomicSEM LDSC output (genetic covariance/correlation matrix, Chr 1-22)
e <- new.env()
load("/path/to/project/GenomicSEM/LDSCoutput_cate_liverMRI_ALL.RData", envir = e)
GenomicSEM_LDSCoutput_ALL <- e$LDSCoutput_ALL

# S_Stand: 30x30 standardized genetic covariance (correlation) matrix (30 features considered initially)
# V_Stand: 465x465 sampling covariance matrix (465 = 30*29/2 + 30 unique elements)
dim(GenomicSEM_LDSCoutput_ALL$S_Stand) # 30 30
dim(GenomicSEM_LDSCoutput_ALL$V_Stand) # 465 465

rg_matrix <- GenomicSEM_LDSCoutput_ALL$S_Stand

# Extract SE matrix from the diagonal of V_Stand (lower triangle including diagonal)
rg_variance_matrix <- matrix(NA_real_, 30, 30)
rg_variance_matrix[lower.tri(rg_variance_matrix, diag = TRUE)] <- diag(GenomicSEM_LDSCoutput_ALL$V_Stand)
colnames(rg_variance_matrix) <- colnames(rg_matrix)
rownames(rg_variance_matrix) <- colnames(rg_matrix)
rownames(rg_matrix) <- colnames(rg_matrix)

# Verification example: cross-trait intercept and genetic covariance for firstorder_90 vs firstorder_TotalEnergy
# Mean Z*Z: -0.0118
# Cross trait Intercept: 0.0143 (0.0053)
# Total Observed Scale Genetic Covariance (g_cov): -0.0353 (0.0126)
# g_cov Z: -2.81
# g_cov P-value: 0.0049776
-0.193493691 / sqrt(0.004746700)  # -2.808477 ≈ -2.81 (confirms SE extraction is correct)


### Subset to the 10 GenomicSEM lead features (selected by |z-score of SNP-heritability| threshold)
GenomicSEM_names <- readRDS("/path/to/project/GenomicSEM/final_GenomicSEM_names.rds")

# Features ordered by GenomicSEM factor membership (F1 through F5)
GenomicSEM_names <- c(
  'firstorder_TotalEnergy', 'glszm_GrayLevelNonUniformity_inp', 'shape_LeastAxisLength',
  'firstorder_TotalEnergy_inp', 'glrlm_RunLengthNonUniformity',  # F1
  'shape_MinorAxisLength',                                         # F5
  'glszm_ZoneVariance', 'shape_SurfaceVolumeRatio',               # F2
  'shape_Sphericity',                                              # F3
  'liver_fat'                                                      # F4
)

rg_matrix_sub        <- rg_matrix[GenomicSEM_names, GenomicSEM_names]
rg_variance_matrix_sub <- rg_variance_matrix[GenomicSEM_names, GenomicSEM_names]
rg_se_matrix_sub     <- sqrt(rg_variance_matrix_sub)

# Symmetrize SE matrix (fill upper triangle from lower triangle)
M <- rg_se_matrix_sub; M[is.na(M)] <- t(M)[is.na(M)]; rg_se_matrix_sub <- M

dim(rg_matrix_sub); dim(rg_se_matrix_sub)  # 10 10

# Build combined rg (SE) matrix for supplementary table
rg_combined <- matrix(
  paste0(round(rg_matrix_sub, 3), " (", round(rg_se_matrix_sub, 3), ")"),
  nrow = 10,
  dimnames = dimnames(rg_matrix_sub)
)
saveRDS(rg_combined, file = "/path/to/project/GWAS_res/rg_combined.rds")


### Corrplot heatmap (lower triangle + Bonferroni significance asterisks) --------------------
library(corrplot)

rg <- rg_matrix_sub; se <- rg_se_matrix_sub
z  <- rg / se; p <- 2 * pnorm(-abs(z))  # Z-scores and two-sided P-values

# Bonferroni correction for unique off-diagonal pairs
m    <- nrow(rg) * (nrow(rg) - 1) / 2
bonf <- 0.05 / m

# Set diagonal to neutral values (no self-correlation highlight)
diag(p)  <- 1
diag(rg) <- 0

corrplot(rg,
         method      = "color",       # color-block heatmap
         type        = "full",        # show full matrix
         diag        = TRUE,
         tl.col      = "black",
         tl.cex      = 0.8,
         tl.pos      = NULL,
         cl.pos      = 'b',
         col.lim     = c(-0.9, 0.9),
         p.mat       = p,             # P-value matrix for significance marking
         sig.level   = bonf,          # Bonferroni-corrected significance threshold
         insig       = "label_sig",   # mark significant cells with a symbol
         pch         = "*",
         pch.cex     = 1.2,
         pch.col     = "black",
         addgrid.col = "black"
)

# Helper: plot feature names broken at underscores for compact display
x <- GenomicSEM_names
n <- length(x)
plot.new()
par(mar = c(0, 0, 0, 0))
plot.window(xlim = c(0, 1), ylim = c(0, 1))
y   <- seq(0.9, 0.1, length.out = n)
lab <- gsub("_", "\n", x, fixed = TRUE)
cex <- pmax(0.65, pmin(1.0, 26 / nchar(x)))
for (i in seq_len(n)) text(0.5, y[i], lab[i], cex = cex[i], adj = c(0.5, 0.5))




### (2) Fuji-Circos GWAS figure --------------------------------------------------------------
### ------------------------------------------------------------------------------------------


### About the Fuji-Circos pipeline:
# Circos is a Perl-based visualization tool. The fujiplot project wraps Circos
# for multi-trait GWAS locus visualization.
# Reference: https://github.com/mkanai/fujiplot
# Command syntax: Rscript fujiplot.R [input.txt] [traitlist.txt] [output_dir]

### Installation notes (macOS via Homebrew):
# Step 1: Install Circos and its Perl module dependencies via Homebrew and cpanm
#   brew install gd
#   brew install cpanminus
#   cpanm Config::General Font::TTF::Font Math::Bezier Math::VecStat Readonly Set::IntSpan Text::Format
#   cpanm --force GD::Polyline
#   brew install circos
#   circos -modules   # verify all required modules are installed

### Step 2: Download fujiplot from GitHub: https://github.com/mkanai/fujiplot/tree/master

### Step 3: Prepare the two key input files — [input.txt] and [traitlist.txt] — (see below)

### Step 4: Run my_fujiplot.R interactively in RStudio rather than via Terminal
# (my_fujiplot.R is a modified version of fujiplot.R with corrected paths and
#  Circos command adjustments to prevent errors)
# Key customizations applied:
#   - chromosome.conf / ideogram.conf: restricted to hs1-hs22 (exclude X and Y)
#   - Pleiotropy defined as any locus appearing in >1 TRAIT (not >1 CATEGORY)
#   - Gene labels drawn for all lead SNPs, with deduplication for visual clarity
#   - LARGE_POINT_SIZE increased from 16 to 24 for pleiotropic loci
#   - data_tracks.conf: auxiliary line color set to vlgrey

### Step 5: Post-processing SVG for font control
# Circos font handling in PNG/SVG/PDF is limited via .conf files alone.
# Direct SVG text manipulation via Perl one-liners is more reliable:
#   # Set all text to CMUBright-Roman italic
#   perl -0777 -i -pe 's/font-family="[^"]*"/font-family="CMUBright-Roman"/g; s/style="([^"]*)"/style="$1;font-style:italic"/g' circos.svg
#   # Remove italic from chr 1-22 labels
#   perl -0777 -i -pe 's{(<text\b[^>]*>\s*\d+\s*</text>)}{ my $t=$1; $t=~s/;?font-style:italic//g; $t }gse' circos.svg
#   # Reduce chr label font size from 62.4px to 48px
#   perl -0777 -i -pe 's/(<text\b[^>]*\bfont-size=")62\.4px("(?=[^>]*>\s*\d+\s*<\/text>))/${1}48px$2/g' circos.svg
#   # Reduce gene label font size from 31.2px to 24px
#   perl -0777 -i -pe 's/(<text\b[^>]*\bfont-size=")31\.2px(")/${1}24px$2/g' circos.svg
# Convert final SVG to PDF via Inkscape:
#   /Applications/Inkscape.app/Contents/MacOS/inkscape circos.svg --export-type=pdf --export-filename circos.pdf


### About the two input files for Fuji-Circos:
# [input.txt] — one row per significant trait-locus association:
#   LOCUS_ID  CATEGORY  TRAIT  CHR  BP(hg19)  MARKER(rsID)  GENE
#   LOCUS_ID is an integer identifier; loci sharing LD (r2 > 0.2, within 1 Mb) share the same ID.
#   Pleiotropic loci (LOCUS_ID appearing for >1 TRAIT) are highlighted in the plot.
#
# [traitlist.txt] — defines trait grouping order and colors:
#   CATEGORY  TRAIT  COLOR


### Prepare [input.txt] -----------------------------------------------------------------------

GenomicSEM_names  # 10 features; only genome-wide significant (5e-8) lead SNPs after LD-clumping

# GWAS summary statistics sources:
#   Liver radiomics features: /path/to/project/new_UVMR/retained_UVMR_summary_1e6_new
#     Columns: CHR  BP  SNP  ALLELE1  A1FREQ  BETA  SE  P
#   Liver fat: /path/to/project/targetMR/liver_fat_leadGWAS.tsv
#     Columns: variant_id  chromosome  base_pair_location  effect_allele  other_allele
#              effect_allele_frequency  beta  standard_error  P_BOLT_LMM_INF
# Note: input files use p < 1e-6 threshold; here we further filter to p < 5e-8 for the Fuji plot.

# Note on LD across merged rsIDs: fujiplot handles pleiotropy detection internally via
# shared LOCUS_ID; no manual LD check across traits is needed.


### Helper function 1: map rsID list to nearest gene within ±0.1 Mb (Ensembl BioMart)
library(biomaRt)
library(data.table)

get_closest_gene <- function(rsID_vector) {

  ## Step 1: retrieve genomic coordinates (GRCh38) for each rsID via Ensembl SNP mart
  m_snp <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp", mirror = "useast")
  v <- as.data.table(getBM(
    attributes = c("refsnp_id", "chr_name", "chrom_start"),
    filters    = "snp_filter", values = rsID_vector, mart = m_snp
  ))
  setnames(v, c("refsnp_id", "chr_name", "chrom_start"), c("rsID", "chr", "pos"))
  v <- v[chr %in% c(as.character(1:22), "X", "Y"), unique(.SD), by = rsID]

  ## Step 2: query gene annotations within ±0.1 Mb windows around each rsID
  win <- 1e5  # 100 kb window
  v[, `:=`(start = pmax(1L, pos - win), end = pos + win)]
  v[, region := paste0(chr, ":", start, ":", end)]
  m_gene <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  g <- as.data.table(getBM(
    attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
    filters    = "chromosomal_region", values = v$region, mart = m_gene
  ))
  setnames(g, c("chromosome_name", "start_position", "end_position", "external_gene_name"),
           c("chr", "gstart", "gend", "gene"))

  ## Step 3: for each rsID, find the nearest named gene
  res <- rbindlist(lapply(seq_len(nrow(v)), function(i) {
    chr <- v$chr[i]; pos <- v$pos[i]; rs <- v$rsID[i]
    gg  <- as.data.frame(g)
    gg  <- gg[gg$chr == chr, ]
    gg  <- gg[(abs(gg$gstart - pos) < win) | (abs(gg$gend - pos) < win), ]
    gg  <- gg[!is.na(gg$gene) & gg$gene != "", ]  # retain only genes with an HGNC-style name
    if (!nrow(gg)) return(data.table(rsID = rs, nearest_gene = NA_character_, distance_bp = NA_real_))
    # distance = 0 if rsID falls within a gene body; otherwise minimum flanking distance
    d <- ifelse(pos < gg$gstart, gg$gstart - pos, ifelse(pos > gg$gend, pos - gg$gend, 0))
    j <- which.min(d)
    data.table(rsID = rs, nearest_gene = gg$gene[j], distance_bp = d[j])
  }))
  return(res)
}

# Example usage
rsID_vector <- c("rs58542926", "rs429358", "rs738409")
gene_res <- get_closest_gene(rsID_vector)
gene_res  # output: rsID | nearest_gene | distance_bp (distance_bp = 0 means rsID is within the gene body)


### Helper function 2: standardize GWAS column names across files with varying headers
library(data.table)

std_gwas_cols <- function(dt) {
  nm <- toupper(names(dt))
  names(dt) <- nm
  get1 <- function(cands) {
    hit <- intersect(cands, nm)
    if (length(hit) == 0) stop("Missing column: ", paste(cands, collapse = ", "))
    hit[1]
  }
  chr_col <- get1(c("CHR", "CHROM", "CHROMOSOME"))
  pos_col <- get1(c("BP", "POS", "BASE_PAIR_LOCATION"))
  snp_col <- get1(c("SNP", "RSID", "VARIANT_ID"))
  p_col   <- get1(c("P", "PVALUE", "P_BOLT_LMM_INF"))
  dt[, .(CHR = get(chr_col), BP = get(pos_col), SNP = get(snp_col), P = get(p_col))]
}


### Collect and merge GWAS results across all 10 features
path_rad <- "/path/to/project/new_UVMR/retained_UVMR_summary_1e6_new"
path_lf  <- "/path/to/project/targetMR/liver_fat_leadGWAS.tsv"

res_list <- lapply(GenomicSEM_names, function(feat) {
  message("Processing: ", feat)
  if (feat == "liver_fat") {
    dt <- fread(path_lf)
  } else {
    f <- file.path(path_rad, paste0(feat, ".tsv"))
    if (!file.exists(f)) stop("File not found: ", f)
    dt <- fread(f)
  }
  out <- std_gwas_cols(dt)
  out[, feature := feat]
  out
})

GWAS_sub <- rbindlist(res_list)

# Filter to genome-wide significance threshold
GWAS_5e8 <- GWAS_sub[GWAS_sub$P < 5 * 10^(-8), ]
dim(GWAS_5e8)              # 133 5
length(unique(GWAS_5e8$SNP))  # 85 unique rsIDs

# Annotate each significant rsID with its nearest gene
gene_res      <- get_closest_gene(GWAS_5e8$SNP)
dim(gene_res)              # 85 3
gene_res_noNA <- na.omit(gene_res)

GWAS_5e8_noNA <- GWAS_5e8[GWAS_5e8$SNP %in% gene_res_noNA$rsID, ]
dim(GWAS_5e8_noNA)         # 131 5 (2 rsIDs without gene annotation removed)
GWAS_5e8_noNA$GENE <- gene_res_noNA$nearest_gene[match(GWAS_5e8_noNA$SNP, gene_res_noNA$rsID)]

# Assign GenomicSEM factor category (F1-F5) to each feature
library(data.table)

map <- data.table(
  feature = c(
    "firstorder_TotalEnergy",             # F1
    "glszm_GrayLevelNonUniformity_inp",   # F1
    "shape_LeastAxisLength",              # F1
    "firstorder_TotalEnergy_inp",         # F1
    "glrlm_RunLengthNonUniformity",       # F1
    "shape_MinorAxisLength",              # F5
    "glszm_ZoneVariance",                 # F2
    "shape_SurfaceVolumeRatio",           # F2
    "shape_Sphericity",                   # F3
    "liver_fat"                           # F4
  ),
  category = c("F1", "F1", "F1", "F1", "F1", "F5", "F2", "F2", "F3", "F4")
)

setDT(GWAS_5e8_noNA)
GWAS_5e8_noNA[map, category := i.category, on = "feature"]


### Assign LOCUS_ID via PLINK LD clumping (r2 > 0.2 within 1 Mb defines a shared locus)

# Write dummy summary statistics for PLINK input (P must be < 0.05 for SP2 output)
dummy_sumstats <- data.table(SNP = unique(GWAS_5e8_noNA$MARKER), P = 0.001)
dummy_file  <- '/path/to/project/GWAS_res/PLINKoutput/dummy_sumstats.txt'
fwrite(dummy_sumstats, dummy_file, sep = "\t", quote = FALSE)

pathname    <- "/path/to/project"
plink_path  <- file.path(pathname, "PheWAS/plink_mac_20250615/plink")
bfile_path  <- file.path(pathname, "PheWAS/1000G_QC")
out_prefix  <- '/path/to/project/GWAS_res/PLINKoutput/locusID'

plink_cmd <- sprintf(
  '"%s" --bfile "%s" --clump "%s" --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 1000 --out "%s"',
  plink_path, bfile_path, dummy_file, out_prefix
)
system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
# Result: 57 clumps formed from 83 top variants

# Parse PLINK .clumped output: expand each clump (lead + secondary SNPs) into a MARKER -> LOCUS_ID map
clumped   <- fread(paste0(out_prefix, ".clumped"))
locus_map <- clumped[, {
  sp2  <- unlist(strsplit(SP2, ","))
  sp2  <- sp2[sp2 != "NONE"]
  sp2  <- sub("\\(.*\\)$", "", sp2)  # remove allele-frequency suffixes e.g. "(1)"
  snps <- unique(c(SNP, sp2))
  data.table(MARKER = snps)
}, by = seq_len(nrow(clumped))]
setnames(locus_map, "seq_len", "LOCUS_ID")

GWAS_5e8_noNA <- merge(GWAS_5e8_noNA, locus_map, by = "MARKER", all.x = TRUE)


### Finalize input.txt for fujiplot
GWAS_5e8_final <- GWAS_5e8_noNA
setnames(GWAS_5e8_final,
         old = c("category", "feature", "SNP"),
         new = c("CATEGORY", "TRAIT",   "MARKER"))
GWAS_5e8_final <- GWAS_5e8_final[, .(LOCUS_ID, CATEGORY, TRAIT, CHR, BP, MARKER, GENE)]

dim(GWAS_5e8_final)  # 131 7
fwrite(GWAS_5e8_final,
       "/path/to/project/GWAS_res/input.txt",
       sep = "\t", quote = FALSE, na = "NA")

# Number of pleiotropic loci (LOCUS_ID appearing across >1 trait)
sum(table(GWAS_5e8_final$LOCUS_ID) > 1)  # 21 pleiotropic loci


### Prepare [traitlist.txt] -------------------------------------------------------------------

# Derive unique CATEGORY-TRAIT pairs in factor-order (F1-F5), then assign colors
traitlist <- unique(GWAS_5e8_final[, .(CATEGORY, TRAIT)], by = c("CATEGORY", "TRAIT"))
traitlist[, ord := match(TRAIT, map$feature)]
setorder(traitlist, ord); traitlist$ord <- NULL

cat_col           <- unique(traitlist[, .(CATEGORY)])
cat_col$CATEGORY  <- paste0('F', 1:5)
cat_col[, COLOR   := c("lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk")]  # F1 to F5

traitlist[, COLOR := cat_col$COLOR[match(CATEGORY, cat_col$CATEGORY)]]

dim(traitlist)  # 10 3
fwrite(traitlist,
       "/path/to/project/GWAS_res/traitlist.txt",
       sep = "\t", quote = FALSE, na = "NA")
