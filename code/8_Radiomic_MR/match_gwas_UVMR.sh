#!/usr/bin/env bash
# =============================================================================
# Author: Haodong Tian
# Description: Shell script to match and align GWAS summary statistics for
#              two-sample UVMR analysis.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================
set -euo pipefail

# This script matches GWAS summary statistics for all UVMR IV variants across
# multiple outcome traits. For each trait, it attempts to find a direct rsID
# match and, for missing SNPs, queries proxy SNPs via LDlink. The output is a
# set of per-trait TSV files containing harmonised summary data ready for UVMR.

# Output: per-trait TSV files under $OUT_DIR, one per outcome GWAS.
# Each file contains the same set of SNPs (from UVMRcoords_new.txt, i.e. the
# 59 retained feature lead variants), facilitating downstream UVMR analyses.

# --- Path configuration ---
PHEWAS_DIR="/path/to/project/PheWAS"
COORDS="/path/to/project/new_UVMR/UVMRcoords_new.txt"  # lead SNPs for 59 retained features
TRAIT_INFO="/path/to/project/trait_info.xlsx"
OUT_DIR="/path/to/project/new_UVMR"
# --------------------------

mkdir -p "$OUT_DIR"
export PHEWAS_DIR COORDS TRAIT_INFO OUT_DIR

echo "[INFO] Reading SNP coordinate list..." >&2
echo "[INFO] Total rsIDs: $(( $(wc -l < "$COORDS") ))" >&2

Rscript - <<'EOF'
library(data.table)
library(readxl)

# Set working directory to PHEWAS_DIR so that helper scripts are accessible
setwd(Sys.getenv("PHEWAS_DIR"))

# Detect whether the coords file has a header by inspecting the third field of the first line
first <- readLines(Sys.getenv("COORDS"), n = 1)
third <- strsplit(first, "[ \t]+")[[1]][3]
skip  <- if (grepl("^rs", third, ignore.case = TRUE)) 0 else 1
coords <- fread(Sys.getenv("COORDS"), skip = skip, header = FALSE)[[3]]  # load column 3 (rsIDs)

# Read trait metadata; skip rows with missing file names
info <- read_excel(Sys.getenv("TRAIT_INFO"))
info <- info[!is.na(info[["file name"]]) & nzchar(info[["file name"]]), ]
cat(sprintf("[R] %d trait files to process.\n", nrow(info)))

for (i in seq_len(nrow(info))) {
  fn       <- info[["file name"]][i]
  infile   <- fn
  base     <- sub("\\..*$", "", fn)
  out_file <- file.path(Sys.getenv("OUT_DIR"), paste0(base, ".tsv"))

  # Skip if output already exists or source file is missing
  if (file.exists(out_file)) { cat(sprintf("[R] Skipping %s (output already exists)\n", fn)); next }
  if (!file.exists(infile))    { cat(sprintf("[R] Skipping %s (source file not found)\n", fn)); next }

  # Build column name mapping from trait metadata dictionary
  cat(sprintf("[R] Processing %s ...\n", fn))
  dt      <- fread(infile)
  rscol   <- info$rsID[i];    EAcol    <- info$EA[i]
  nonEAcol <- info$nonEA[i];  EAFcol   <- info$EAF[i]
  betacol <- info$beta[i];    secol    <- info$se[i]
  pcol    <- info$pvalue[i]

  # Direct rsID matching
  matched <- dt[get(rscol) %chin% coords,
                .(rsID    = get(rscol),
                  proxyID = NA_character_,
                  EA      = get(EAcol),
                  nonEA   = get(nonEAcol),
                  EAF     = get(EAFcol),
                  beta    = get(betacol),
                  se      = get(secol),
                  pvalue  = get(pcol),
                  LDlink  = NA_character_,
                  trueEA  = toupper(get(EAcol)))]

  # For unmatched rsIDs, search for the best LD proxy
  missing <- setdiff(coords, matched$rsID)
  for (rs in missing) {
    cat(sprintf("[R] %s: %s not matched, querying proxy...\n", base, rs))
    proxy <- tryCatch(
      system2("bash", args = c("./get_mostLD_rsID.sh", fn, rs), stdout = TRUE),
      error = function(e) character(0))
    if (length(proxy) != 1 || !startsWith(proxy, "rs")) {
      cat(sprintf("[R] %s: proxy not found (%s), skipping\n", base, paste(proxy, collapse = " ")))
      next
    }
    cat(sprintf("[R] %s: proxy found: %s\n", base, proxy))
    row_p <- dt[get(rscol) == proxy]
    if (!nrow(row_p)) next

    EA_   <- row_p[[EAcol]]; nonEA_ <- row_p[[nonEAcol]]
    EAF_  <- row_p[[EAFcol]]; beta_  <- row_p[[betacol]]
    se_   <- row_p[[secol]];  pval_  <- row_p[[pcol]]

    # Query LD relationship between proxy and target SNP for trueEA determination
    cat(sprintf("[R] %s: calling get_LDpair.sh %s %s\n", base, proxy, rs))
    ld <- tryCatch(
      system2("bash", args = c("get_LDpair.sh", proxy, rs), stdout = TRUE),
      error = function(e) "")
    cat(sprintf("[R] %s: LDlink result: %s\n", base, paste(ld, collapse = ";")))

    # Parse trueEA: determine which allele of the proxy corresponds to the target rsID EA
    te    <- ""
    ea_up <- toupper(EA_)
    for (m in strsplit(ld, ";")[[1]]) {
      pr <- strsplit(trimws(m), "->")[[1]]
      if (length(pr) == 2 && toupper(pr[1]) == ea_up) {
        te <- pr[2]
        cat(sprintf("[R] %s: trueEA determined as %s\n", base, te))
        break
      }
    }

    matched <- rbind(matched,
                     data.table(rsID    = rs,
                                proxyID = proxy,
                                EA      = EA_,   nonEA  = nonEA_,
                                EAF     = EAF_,  beta   = beta_,
                                se      = se_,   pvalue = pval_,
                                LDlink  = ld,    trueEA = te),
                     fill = TRUE)
  }

  # Write final harmonised output
  setcolorder(matched, c("rsID", "proxyID", "EA", "nonEA", "EAF",
                          "beta", "se", "pvalue", "LDlink", "trueEA"))
  fwrite(matched, out_file, sep = "\t")
  cat(sprintf("[R] Written %s (%d rows)\n", out_file, nrow(matched)))
}
EOF

echo "[INFO] All done. Results saved to: $OUT_DIR" >&2
