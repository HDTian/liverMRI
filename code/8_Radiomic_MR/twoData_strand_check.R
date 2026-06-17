# =============================================================================
# Author: Haodong Tian
# Description: Utility function to verify strand alignment between two GWAS
#              summary datasets before harmonization. Identifies palindromic
#              vs. non-palindromic SNPs and reports same-strand rates.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

## Input:  Two data frames containing rsID, EA, EAF, and optionally nonEA columns.
## Output:
##   - A count table: palindromic/non-palindromic SNP x same/different strand
##   - Number of ambiguous palindromic SNPs (0.45 < EAF < 0.55)

twoData_strand_check <- function(Data1, Data2) {

  # Auto-detect column names from common naming conventions
  find_col <- function(df, patterns) {
    cols <- colnames(df)
    for (pat in patterns) {
      hit <- grep(pat, cols, ignore.case = TRUE, value = TRUE)
      if (length(hit) > 0) return(hit[1])
    }
    return(NULL)
  }

  col_patterns <- list(
    rsID  = c("rsid", "snp", "rs_id", "snpid", "markername", "id"),
    EA    = c("EA", "effect_allele", "effectallele", "ALLELE1", "a1"),
    nonEA = c("nonEA", "non_effect", "noneffect", "other_allele", "otherallele", "a2"),
    EAF   = c("EAF", "effect_allele_freq", "freq", "frq", "maf", "a1freq")
  )

  find_cols <- function(df) lapply(col_patterns, find_col, df = df)
  c1 <- find_cols(Data1); c2 <- find_cols(Data2)

  # Require rsID, EA, EAF; nonEA is optional
  for (field in c("rsID", "EA", "EAF")) {
    if (is.null(c1[[field]])) stop(paste("Data1 missing column:", field))
    if (is.null(c2[[field]])) stop(paste("Data2 missing column:", field))
  }

  has_nonEA1 <- !is.null(c1$nonEA)
  has_nonEA2 <- !is.null(c2$nonEA)
  can_check_nonpal <- has_nonEA1 | has_nonEA2

  # Extract standardised columns; fill nonEA with NA if absent
  extract <- function(df, cols) data.frame(
    rsID  = tolower(df[[cols$rsID]]),
    EA    = toupper(df[[cols$EA]]),
    nonEA = if (!is.null(cols$nonEA)) toupper(df[[cols$nonEA]]) else NA_character_,
    EAF   = df[[cols$EAF]]
  )
  D1 <- extract(Data1, c1); D2 <- extract(Data2, c2)

  # Deduplicate on rsID, keeping first occurrence
  if (any(duplicated(D1$rsID))) { warning("Data1 has duplicate rsIDs; keeping first."); D1 <- D1[!duplicated(D1$rsID), ] }
  if (any(duplicated(D2$rsID))) { warning("Data2 has duplicate rsIDs; keeping first."); D2 <- D2[!duplicated(D2$rsID), ] }

  shared_rsID <- intersect(D1$rsID, D2$rsID)
  if (length(shared_rsID) == 0) stop("No shared rsIDs between the two datasets!")

  D1 <- D1[D1$rsID %in% shared_rsID, ]
  D2 <- D2[D2$rsID %in% shared_rsID, ]
  D2 <- D2[match(D1$rsID, D2$rsID), ]

  # Identify palindromic SNPs using whichever dataset has nonEA (prefer Data1)
  is_palindromic <- if (has_nonEA1) {
    paste(D1$EA, D1$nonEA) %in% c("A T", "T A", "C G", "G C")
  } else {
    paste(D2$EA, D2$nonEA) %in% c("A T", "T A", "C G", "G C")
  }

  # Non-palindromic SNPs: same strand if all four alleles (EA/nonEA from both datasets) use ≤2 distinct letters
  non_pal <- !is_palindromic
  if (!can_check_nonpal) {
    warning("Neither dataset has nonEA — cannot check strand for non-palindromic SNPs.")
    non_pal_same <- non_pal_diff <- NA
  } else {
    is_same <- mapply(function(...) length(unique(na.omit(c(...)))) <= 2,
                      D1$EA[non_pal], D1$nonEA[non_pal], D2$EA[non_pal], D2$nonEA[non_pal])
    non_pal_same <- sum(is_same)
    non_pal_diff <- sum(non_pal) - non_pal_same
  }

  # Palindromic SNPs: ambiguous if EAF near 0.5; for clear ones, check EAF direction agreement
  pal_ambig <- sum(is_palindromic & D1$EAF >= 0.45 & D1$EAF <= 0.55)
  pal_clear <- is_palindromic & !(D1$EAF >= 0.45 & D1$EAF <= 0.55)

  # Align EAF: ensure both datasets refer to the same allele before comparing direction
  D2_EAF_aligned <- ifelse(D2$EA == D1$EA, D2$EAF, 1 - D2$EAF)

  pal_same <- sum(pal_clear & ((D1$EAF > 0.5) == (D2_EAF_aligned > 0.5)))
  pal_diff <- sum(pal_clear) - pal_same


  result_table <- matrix(
    c(non_pal_same, non_pal_diff, pal_same, pal_diff),
    nrow = 2, dimnames = list(
      c("Same strand", "Different strand"),
      c("Non-palindromic SNP", "Palindromic SNP (EAF not close to 0.5)")
    )
  )

  cat("=== same strand check results ===\n")
  cat("Data1:", paste(names(c1), sapply(c1, function(x) if(is.null(x)) "?" else x), sep = "=", collapse = ", "), "\n")
  cat("Data2:", paste(names(c2), sapply(c2, function(x) if(is.null(x)) "?" else x), sep = "=", collapse = ", "), "\n")
  cat(sprintf("Shared rsIDs: %d\n", length(shared_rsID)))
  print(result_table)
  cat(sprintf("\n  Ambiguous palindromic SNPs (EAF in [0.45, 0.55]): %d", pal_ambig))
  cat(paste0(" (same-strand rate for remaining SNPs: ",
             sum(result_table[1,]) / sum(result_table) * 100, "%\n"))
  invisible(list(table = result_table, n_ambiguous_palindromic = pal_ambig))
}
