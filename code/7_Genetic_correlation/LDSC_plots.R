# =============================================================================
# Author: Haodong Tian
# Description: Visualization of LDSC genetic correlation results — generates
#              heatmaps and forest plots of rg estimates across traits.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

### LDSC plots

# The heatmap visualization is also implemented in LDSC.R.
# This script focuses on bar/forest plots of rg^2 for individual traits.

# Q: Should we use unadjusted or BMI-adjusted LDSC rg for visualization?
# A: Use BMI-adjusted LDSC rg. It gives a cleaner correlation structure,
#    since most traits of interest are strongly influenced by BMI.

# Q: If BMI-adjusted MRI is used for LDSC, why not use it for GenomicSEM and MR?
# A: Adjusting for heritable covariates can induce severe collider bias.
#    Although post-hoc MR methods exist to address this, we use unadjusted MRI
#    features and account for BMI effects at the MR stage.

# Q: The 59 retained liver MRI features were selected using unadjusted MRI rg.
#    Why not use BMI-adjusted rg for feature selection?
# A: Because the main MR analyses use unadjusted MRI, feature selection is also
#    based on the unadjusted results for consistency.




### Collect LDSC results  -----------------------------------------------------------------------------------
### --------------------------------------------------------------------------------------------------------

### Utility functions

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


trait_names <- sub("\\.(tsv|txt|csv)(\\.gz)?$", "", trait_info_table$`file name`, ignore.case = TRUE)

## Load non-BMI-adjusted LDSC results ------------------
logs_dir <- file.path(pathname, "LDSC/rg_all")

# Initialise result matrices
rg_mat <- matrix(NA_real_, length(liver_feature_names), length(trait_names),
                 dimnames = list(liver_feature_names, trait_names))
se_mat <- rg_mat

# List all log files
logs <- list.files(logs_dir, pattern = "_ldsc\\.log$", full.names = TRUE)
length(logs) == length( trait_names ) * length( liver_feature_names )  # must be TRUE

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

rg_unadjusted <-rg_mat
se_unadjusted <- se_mat




## Load BMI-adjusted LDSC results (instance 0) --------------------------
logs_dir <- file.path(pathname, "LDSC/rg_all_BMI_instance0")

# Initialise result matrices
rg_mat <- matrix(NA_real_, length(liver_feature_names), length(trait_names),
                 dimnames = list(liver_feature_names, trait_names))
se_mat <- rg_mat

# List all log files
logs <- list.files(logs_dir, pattern = "_ldsc\\.log$", full.names = TRUE)
length(logs) == length( trait_names ) * length( liver_feature_names )  # must be TRUE

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

rg_BMIadjust <- rg_mat
se_BMIadjust <- se_mat



### Scatter plots: unadjusted vs BMI-adjusted rg, per trait -------------------------
### ---------------------------------------------------------------------------------
trait_names <- colnames(rg_unadjusted)

par(mfrow = c(4,4), mar = c(3.5, 3.5, 2.2, 1.0), oma = c(1,0,3,0))

for (kk in seq_along(trait_names)) {
  current_trait_name <- trait_names[kk]
  x <- rg_unadjusted[, kk]
  y <- rg_BMIadjust[, kk]
  ok <- is.finite(x) & is.finite(y)

  plot(x[ok], y[ok],
       main = current_trait_name,
       xlab = " ", ylab = " ",
       pch = 16, cex = 0.8)
  abline(0, 1, col = "red", lwd = 2)
  abline(h = 0, v = 0, lty = 2, col = "grey60")
}

mtext("rg: unadjusted (x-axis) vs BMI-adjusted (y-axis)", outer = TRUE, cex = 1.0)
par(mfrow = c(1,1))




### Bar/forest plots per trait  -------------------------------------------------------------------------------
### -----------------------------------------------------------------------------------------------------------

dim(rg_BMIadjust ); dim(se_BMIadjust)  # 59 13
dim(rg_unadjusted); dim(se_unadjusted) # 19 13
colnames( rg_BMIadjust ) == colnames( rg_unadjusted )
colnames( rg_BMIadjust ) == colnames( se_BMIadjust )
colnames( rg_BMIadjust ) == colnames( se_unadjusted )
trait_names <- colnames(rg_unadjusted)

# For each trait: show the top-K MRI features most correlated in the BMI-adjusted results,
# with bars colored by direction and 95% CI overlaid for both unadjusted and adjusted estimates.

plot_rg_bar_single <- function(tr,
                               main_used =NA,
                               topK = 5,
                               xlim = c(0,1)) {
  op <- par(mar = c(5, 25, 4, 2))
  on.exit(par(op))


  rA  <- rg_BMIadjust[, tr];  seA <- se_BMIadjust[, tr]
  rU  <- rg_unadjusted[, tr]; seU <- se_unadjusted[, tr]

  idx <- order(-abs(rA))[1:topK]   # top-K by |adjusted rg|
  idx <- rev(idx)                  # display top at top of plot

  # row1 = unadjusted (light), row2 = BMI-adjusted (dark)
  mat <- rbind(abs(rU[idx]), abs(rA[idx]))

  # Direction determined from adjusted rg
  dir <- ifelse(rA[idx] > 0, "pos", "neg")

  col_U <- ifelse(dir == "pos", "pink",      "lightblue")  # unadjusted (light shade)
  col_A <- ifelse(dir == "pos", "red3",      "blue3")      # adjusted (dark shade)
  cols  <- as.vector(rbind(col_U, col_A))

  if(is.na(main_used)){main_used=tr  }
  bp <- barplot(mat, beside = TRUE, horiz = TRUE, col = cols, names.arg = rownames(rg_BMIadjust)[idx],
                las = 1, xlim = xlim, main = main_used, xlab = " ", border = "black", lwd = 0.6)

  # 95% CI: unadjusted (row 1)
  arrows(sign(rU[idx])*(rU[idx] - 1.96 * seU[idx]), bp[1,],
         sign(rU[idx])*(rU[idx] + 1.96 * seU[idx]), bp[1,], angle = 90, code = 3, length = 0.04)
  # 95% CI: adjusted (row 2)
  arrows(sign(rU[idx])*(rA[idx] - 1.96 * seA[idx]), bp[2,],
         sign(rU[idx])*(rA[idx] + 1.96 * seA[idx]), bp[2,], angle = 90, code = 3, length = 0.04)

  abline(v = 0, lty = 2, col = "grey60")
  box()

  invisible(idx)
}

trait_names


plot_rg_bar_single('LDL_C')
plot_rg_bar_single('CAD')
plot_rg_bar_single('Lp_a')

plot_rg_bar_single('TG')
plot_rg_bar_single('TG_HDL_C')
plot_rg_bar_single('liver_fat')

plot_rg_bar_single('T2D')


outdir <- "/path/to/project/LDSC/rg_barplots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


pdf(file.path(outdir, paste0('liver_cT1', "_rg_barplot.pdf")),width = 7.2 , height = 4.0)
plot_rg_bar_single('liver_iron',main_used='liver_cT1')
dev.off()

pdf(file.path(outdir, paste0('TG:HDL_C', "_rg_barplot.pdf")),width = 7.2 , height = 4.0)
plot_rg_bar_single('TG_HDL_C',main_used='TG:HDL_C')
dev.off()



### Save bar plots for all traits
for(tr in trait_names){

  pdf(file.path(outdir, paste0(tr, "_rg_barplot.pdf")),
      width = 7.2 , height = 4.0)

  plot_rg_bar_single(tr)

  dev.off()
}

# Suggested 3x3 panel layout for figures:
# liver_fat    TG     T2D
# liver_iron   Cirrhosis  liver_cancer
# LDL-C        Lp(a)      CAD
