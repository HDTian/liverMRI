# =============================================================================
# Author: Haodong Tian
# Description: LDSC genetic correlation analysis specifically for WHR adjusted
#              for BMI and related anthropometric/metabolic traits against liver
#              MRI features.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

### LDSC analysis: comparing genetic correlations with BMI vs WHRadjBMI

## Pipeline overview (see LDSC.R for full implementation):
## Step 0: Standardise external trait GWAS files into a clean format for munging
## Step 1: Munge all traits/phenotypes from local PheWAS GWAS results into *.sumstats
## Step 2: Bivariate LDSC for liver_features x traits to estimate genetic correlations


trait_info_table <- read_excel( paste0(pathname , '/trait_info.xlsx' ) );setDT(trait_info_table)
trait_info_table <- trait_info_table[!is.na(`file name`)]

feature_paths <- paste0(pathname , '/PheWAS/' ,  trait_info_table$`file name`    )

### Collect LDSC results  -----------------------------------------------------------------------------------
### --------------------------------------------------------------------------------------------------------
trait_names <- sub("\\.(tsv|txt|csv)(\\.gz)?$", "", trait_info_table$`file name`, ignore.case = TRUE)
logs_dir <- file.path(pathname, "LDSC/rg_all")      ## LDSC with non-BMI-adjusted liver feature GWAS


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

dim(rg_mat); dim(se_mat)  # 59 14



### Scatter plot: |rg with BMI| vs |rg with WHRadjBMI| across 59 MRI features

rg_mat<-as.data.frame(rg_mat)
se_mat<-as.data.frame(se_mat)
plot( abs(rg_mat$BMI/se_mat$BMI),  abs(rg_mat$WHRadjBMI/se_mat$WHRadjBMI),xlim=c(0,15),ylim=c(0,15) ) # z-score scatterplot
abline(0,1,col='blue')


x <- abs(rg_mat$BMI);y <- abs(rg_mat$WHRadjBMI)
plot(x, y, pch = 16, xlim=c(0,0.6),ylim=c(0,0.6) , xlab = "|rg| with BMI", ylab = "|rg| with WHRadjBMI")
abline(0, 1, lty = 2, col = "grey60")
points(x, y, pch = 16)




library(ggplot2)

thres <- 0.05 / nrow(rg_mat)

x <- abs(rg_mat$BMI); y <- abs(rg_mat$WHRadjBMI)
p_bmi <- 2 * pnorm(-abs(rg_mat$BMI / se_mat$BMI))
p_whr <- 2 * pnorm(-abs(rg_mat$WHRadjBMI / se_mat$WHRadjBMI))

sig_group <- ifelse(p_bmi < thres & p_whr < thres, "Both significant",
                    ifelse(p_bmi < thres & p_whr >= thres, "BMI only",
                           ifelse(p_bmi >= thres & p_whr < thres, "WHRadjBMI only", "Neither")))
sig_group <- factor(sig_group, levels = c("BMI only", "WHRadjBMI only", "Both significant", "Neither"))

plot_dt <- data.frame(MRI = rownames(rg_mat), x = x, y = y, p_bmi = p_bmi, p_whr = p_whr, sig_group = sig_group, stringsAsFactors = FALSE)
plot_dt <- plot_dt[is.finite(plot_dt$x) & is.finite(plot_dt$y), ]

# Select labels: top 2 by BMI rg and top 2 by WHRadjBMI rg from "Both significant";
# plus top 2 from each of the single-trait-significant groups.
lab_both_x <- head(plot_dt[plot_dt$sig_group == "Both significant", ][order(-plot_dt[plot_dt$sig_group == "Both significant", ]$x), ], 2)
lab_both_y <- head(plot_dt[plot_dt$sig_group == "Both significant", ][order(-plot_dt[plot_dt$sig_group == "Both significant", ]$y), ], 2)
lab_bmi    <- head(plot_dt[plot_dt$sig_group == "BMI only", ][order(-plot_dt[plot_dt$sig_group == "BMI only", ]$x), ], 2)
lab_whr    <- head(plot_dt[plot_dt$sig_group == "WHRadjBMI only", ][order(-plot_dt[plot_dt$sig_group == "WHRadjBMI only", ]$y), ], 2)

label_dt <- unique(rbind(lab_both_x, lab_both_y, lab_bmi, lab_whr))
rownames(label_dt) <- NULL

# Place all labels directly to the right of their respective points
label_dt$x_lab <- pmin(label_dt$x + 0.015, 0.695)

sig_cols <- c("BMI only" = "#4A6FA5", "WHRadjBMI only" = "#4E7F5A", "Both significant" = "#B05A5A", "Neither" = "#BDBDBD")

p <- ggplot(plot_dt, aes(x = x, y = y, color = sig_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey80") +
  geom_point(size = 2) +
  geom_text(data = label_dt, aes(x = x_lab, y = y, label = MRI), inherit.aes = FALSE, hjust = 0, size = 3) +
  scale_color_manual(values = sig_cols, breaks = c("BMI only", "WHRadjBMI only", "Both significant", "Neither")) +
  coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 0.6), clip = "off") +
  labs(x = "Absolute genetic correlation with BMI", y = "Absolute genetic correlation with WHRadjBMI", color = NULL)+
  theme_classic() +
  theme(axis.line = element_line(color = "black"),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", fill = "white"),
        legend.box.background = element_rect(color = "black", fill = "white"),
        plot.margin = margin(5.5, 140, 5.5, 5.5))

p

saveRDS(plot_dt, file = "/path/to/project/Metabolites/rg_BMI_WHRadjBMI.rds")


table(plot_dt$sig_group)
# BMI only   WHRadjBMI only Both significant          Neither
#  22                7               19               11
