# =============================================================================
# Author: Haodong Tian
# Description: Generates a horizontal-layout heatmap for displaying genetic
#              correlation matrices across multiple trait groups.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

### Horizontal heatmap of LDSC genetic correlations
### (horizontal layout variant of the heatmap in LDSC.R, with features on the x-axis)


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

## Read liver feature clustering information (from GenomicSEM.R)
cluster_df <- readRDS(paste0(pathname, "/LDSC/Rg_cluster_df.rds")); names(cluster_df)[1] <- "feature"

# 1) Sort features by cluster (preserving original order within each cluster)
cf <- data.table(feature = rownames(rg_mat), idx = seq_len(nrow(rg_mat)))[cluster_df, on = "feature"]
if (anyNA(cf$cluster)) {
  mx <- suppressWarnings(max(cf$cluster, na.rm = TRUE)); if (!is.finite(mx)) mx <- 0
  cf[is.na(cluster), cluster := mx + 1L]
}
setorder(cf, cluster, idx); feature_levels <- cf$feature

## Set trait display order
trait_levels <- c('BMI','TG','HDL_C','TG_HDL_C','Lp_a','LDL_C','CAD','HbA1c','T2D','liver_fat','liver_iron','Cirrhosis','liver_cancer')

# Filter and reorder trait columns
rg_mat <- rg_mat[, colnames(rg_mat) %in% trait_levels, drop = FALSE]
se_mat <- se_mat[, colnames(se_mat) %in% trait_levels, drop = FALSE]
rg_mat <- rg_mat[, trait_levels[trait_levels %in% colnames(rg_mat)], drop = FALSE]
se_mat <- se_mat[, trait_levels[trait_levels %in% colnames(se_mat)], drop = FALSE]
length(colnames(rg_mat)) == length(trait_levels)

# 2) Wide to long; compute p-values and significance flags; set tile sizes for the heatmap
dt_rg <- melt(as.data.table(rg_mat, keep.rownames = "feature"), id = "feature", variable.name = "trait", value.name = "rg")
dt_se <- melt(as.data.table(se_mat, keep.rownames = "feature"), id = "feature", variable.name = "trait", value.name = "se")
dt <- merge(dt_rg, dt_se, by = c("feature","trait"), all = TRUE)
data.table::setDT(dt)

dt[, z := rg / se]
dt[, p := fifelse(is.finite(z), 2 * pnorm(abs(z), lower.tail = FALSE), NA_real_)]
dt[, rg_nom := fifelse(!is.na(p) & p < 0.05, rg, NA_real_)]
dt[, q := p.adjust(p, method = "fdr")]                  # global FDR across all pairs
dt[, rg_fdr := fifelse(!is.na(q) & q < 0.05, rg, NA_real_)]
alpha_star <- 0.05 / (dim(rg_mat)[1] * dim(rg_mat)[2]) # global Bonferroni threshold
dt[, star := fifelse(!is.na(p) & p < alpha_star, "*", "")]

dt[, z_abs := abs(rg / se)]
z_max <- min(abs(dt$z_abs)[!is.na(dt$rg_fdr)])           # minimum |z| among FDR-significant pairs
dt[, tile_size := pmin(z_abs / z_max, 1)]
z_min <- min(dt$tile_size[(dt$tile_size != 1) & (!is.na(dt$rg_nom))])
dt[, tile_size := 0.5 + 0.3 * (tile_size - z_min) / (1 - z_min)]  # map to [0.5, 0.8]
dt$tile_size[dt$tile_size < 0.5] <- 0

# Apply factor levels for display order (features on x-axis in horizontal layout)
dt[, feature := factor(feature, levels = feature_levels)]
dt[, trait := factor(trait, levels = trait_levels)]

# 3) Compute cluster bounding boxes (columns, since features are on x-axis) and labels
pos <- setNames(seq_along(feature_levels), feature_levels)
cb <- rbindlist(lapply(split(cf$feature, cf$cluster, drop = TRUE), function(v){
  r <- range(pos[v])
  data.table(xmin = r[1] - 0.5, xmax = r[2] + 0.5, ymin = 0.5, ymax = length(trait_levels) + 0.5)
}), fill = TRUE)

lab_df <- rbindlist(lapply(split(cf$feature, cf$cluster, drop = TRUE), function(v){
  r <- range(pos[v])
  data.table(cluster = unique(cf[feature %in% v, cluster]), x = mean(r))
}), fill = TRUE)
if (nrow(lab_df)) {
  lab_df[, `:=`(label = paste0("C", cluster), y = length(trait_levels) + 0.65)]
}

dt[, rg_show := fifelse(!is.na(p) & p <= alpha_star, rg, NA_real_)]
dt[, star_col := ifelse(star == "", NA, "black")]

# Rename traits for display
dt$trait <- as.character(dt$trait)
dt$trait[dt$trait == "liver_iron"] <- "liver_cT1"
dt$trait[dt$trait == "TG_HDL_C"] <- "TG:HDL_C"
trait_levels <- c('BMI','TG','HDL_C','TG:HDL_C','Lp_a','LDL_C','CAD','HbA1c','T2D','liver_fat','liver_cT1','Cirrhosis','liver_cancer')
dt$trait <- factor(dt$trait, levels = rev(trait_levels))

# Display labels for y-axis (trait names)
outcome_labels <- c(
  "TG"="TG", "HDL_C"="HDL-C", "TG:HDL_C"="TG:HDL-C", "Lp_a"="Lp(a)",
  "LDL_C"="LDL-C", "CAD"="CAD", "HbA1c"="HbA1c", "T2D"="T2D",
  "liver_fat"="Liver fat", "liver_iron"="Liver iron",
  "Cirrhosis"="Cirrhosis", "liver_cancer"="Liver cancer", "liver_cT1"="Liver cT1"
)

# 4) Horizontal heatmap: x = feature, y = trait
p_rg <- ggplot(dt, aes(feature, trait)) +
  geom_tile(fill = "white", color = "grey90", linewidth = 0.2) +
  # Background grid
  geom_tile(aes(fill = rg_fdr), width = 0.98, height = 0.98, color = NA) +
  # FDR-significant tiles: near-full size
  geom_tile(aes(fill = rg_nom, width = tile_size, height = tile_size), color = NA) +
  # Nominally significant tiles: sized proportionally to |z|
  geom_tile(fill = NA, color = "grey70", linewidth = 0.2) +
  # Overlay grid lines
  geom_text(aes(label = star, color = star_col), size = 3.6, vjust = 0.7, show.legend = FALSE) +
  # Star annotation for Bonferroni-significant pairs
  { if (nrow(cb)) geom_rect(data = cb, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = NA, color = "grey40", linewidth = 0.4) } +
  # Cluster bounding boxes (vertical strips for horizontal layout)
  { if (nrow(lab_df)) geom_text(data = lab_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, vjust = 0, size = 3.3, color = "grey20") } +
  # Cluster labels above the heatmap
  scale_color_identity() +
  scale_fill_gradient2(name = NULL, midpoint = 0, low = "blue", mid = "white", high = "red",
                       limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5), na.value = "white") +
  scale_x_discrete(position = "bottom") + scale_y_discrete(labels = outcome_labels) +
  coord_fixed(ratio = 1, xlim = c(0.5, length(feature_levels) + 0.5), clip = "off") +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colorbar(direction = "vertical", barheight = grid::unit(2.5, "cm"),
                               barwidth = grid::unit(0.35, "cm"), ticks = FALSE,
                               frame.colour = "black", title.position = "top",
                               title.hjust = 0.5, label.position = "right")) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0, face = "bold"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.ticks = element_blank(),
    plot.margin = margin(18, 10, 10, 80)
  )


p_rg

## Heatmap interpretation:
# white/null: nominal p-value > 0.05
# small tile: nominal p-value < 0.05
# full tile:  global FDR q < 0.05
# star (*):   global Bonferroni p < 0.05
