# =============================================================================
# Author: Haodong Tian
# Description: Generates Fuji plots (multi-trait Manhattan plots) to visualize
#              GWAS associations across multiple liver MRI features simultaneously.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================


### Terminal usage:
# Rscript fujiplot.R [input.txt] [traitlist.txt] [output_dir]


# ---- RStudio mode: set paths manually ----
# Set script_dir to the local clone of the fujiplot GitHub project
# (the directory containing the config/ subdirectory)
script_dir      <- '/path/to/project/GWAS_res/fujiplot-master'

input_fname     <- '/path/to/project/GWAS_res/input.txt'     # [input.txt] path
traitlist_fname <- '/path/to/project/GWAS_res/traitlist.txt' # [traitlist.txt] path
output_dir      <- '/path/to/project/GWAS_res'               # desired output directory

args <- c(input_fname, traitlist_fname, output_dir)


#### Main pipeline ------------------------------------------------------------

library(stringr)

# Constants
.VERSION                      = "1.0.3"
CIRCOS_CONF                   = file.path(script_dir, "config", "circos.conf")
CIRCOS_PATH                   = "/opt/homebrew/bin/circos"  # verify with: which circos
CIRCOS_DEBUG_GROUP            = "summary"
OUTPUT_BARPLOT                = TRUE
SCATTER_BACKGROUND_COLOR_ALPHA = 0.3
LARGE_POINT_SIZE              = 24  # size of pleiotropic locus markers (increased from default 16)
SMALL_POINT_SIZE              = 8

# Intermediate data and config file paths (relative to fujiplot project directory)
COLOR_CONF              = file.path(script_dir, "config",      "color.conf")
SCATTER_BACKGROUND_CONF = file.path(script_dir, "config",      "scatter_background.conf")
HIGHLIGHT_DATA          = file.path(script_dir, "data_tracks", "highlights.txt")
SCATTER_DATA            = file.path(script_dir, "data_tracks", "scatter.txt")
STACKED_DATA            = file.path(script_dir, "data_tracks", "stacked.txt")
LABEL_DATA              = file.path(script_dir, "data_tracks", "label.txt")

################################################################################
writeLines(c(
  "*********************************************************************",
  "* Fuji plot -- a circos representation of multiple GWAS results",
  sprintf("* Version %s", .VERSION),
  "* Masahiro Kanai (mkanai@g.harvard.edu)",
  "* Harvard Medical School / RIKEN IMS / Osaka Univerisity",
  "* GNU General Public License v3",
  "*********************************************************************"
))

if (identical(args, character(0))) {
  args = file.path(script_dir, "input_example", c("input.txt", "traitlist.txt"))
}
input_fname     = normalizePath(args[1])
traitlist_fname = normalizePath(args[2])
if (length(args) > 2) {
  output_dir = normalizePath(args[3])
} else {
  output_dir = file.path(script_dir, 'output_example')
}

################################################################################
# Helper function
most_common = function(x) { tail(names(sort(table(x))), 1) }

################################################################################
# Load input data
message("Loading input files...")

df        = read.table(input_fname,     T, sep = '\t', quote = '', comment.char = '')
traitlist = read.table(traitlist_fname, T, sep = '\t', quote = '', comment.char = '', fileEncoding = 'utf-8')

n_loci = length(unique(df$LOCUS_ID))
writeLines(c(
  sprintf("* Input data: %s",          input_fname),
  sprintf("* Number of significant SNPs: %d", nrow(df)),
  sprintf("* Number of unique loci: %d",      n_loci),
  "",
  sprintf("* Trait list: %s",          traitlist_fname),
  sprintf("* Number of traits: %d (%s)", nrow(traitlist), str_c(traitlist$TRAIT, collapse = ',')),
  sprintf("* Number of categories: %d (%s)", length(unique(traitlist$CATEGORY)), str_c(unique(traitlist$CATEGORY), collapse = ',')),
  "",
  sprintf("* Output dir: %s",          output_dir)
))

if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

input_traits = traitlist$TRAIT
if (!all(df$TRAIT %in% input_traits)) {
  missing_traits = setdiff(df$TRAIT, input_traits)
  n_missing = length(missing_traits)
  stop(sprintf("TRAIT columns mismatch.\n%d trait%s in %s %s missing from %s (%s).",
               n_missing, ifelse(n_missing > 1, "s", ""), input_fname,
               ifelse(n_missing > 1, "are", "is"), traitlist_fname,
               str_c(missing_traits, collapse = ",")))
}

traitlist = traitlist %>%
  filter(TRAIT %in% df$TRAIT) %>%
  mutate(idx = 1:n(),
         category_lower = str_replace_all(str_to_lower(CATEGORY), '[^a-z0-9_]', '_')) %>%
  mutate(parameters = str_c('fill_color=', category_lower))

excluded_traits = setdiff(input_traits, traitlist$TRAIT)
writeLines(c(
  "",
  sprintf("Excluded %d traits because of no significant SNPs (%s).",
          length(excluded_traits), str_c(excluded_traits, collapse = ',')),
  ""
))

################################################################################
message("Generating configuration and data files for circos...")

# Write color configuration for each category
str_c_comma = function(x) { str_c(x, collapse = ",") }

cols <- traitlist %>%
  dplyr::select(category_lower, COLOR) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    rgb  = apply(t(grDevices::col2rgb(COLOR)), 1, str_c_comma),
    rgba = apply(
      floor(t((1 - SCATTER_BACKGROUND_COLOR_ALPHA) * 255 +
                SCATTER_BACKGROUND_COLOR_ALPHA * grDevices::col2rgb(COLOR))),
      1, str_c_comma
    )
  )

writeLines(c("<colors>",
             sprintf("\t%s = %s\n\talpha_%s = %s",
                     cols$category_lower, cols$rgb, cols$category_lower, cols$rgba),
             "</colors>"), COLOR_CONF)
message(sprintf("* Color configuration: %s", COLOR_CONF))


################################################################################
# Write scatter plot background configuration (shaded bands per trait category)
colsep = table(factor(traitlist$category_lower, levels = unique(traitlist$category_lower)))
bg = data.frame(
  category_lower = unique(traitlist$category_lower),
  y0 = nrow(traitlist) - cumsum(colsep) - 0.5,
  y1 = nrow(traitlist) - cumsum(c(0, colsep))[1:length(colsep)] - 0.5
)
writeLines(c("<backgrounds>",
             sprintf("<background>\n\tcolor = alpha_%s\n\ty0 = %.1f\n\ty1 = %.1f\n</background>",
                     bg$category_lower, bg$y0, bg$y1),
             "</backgrounds>"), SCATTER_BACKGROUND_CONF)
message(sprintf("* Scatter background configuration: %s", SCATTER_BACKGROUND_CONF))


################################################################################
# Write highlight data: pleiotropic loci (defined as loci associated with >1 TRAIT)
nsnps_per_locus = df %>% group_by(LOCUS_ID) %>% summarize(n = n())
df = df %>% mutate(
  CHR   = str_c("hs", CHR),
  nsnps = nsnps_per_locus$n[match(LOCUS_ID, nsnps_per_locus$LOCUS_ID)]
)

# Pleiotropy is defined as association with >1 TRAIT (not >1 CATEGORY)
inter_categorical = df %>%
  group_by(LOCUS_ID) %>%
  summarize(
    CHR = most_common(CHR),
    BP  = most_common(BP),
    n   = length(unique(TRAIT))
  ) %>%
  filter(n > 1)

write.table(inter_categorical[c("CHR", "BP", "BP")],
            HIGHLIGHT_DATA, sep = "\t", row.names = F, col.names = F, quote = F)
message(sprintf("* Highlights data (inter-categorical pleiotropic loci): %s", HIGHLIGHT_DATA))


################################################################################
# Write outer scatter plot data
scatter            = merge(df, traitlist, by = "TRAIT", all.x = T)
scatter$value      = nrow(traitlist) - scatter$idx
scatter$parameters = str_c(scatter$parameters,
                            str_c('z=', scatter$nsnps),
                            str_c('glyph_size=', ifelse(scatter$nsnps > 1, LARGE_POINT_SIZE, SMALL_POINT_SIZE)),
                            sep = ",")
scatter = scatter[order(scatter$nsnps, decreasing = T), ]
write.table(scatter[c("CHR", "BP", "BP", "value", "parameters")],
            SCATTER_DATA, sep = "\t", row.names = F, col.names = F, quote = F)
message(sprintf("* Scatter plot data (significant loci): %s", SCATTER_DATA))


################################################################################
# Write inner stacked scatter plot data (trait count per locus)
stacked   = list()
stacked_y = rep(0, n_loci)
names(stacked_y) = sort(unique(scatter$LOCUS_ID))

for (i in 1:nrow(traitlist)) {
  x        = subset(scatter, idx == i)
  x$value  = stacked_y[x$LOCUS_ID]
  stacked_y[x$LOCUS_ID] = stacked_y[x$LOCUS_ID] + 1
  stacked[[i]] = x
}

stacked            = do.call(rbind, stacked)
stacked$parameters = str_c(stacked$parameters, str_c('z=', stacked$nsnps), sep = ",")
stacked = stacked[order(stacked$nsnps, decreasing = T), ]
write.table(stacked[c("CHR", "BP", "BP", "value", "parameters")],
            STACKED_DATA, sep = "\t", row.names = F, col.names = F, quote = F)
message(sprintf("* Stacked bar plot data (# significant SNPs per locus): %s", STACKED_DATA))


################################################################################
# Write label data (gene name annotations)
# The label data frame determines which gene name labels appear on the plot;
# it is independent of the scatter/highlight point data.
# Each LOCUS_ID maps to one CHR, BP, and GENE name; duplicates in LOCUS_ID produce
# multiple labels per locus, which we prevent by deduplication on GENE below.

# Display gene labels for ALL loci with a lead SNP (not only pleiotropic ones)
label = df %>%
  group_by(LOCUS_ID) %>%
  summarize(
    CHR  = most_common(CHR),
    BP   = most_common(BP),
    GENE = most_common(GENE)
  )

# Deduplicate on GENE: if multiple nearby loci share the same nearest gene name,
# retain only the first occurrence to avoid redundant labels on the plot
label <- label %>% dplyr::distinct(GENE, .keep_all = TRUE)

write.table(label[c("CHR", "BP", "BP", "GENE")],
            LABEL_DATA, sep = "\t", row.names = F, col.names = F, quote = F)
message(sprintf("* Label data (gene name labels for all lead SNP loci): %s", LABEL_DATA))


################################################################################
# Call Circos to render the final plot
# The -conf flag specifies the main Circos configuration file (circos.conf),
# which controls all visual aspects of the plot.
# Working directory must be set to script_dir because config files use relative paths.
cmd = sprintf("%s -conf '%s'", CIRCOS_PATH, CIRCOS_CONF)
setwd(script_dir)
system(cmd)

# Output is written to the 'output/' subdirectory of the current working directory
# (controlled by Circos, not by circos.conf).
# To run directly from Terminal:
#   cd "/path/to/project/GWAS_res/fujiplot-master"
#   /opt/homebrew/bin/circos -conf "/path/to/project/GWAS_res/fujiplot-master/config/circos.conf"


# Note: post-processing blocks (file move, cleanup, barplot) are retained below
# as commented-out code for reference but are not needed when running interactively.

# # move the output file from circos to the specified location
# if (output_dir != file.path(script_dir, 'output')){
#   for(ext in c('png', 'svg')){
#     cmd <- sprintf("mv %s %s",
#                    shQuote(file.path(script_dir, "output", sprintf("circos.%s", ext))),
#                    shQuote(file.path(output_dir, sprintf("circos.%s", ext))))
#     system(cmd)
#   }
# }
#
# # clean-up the intermediate files
# for(f in c(
#   COLOR_CONF,
#   SCATTER_BACKGROUND_CONF,
#   HIGHLIGHT_DATA,
#   SCATTER_DATA,
#   STACKED_DATA,
#   LABEL_DATA
# )){
#   if (file.exists(f)) file.remove(f)
# }
#
# ################################################################################
# # output bar plot
# if (OUTPUT_BARPLOT) {
#   bar = df %>% group_by(TRAIT) %>%
#                summarize(total = length(MARKER),
#                          pleiotropic = length(MARKER[nsnps > 1]),
#                          inter_categorical = length(MARKER[LOCUS_ID %in% inter_categorical$LOCUS_ID])) %>%
#                mutate(single = total - pleiotropic,
#                       intra_categorical = pleiotropic - inter_categorical)
#   bar = merge(bar, traitlist, by = "TRAIT")
#   bar = bar[order(bar$idx, decreasing=T),]
#   rownames(bar) = bar$TRAIT
#
#   cairo_pdf(file.path(output_dir, "barplot.pdf"), width = 8, height = 8, family = "Helvetica")
#     barplot(t(bar[,c("inter_categorical", "intra_categorical", "single")]), ylim = c(0, 100), space = 0, col = c("black", "grey50", "white"))
#   . = dev.off()
# }
#
# writeLines(c("", "",
#              sprintf("* Final circos outputs: %s.{png,svg}.", file.path(output_dir, 'circos')),
#              "",
#      sprintf("Finished at %s.", Sys.time())
#            ))
