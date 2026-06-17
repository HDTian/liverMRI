# =============================================================================
# Author: Haodong Tian
# Description: LDSC-based genetic correlation analysis between NMR metabolomics
#              traits and liver MRI radiomics features.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

### Genetic correlation between UKB NMR metabolomics and liver MRI features (LDSC)

library(data.table); library(ggplot2); library(readxl)



### Load complete rg results
rg_long <- fread('/path/to/project/Metabolites/all_rg_results.tsv')
dim(rg_long)  # 14632 x 15
length(unique(rg_long$MRI)); length(unique(rg_long$NMR))  # 59 MRI traits, 248 NMR metabolites
rg_long <- rg_long[, 1:9]
# Columns: MRI, NMR, rg, se, z, p, h2_obs, h2_obs_se, h2_int
# h2_obs: NMR SNP-heritability; h2_int: NMR univariate LDSC intercept

rg_long[, fdr := p.adjust(p, method = "fdr")]  # global FDR across all MRI x NMR pairs



### Visualization: scatter plot summarizing each NMR across all MRI features
# x-axis: number of significantly correlated MRI traits (FDR<0.05)
# y-axis: median absolute genetic correlation

library(data.table); library(ggplot2); library(readxl); library(ggrepel)

# NMR category annotations from UKB-NMR 2025 Nature Genetics paper
nmr_annot <- as.data.table(read_excel('/path/to/project/Metabolites/NMR_supp_table.xlsx', sheet = "ST1", skip = 2))
nmr_annot <- nmr_annot[, .(NMR = Biomarker, Group)]

# Summarize each NMR across all MRI features
nmr_sum <- rg_long[, .(num_of_sig_MRI = sum(fdr < 0.05, na.rm = TRUE), median_abs_rg = median(abs(rg), na.rm = TRUE)), by = NMR]
nmr_sum <- merge(nmr_sum, nmr_annot, by = "NMR", all.x = TRUE)

# Collapse into 4 broad categories for main-text figure
nmr_sum[, Group2 := fcase(
  Group %in% c("Lipoprotein subclasses","Lipoprotein particle concentrations","Lipoprotein particle sizes","Relative lipoprotein lipid concentrations","Cholesterol","Cholesteryl esters","Free cholesterol","Triglycerides","Phospholipids","Total lipids","Other lipids","Apolipoproteins"), "Lipoprotein / lipid measures",
  Group %in% c("Fatty acids"), "Fatty acids",
  Group %in% c("Amino acids","Glycolysis related metabolites","Ketone bodies"), "Small-molecule metabolites",
  Group %in% c("Inflammation","Fluid balance"), "Other systemic markers",
  default = "Other"
)]

# Top 10 NMR metabolites by number of significantly correlated MRI traits
top_nmr   <- rg_long[fdr < 0.05, .N, by = NMR][order(-N, NMR)][1:10]
label_dt  <- nmr_sum[NMR %in% top_nmr$NMR]

# Color palette by metabolite group
grp_cols4 <- c(
  "Lipoprotein / lipid measures" = "black",
  "Fatty acids"                  = "#4E7F5A",
  "Small-molecule metabolites"   = "#C97B00",
  "Other systemic markers"       = "#7A5C99",
  "Other"                        = "#9A9A9A"
)

# Final plot
p <- ggplot(nmr_sum, aes(x = num_of_sig_MRI, y = median_abs_rg, color = Group2)) +
  geom_point(size = 2, alpha = 1) +
  geom_text_repel(data = label_dt, aes(label = NMR), size = 3, nudge_x = 4,
                  direction = "y", hjust = 0, box.padding = 0.3, point.padding = 0.2,
                  segment.color = "grey60", min.segment.length = 0, show.legend = FALSE) +
  scale_color_manual(values = grp_cols4, breaks = c("Lipoprotein / lipid measures", "Fatty acids", "Small-molecule metabolites", "Other systemic markers", "Other")) +
  labs(x = "Number of genetically correlated MRI traits (FDR<0.05)", y = "Median absolute genetic correlation", color = NULL) +
  coord_cartesian(xlim = c(0, 60)) +
  theme_classic() +
  theme(axis.line = element_line(color = "black"),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.direction = "vertical",
        legend.background = element_rect(color = "black", fill = "white"),
        legend.box.background = element_rect(color = "black", fill = "white"))
p  # export at 560x530 px

saveRDS(nmr_sum, file = "/path/to/project/Metabolites/NMR_rg_results.rds")
