# =============================================================================
# Author: Haodong Tian
# Description: External pleiotropy analysis with PLATLAS
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================


library(httr2)
library(jsonlite)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(data.table)
library(readxl)
library(UpSetR)   # UpSet figure
library(ggtext)
library(patchwork) # for ggplot


### Part 1: Prepare SNP identifiers in CHR:POS:REF:ALT format (hg38) for PLATLAS query ======

# table5 contains the 9 lead MRI GWAS features: locus ID, CHR, BP, rsID, nearest gene, PheWAS results
table5 <- readRDS('/path/to/project/GWAS_res/PheWAS/the9leadMRI.rds')
table5 <- table5[order(table5$LOCUS_ID), ]

# Truncate long PheWAS trait name strings for display (retain first 30 trait names)
Nkeep <- 30
table5[, PheWASnames := sapply(strsplit(PheWASnames, "\\s*\\|\\s*"), function(v) {
  if (length(v) <= Nkeep) paste(v, collapse = " | ")
  else paste0(paste(v[1:Nkeep], collapse = " | "), " | ... (", length(v), " traits)")
})]

rsids <- unique(table5$MARKER)
length(rsids)  # 77 unique rsIDs


### Batch query Ensembl REST API to retrieve hg38 coordinates for each rsID
get_rsids_batch <- function(rsids, max_try = 5, sleep_sec = 2) {

  for (i in 1:max_try) {

    out <- try({
      x <- request("https://rest.ensembl.org/variation/human") |>
        req_headers(
          `Content-Type` = "application/json",
          Accept         = "application/json"
        ) |>
        req_body_json(list(ids = as.list(rsids))) |>
        req_method("POST") |>
        req_perform() |>
        resp_body_string() |>
        fromJSON(flatten = TRUE)

      res_list <- lapply(names(x), function(my_rsid) {
        one <- x[[my_rsid]]

        if (is.null(one$mappings) || nrow(as.data.frame(one$mappings)) == 0) {
          return(data.frame(
            rsid = my_rsid, CHR = NA, POS = NA, REF = NA, ALT = NA,
            allele_string = NA, variant_hg38 = NA, stringsAsFactors = FALSE
          ))
        }

        as.data.frame(one$mappings, stringsAsFactors = FALSE) |>
          transmute(
            rsid          = my_rsid,
            CHR           = seq_region_name,        # chromosome
            POS           = start,                   # position (hg38)
            allele_string = allele_string,
            REF           = sapply(strsplit(allele_string, "/"), `[`, 1),
            ALT           = sapply(
              strsplit(allele_string, "/"),
              function(z) if (length(z) <= 1) NA_character_ else paste(z[-1], collapse = ",")
            ),
            variant_hg38  = paste(CHR, POS, REF, ALT, sep = ":")
          )
      })
      bind_rows(res_list)
    }, silent = TRUE)

    if (!inherits(out, "try-error")) return(out)
    message("Batch query attempt ", i, "/", max_try, " failed; retrying...")
    Sys.sleep(sleep_sec * i)  # incremental back-off
  }

  stop("Batch query failed after multiple retries; try splitting rsIDs into smaller batches.")
}

res       <- get_rsids_batch(rsids)
dim(res)        # 87 7
res_clean <- res[res$CHR %in% as.character(1:22), ]
dim(res_clean)  # 77 7

# Note on REF/ALT vs. effect/non-effect allele:
# REF = reference allele in the genome build; ALT = alternate allele(s).
# These are independent of GWAS effect/non-effect allele assignment.
# For a bi-allelic SNP where GWAS EA = T and nonEA = C, the correct hg38 notation is CHR:POS:C:T
# (i.e., REF:ALT, where REF corresponds to the non-effect allele and ALT to the effect allele).
# Multi-allelic SNPs (e.g., REF=C, ALT=T/G) are treated as separate bi-allelic variants in GWAS QC.


### Retrieve EA and nonEA from QC-filtered liver MRI GWAS (MAF > 0.01) ----------------------

length(rsids)  # 77 unique rsIDs

gwas <- fread("/path/to/server/project/GWAS_regenie/MR_GWAS_new/firstorder_10Percentile_inp.regenie")  # ~2 min
dat  <- gwas[SNP %in% rsids, .(SNP, CHR, BP, nonEA = ALLELE0, EA = ALLELE1)]
dim(dat)  # 77 5 — QC-filtered GWAS has no multi-allelic entries for these rsIDs

# Merge hg19 allele info with hg38 coordinates
setDT(dat); setDT(res_clean)
setnames(dat,       c("CHR", "BP"),  c("CHR_hg19", "BP_hg19"))
setnames(res_clean, c("CHR", "POS"), c("CHR_hg38", "POS_hg38"))
dat_merge <- merge(dat, res_clean, by.x = "SNP", by.y = "rsid", all.x = TRUE)
dim(dat_merge)

# Construct final CHR:POS:REF:ALT string (hg38), matching REF/ALT to GWAS EA/nonEA
dat_final <- dat_merge[, {
  alt_vec   <- unlist(strsplit(ALT, ",", fixed = TRUE))
  hit       <- alt_vec[sapply(alt_vec, function(a) setequal(c(EA, nonEA), c(REF, a)))]
  .(CHR_hg19, BP_hg19, nonEA, EA, CHR_hg38, POS_hg38, allele_string, REF, ALT,
    ALT_final = if (length(hit) == 1) hit else NA_character_)
}, by = SNP]

dat_final[, variant_hg38_final := paste(CHR_hg38, POS_hg38, REF, ALT_final, sep = ":")]

dim(dat_final)
fwrite(dat_final, file = "/path/to/project/Pleiotropy/dat_final_hg38.tsv", sep = "\t")


### Part 2: Read PLATLAS results for the 77 variants (EUR GWAS meta-analysis) ================
# PLATLAS output file naming convention: Phe_XXX_Y.EUR.gwama.sumstats.txt.gz
# where XXX = phenotype code, Y = subtype, EUR = ancestry, GWAMA = meta-analysis method.
# ALL-ancestry files use MR-MEGA (meta-regression with PCs as effect modifiers).


### Figure 1: UpSet plot — number of variants per GenomicSEM factor category (F1-F5) ---------

# Using F-category (GenomicSEM-derived) rather than MRI-defined category because:
# 1. Genetically derived factors (F1-F5) better capture shared genetic architecture.
# 2. F-category is based on SNP-heritability |z| > 7.5, limiting to robustly heritable traits.

Fcate_info <- data.frame(
  LOCUS  = table5$LOCUS_ID,
  Fcate  = table5$CATEGORY,
  rsID   = table5$MARKER
) %>% distinct()
dim(Fcate_info)  # 92 3

# rs529565 is absent from the PLATLAS database; exclude it
Fcate_info <- Fcate_info[Fcate_info$rsID != 'rs529565', ]
dim(Fcate_info)  # 91 3

length(unique(Fcate_info$LOCUS))  # 51 unique loci (liver fat excluded from this analysis)
length(unique(Fcate_info$rsID))   # 76 unique rsIDs

# For each locus, select a single representative rsID (the one appearing most frequently across traits)
# so that internal MRI pleiotropy information from all rsIDs within a locus is consolidated.
rep_rsid <- table5[, .N, by = .(LOCUS_ID, MARKER)][
  order(LOCUS_ID, -N, MARKER)][, .SD[1], by = LOCUS_ID][, .(LOCUS = LOCUS_ID, rep_rsID = MARKER)]

Fcate_info2          <- merge(Fcate_info, rep_rsid, by = "LOCUS", all.x = TRUE)
Fcate_info2$rsID     <- Fcate_info2$rep_rsID
Fcate_info2$rep_rsID <- NULL
length(unique(Fcate_info2$rsID))  # 51 unique rsIDs (one per locus)
dim(Fcate_info2)                  # 91 3
Fcate_info2 <- unique(Fcate_info2)
dim(Fcate_info2)                  # 70 3

# Pivot to binary matrix for UpSet plot
upset_data <- Fcate_info2 %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Fcate, values_from = value, values_fill = 0) %>%
  as.data.frame()
dim(upset_data)  # 51 6

# Draw UpSet plot grouped by F-category
fcate_sets <- unique(Fcate_info2$Fcate)
upset(upset_data, sets = fcate_sets, order.by = "freq", point.size = 3, line.size = 1,
      mainbar.y.label = "Variant count", sets.x.label = "Total variants per category",
      text.scale = 1.5)

# Save as PDF (EPS cannot render UpSet shading correctly)
pdf("/path/to/project/UpSet.pdf", width = 8, height = 6, onefile = FALSE)
upset(upset_data, sets = fcate_sets, order.by = "freq", point.size = 3, line.size = 1,
      mainbar.y.label = "Variant count", sets.x.label = "Total variants per category",
      text.scale = 1.5)
dev.off()


### Figure 2: Pleiotropy bubble heatmap -------------------------------------------------------

# Analysis strategy:
# 0. Use the 51 representative rsIDs (one per locus), inheriting internal MRI pleiotropy information.
# 1. Restrict to EUR GWAS meta-analysis results.
# 2. Significance threshold: p < 0.05 / (J * K), where J = 76 variants, K = 1119 traits (Bonferroni).

Phe_dat <- read.table(
  "/path/to/project/Pleiotropy/dat_EUR_gwama.txt",
  header = TRUE, sep = "\t"
)
length(unique(Phe_dat$SUMSTATS))  # 1119 traits
dim(Phe_dat)                       # ~84103 rows ≈ 1119 traits × 76 variants (rs529565 missing)

# Reference panel for phenotype code definitions and disease categories
phe_info <- read_excel(
  "/path/to/project/Pleiotropy/phe_code_info.xlsx",
  sheet = 1
)
dim(phe_info)         # 1912 3
table(phe_info$Category)  # 17 non-null categories

# Extract trait identifier from SUMSTATS filename and merge with category labels
Phe_dat <- Phe_dat %>% mutate(Trait = str_extract(SUMSTATS, "Phe_[^.]+"))
Phe_dat <- Phe_dat %>% left_join(phe_info %>% dplyr::select(Trait, Category), by = "Trait")
length(unique(Phe_dat$Trait))     # 1119 traits
length(unique(Phe_dat$Category))  # 18 categories

# Build supplementary table: trait code -> human-readable description + disease category
Supp_phe_info       <- unique(data.frame(Trait = Phe_dat$Trait, Category = Phe_dat$Category))
dim(Supp_phe_info)  # 1119 2
Supp_phe_info$Trait <- phe_info$Description[match(Supp_phe_info$Trait, phe_info$Trait)]
saveRDS(Supp_phe_info,
        file = "/path/to/project/Pleiotropy/Supp_phe_info.rds")


# Construct the final data matrix:
# rsID | CHR:BP:REF:ALT | Fcate info | Category1 (sig/total) | Category2 | ...
upset_data
dim(upset_data)  # 51 6

final_data <- upset_data[, -1] %>%
  pivot_longer(-rsID, names_to = "Fcate") %>%
  filter(value == 1) %>%
  group_by(rsID) %>%
  summarise(Fcate = paste(Fcate, collapse = ","))
dim(final_data)  # 51 2
final_data$ID <- dat_final$variant_hg38_final[match(final_data$rsID, dat_final$SNP)]

# Bonferroni threshold: 0.05 / (1119 traits * 51 loci)
p_thresh <- 0.05 / (1119 * 51)

# Summarize significant associations per variant × disease category
phewas_summary <- Phe_dat %>%
  group_by(ID, Category) %>%
  summarise(sig = sum(P < p_thresh, na.rm = TRUE), total = n(), .groups = "drop") %>%
  mutate(ratio = paste0(sig, "/", total)) %>%
  dplyr::select(ID, Category, ratio) %>%
  pivot_wider(names_from = Category, values_from = ratio, values_fill = "0/0")

result <- final_data %>%
  left_join(phewas_summary, by = "ID") %>%
  dplyr::rename(`CHR:BP:REF:ALT` = ID, `Fcate info` = Fcate)

head(result)


### Annotate variants with nearest gene name (±0.1 Mb, Ensembl BioMart) ---------------------

get_closest_gene <- function(rsID_vector) {

  ## Step 1: retrieve GRCh38 coordinates
  m_snp <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp", mirror = "useast")
  v <- as.data.table(getBM(
    attributes = c("refsnp_id", "chr_name", "chrom_start"),
    filters    = "snp_filter", values = rsID_vector, mart = m_snp
  ))
  setnames(v, c("refsnp_id", "chr_name", "chrom_start"), c("rsID", "chr", "pos"))
  v <- v[chr %in% c(as.character(1:22), "X", "Y"), unique(.SD), by = rsID]

  ## Step 2: query genes within ±0.1 Mb windows
  win <- 1e5
  v[, `:=`(start = pmax(1L, pos - win), end = pos + win)]
  v[, region := paste0(chr, ":", start, ":", end)]
  m_gene <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  g <- as.data.table(getBM(
    attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
    filters    = "chromosomal_region", values = v$region, mart = m_gene
  ))
  setnames(g, c("chromosome_name", "start_position", "end_position", "external_gene_name"),
           c("chr", "gstart", "gend", "gene"))

  ## Step 3: find nearest gene for each rsID
  res <- rbindlist(lapply(seq_len(nrow(v)), function(i) {
    chr <- v$chr[i]; pos <- v$pos[i]; rs <- v$rsID[i]
    gg  <- as.data.frame(g)
    gg  <- gg[gg$chr == chr, ]
    gg  <- gg[(abs(gg$gstart - pos) < win) | (abs(gg$gend - pos) < win), ]
    gg  <- gg[!is.na(gg$gene) & gg$gene != "", ]  # retain genes with HGNC-style names only
    if (!nrow(gg)) return(data.table(rsID = rs, nearest_gene = NA_character_, distance_bp = NA_real_))
    d        <- ifelse(pos < gg$gstart, gg$gstart - pos, ifelse(pos > gg$gend, pos - gg$gend, 0))
    min_dist <- min(d)

    if (min_dist == 0) {
      # rsID falls within one or more gene bodies; concatenate all overlapping gene names
      target_genes <- unique(gg$gene[d == 0])
      final_gene   <- paste(target_genes, collapse = "/")
    } else {
      # rsID is intergenic; report the single nearest gene (or concatenate ties)
      target_genes <- unique(gg$gene[d == min_dist])
      final_gene   <- paste(target_genes, collapse = "/")
    }
    data.table(rsID = rs, nearest_gene = final_gene, distance_bp = min_dist)
  }))
  return(res)
}

gene_res       <- get_closest_gene(result$rsID)
result$Gene    <- gene_res$nearest_gene[match(result$rsID, gene_res$rsID)]


### Bubble plot: PLATLAS pleiotropy summary (variant × disease category) ---------------------

# Build plotting data (include all variants; zero-significant cells are rendered as absent points)
plot_dat <- result %>%
  pivot_longer(
    -c(rsID, `Fcate info`, `CHR:BP:REF:ALT`, Gene),
    names_to  = "Category",
    values_to = "ratio"
  ) %>%
  separate(ratio, into = c("sig", "total"), sep = "/", convert = TRUE) %>%
  mutate(
    rate       = sig / total,
    rate_group = cut(rate, breaks = c(0, 0.02, 0.05, 0.10, 1),
                     labels = c("<2%", "2-5%", "5-10%", ">10%")),
    label      = paste0(rsID, " (*", Gene, "*)")
  )

# Order labels by Fcate group for a structured display
label_levels <- plot_dat %>%
  dplyr::select(label, `Fcate info`) %>% distinct() %>%
  arrange(`Fcate info`) %>% pull(label)
plot_dat <- plot_dat %>% mutate(label = factor(label, levels = label_levels))

# UpSet-style Fcate annotation strip (lower panel)
fcate_dat <- plot_dat %>%
  dplyr::select(label, `Fcate info`) %>% distinct() %>%
  separate_rows(`Fcate info`, sep = ",") %>%
  mutate(`Fcate info` = str_trim(`Fcate info`), label = factor(label, levels = label_levels))

# Total significant trait counts per variant (for top histogram)
hist_dat <- plot_dat %>%
  filter(sig > 0) %>%
  group_by(label) %>%
  summarise(total_sig = sum(sig), .groups = "drop") %>%
  mutate(label = factor(label, levels = label_levels))

# Top panel: total significant association count per variant
p_hist <- ggplot(hist_dat, aes(x = label, y = total_sig)) +
  geom_col(fill = "#555555", width = 0.7) +
  scale_x_discrete(limits = label_levels) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title  = element_text(hjust = 0.5, size = 10)) +
  labs(x = NULL, y = NULL, title = "Total sig count")

# Middle panel: bubble plot (size = count, color = rate group)
p_main <- ggplot(plot_dat, aes(x = label, y = Category)) +
  geom_point(data = plot_dat %>% filter(sig > 0),
             aes(size = sig, color = rate_group), alpha = 1.0) +
  scale_size_continuous(range = c(2, 10), name = "Sig count") +
  scale_color_manual(
    values = c("<2%" = "#d0e8f1", "2-5%" = "#4a9eca", "5-10%" = "#e07b39", ">10%" = "#c0392b"),
    name   = "Rate"
  ) +
  scale_x_discrete(limits = label_levels) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, y = NULL)

# Bottom panel: UpSet-style Fcate membership strip
p_fcate <- ggplot(fcate_dat, aes(x = label, y = `Fcate info`)) +
  geom_line(
    data = fcate_dat %>% group_by(label) %>% filter(n() > 1),
    aes(group = label), color = "black", linewidth = 0.8
  ) +
  geom_point(size = 3, color = "black") +
  scale_x_discrete(limits = label_levels) +
  theme_bw() +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL)

# Combine panels (suggested export size: 1150 × 700 px)
p_hist / p_main / p_fcate + plot_layout(heights = c(2, 10, 2))


### Supplementary table: per-variant PLATLAS pleiotropy summary (51 locus-representative variants)
Supp_table <- result
dim(Supp_table)
Supp_table <- Supp_table %>%
  mutate(total_pleiotropy = rowSums(
    across(-(c(rsID, `Fcate info`, `CHR:BP:REF:ALT`, Gene)),
           ~ as.numeric(sub("/.*", "", .))),
    na.rm = TRUE
  ))
saveRDS(Supp_table, file = "/path/to/project/Pleiotropy/Supp_table.rds")

# Variants with no pleiotropy in PLATLAS
sum(Supp_table$total_pleiotropy == 0)  # 21

# Cross-check: compare PLATLAS non-pleiotropic variants with GWAS Catalog PheWAS results
PLATLAS_noplei_rsID <- Supp_table$rsID[Supp_table$total_pleiotropy == 0]
table5$num_of_PheWAS[match(PLATLAS_noplei_rsID, table5$MARKER)]
