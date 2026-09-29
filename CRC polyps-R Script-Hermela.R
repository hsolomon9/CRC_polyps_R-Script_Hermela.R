# ==============================================================================
# Driver Mutation Cases and Polyp Location Analysis in Adenoma Pathways to
# Colorectal Cancer
#
# Reproducible analysis script — Appendix
# Author: Hermela Solomon
# INCISE cohort | Patrick G. Johnston Centre for Cancer Research, QUB
#
# This script reproduces every analysis and figure reported in the Methods
# and Results sections, organised into the same order the thesis presents
# them. Analyses include the six dualGSEA comparisons, their accompanying
# boxplots, and the decoupleR TF comparisons. 
#
# NOTE ON REPRODUCIBILITY: file paths below point to this project's local
# directory structure and will need updating to run on another machine.
# The CollecTRI network file and the stem cell signature file (CBC/RSC/
# proCSC/revCSC) are supervisor-provided resources 
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. SESSION SETUP
# ------------------------------------------------------------------------------

set.seed(121)  # matches the seed used for decoupleR run_ulm (Section 10)

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(gtsummary)
library(rstatix)
library(PCAtools)
library(GSVA)
library(msigdbr)
library(fgsea)
library(decoupleR)
library(ggvenn)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)

# /Users/hermelasolomon/Documents/Dissertation/dualGSEA_function.R

source(file.path(D, "dualGSEA.R"))

# Confirm the package versions actually used, for the Methods section table.

pkgs_to_check <- c("dplyr", "tidyr", "purrr", "ggplot2", "gridExtra", "ggpubr",
                   "PCAtools", "GSVA", "msigdbr", "fgsea", "decoupleR",
                   "ggvenn", "ComplexHeatmap", "circlize", "ggrepel",
                   "gtsummary", "rstatix")
print(data.frame(
  Package = pkgs_to_check,
  Version = sapply(pkgs_to_check, function(p)
    as.character(tryCatch(packageVersion(p), error = function(e) NA)))
))
cat("R version:", R.version.string, "\n")

# Project directories
D      <- "/Users/hermelasolomon/Documents/Dissertation/"
PLOTS  <- file.path(D, "Plots")
RES    <- file.path(D, "Results")
dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
dir.create(RES,   recursive = TRUE, showWarnings = FALSE)

# Shared plotting theme — used by every figure in this script.
# base_size = 15 keeps axis/body text legible at both screen and print size;
# titles and strip labels are stepped up from that base rather than fixed

theme_thesis <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(size = base_size + 3, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size,     hjust = 0.5, colour = "grey30"),
      axis.title    = element_text(size = base_size + 1, face = "bold"),
      axis.text     = element_text(size = base_size,     colour = "black"),
      legend.title  = element_text(size = base_size + 1, face = "bold"),
      legend.text   = element_text(size = base_size),
      strip.text    = element_text(size = base_size + 1, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Two-colour direction palette used across all barplots/boxplots comparing
# two groups. Matched to the palette used consistently throughout the
# original analysis script — "steelblue"/"palevioletred" for Left/Right
# colon (used identically in >15 figures across that script), and the
# #2166AC/#B2182B pairing that the APC-only vs APC+KRAS volcano plot

COL_LEFT_RIGHT <- c("Left_colon" = "steelblue",  "Right_colon" = "palevioletred")
COL_RECTUM     <- "mediumpurple"  # third colour when rectum is included (Section 9's Figure 8)
COL_MUT        <- c("APC_only" = "#2166AC", "APC_KRAS" = "#B2182B")

# Pretty display text for the group names above. 
DISPLAY_LABELS <- c(Left_colon = "Left colon", Right_colon = "Right colon",
                    APC_only = "APC-only", APC_KRAS = "APC+KRAS")


# ------------------------------------------------------------------------------
# 1. DATA LOADING
# ------------------------------------------------------------------------------

metadata <- read.table(file.path(D, "INCISE_Eman_March2026_combined_metadata.csv"),
                       sep = ",", header = TRUE)
expr_data <- readRDS(file.path(D, "INCISE_Eman_March2026_combined_expression.rds"))
mut_data  <- readRDS(file.path(D, "INCISE_Eman_March2026_mutation_binary_calls.rds"))

dim(expr_data)  # genes x samples (expect 14,993 x ~2,494 pre-matching)
dim(mut_data)   # samples x 397 mutation columns (+ Study.ID before cleanup)

sum(duplicated(rownames(expr_data)))
sum(duplicated(metadata$Study.ID))


# ------------------------------------------------------------------------------
# 2. DATA CLEANING & MATCHING
# ------------------------------------------------------------------------------

rownames(mut_data) <- mut_data$Study.ID
mut_data$Study.ID  <- NULL

common_samples <- Reduce(intersect, list(
  colnames(expr_data),
  rownames(mut_data),
  metadata$Study.ID
))
length(common_samples)  # 760

expr_clean <- expr_data[, common_samples]
mut_clean  <- mut_data[common_samples, ]
meta_clean <- metadata[metadata$Study.ID %in% common_samples, ]

dim(expr_clean)   # 14,993 x 760
dim(mut_clean)    # 760 x 397
dim(meta_clean)   # 760 x 25

saveRDS(expr_clean, file.path(D, "expr_clean.rds"))
saveRDS(mut_clean,  file.path(D, "mut_clean.rds"))
saveRDS(meta_clean, file.path(D, "meta_clean.rds"))

# Adenoma subset — all downstream analysis uses this exclusively (n=744),
# following the rationale in Methods §2 (serrated n=16 too small for
# reliable comparative statistics).
adenoma_ids  <- meta_clean$Study.ID[meta_clean$Adenoma_vs_Serrated == "Adenoma"]
expr_adenoma <- expr_clean[, adenoma_ids]
meta_adenoma <- meta_clean[meta_clean$Study.ID %in% adenoma_ids, ]
mut_adenoma  <- mut_clean[adenoma_ids, ]

dim(expr_adenoma)  # 14,993 x 744

saveRDS(expr_adenoma, file.path(D, "expr_adenoma.rds"))
saveRDS(meta_adenoma, file.path(D, "meta_adenoma.rds"))
saveRDS(mut_adenoma,  file.path(D, "mut_adenoma.rds"))

# Serrated subset — retained only for the three-way mutation-frequency
# comparison in Figure 4 / Table 1, not used in any pathway-level analysis.
serrated_ids  <- meta_clean$Study.ID[meta_clean$Adenoma_vs_Serrated == "Serrated"]
mut_serrated  <- mut_clean[serrated_ids, ]

cat("Full cohort:", nrow(meta_clean),
    "| Adenoma:", nrow(meta_adenoma),
    "| Serrated:", length(serrated_ids), "\n")


# ------------------------------------------------------------------------------
# 3. MUTATION COMBINATION LABELS
# ------------------------------------------------------------------------------
# Defined once here and reused everywhere a combo subgroup is needed

add_combo <- function(mut_df) {
  mut_df$combo <- paste0(
    ifelse(mut_df$APC  == 1, "A", ""),
    ifelse(mut_df$KRAS == 1, "K", ""),
    ifelse(mut_df$BRAF == 1, "B", "")
  )
  mut_df$combo[mut_df$combo == ""] <- "Other"
  mut_df
}

mut_clean    <- add_combo(mut_clean)
mut_adenoma  <- add_combo(mut_adenoma)
mut_serrated <- add_combo(mut_serrated)


# ------------------------------------------------------------------------------
# 4. CLINICOPATHOLOGICAL SUMMARY — FULL COHORT (Table 1 / Figure 3)
# ------------------------------------------------------------------------------

meta_clean$Location <- factor(meta_clean$Location,
                              levels = c("Left_colon", "Rectum", "Right_colon"))

table1_full <- meta_clean %>%
  select(Sex, Age, Location, Adenoma_vs_Serrated) %>%
  tbl_summary(
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 1,
    missing = "no"
  ) %>%
  bold_labels()

table1_full %>% as_tibble() %>%
  write.csv(file.path(RES, "Table1_full_cohort.csv"), row.names = FALSE)

# Figure 3B — sex, polyp type, location as bar plots; age as histogram.
# Built as four independent panels combined with ggarrange, matching the
# reporting style used throughout (panel-lettered, common theme, shared font size).

make_pct_bar <- function(var, var_label, fill_col = "steelblue") {
  df <- as.data.frame(table(var))
  colnames(df) <- c("Group", "Count")
  df$Percentage <- round(df$Count / sum(df$Count) * 100, 1)
  ggplot(df, aes(x = Group, y = Percentage)) +
    geom_col(fill = fill_col, width = 0.7) +
    geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5, size = 5) +
    labs(x = var_label, y = "Percentage (%)") +
    theme_thesis() +
    theme(legend.position = "none")
}

p_sex      <- make_pct_bar(meta_clean$Sex, "Sex")
p_type     <- make_pct_bar(meta_clean$Adenoma_vs_Serrated, "Polyp type")
p_location <- make_pct_bar(meta_clean$Location, "Polyp location")

p_age <- ggplot(meta_clean, aes(x = Age)) +
  geom_histogram(fill = "steelblue", colour = "white", bins = 20) +
  labs(x = "Age (years)", y = "Count") +
  theme_thesis()

fig3B <- ggarrange(p_sex, p_type, p_location, p_age,
                   ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(file.path(PLOTS, "Figure3_cohort_overview.png"), fig3B,
       width = 12, height = 10, dpi = 300)


# ------------------------------------------------------------------------------
# 5. ADENOMA vs SERRATED COMPARISON (Table 1 — two polyp types)
# ------------------------------------------------------------------------------

table1_by_type <- meta_clean %>%
  mutate(
    APC  = factor(mut_clean$APC,  levels = c(0, 1), labels = c("Wild-type", "Mutant")),
    KRAS = factor(mut_clean$KRAS, levels = c(0, 1), labels = c("Wild-type", "Mutant")),
    BRAF = factor(mut_clean$BRAF, levels = c(0, 1), labels = c("Wild-type", "Mutant"))
  ) %>%
  select(Sex, Age, APC, KRAS, BRAF, Adenoma_vs_Serrated) %>%
  tbl_summary(
    by = Adenoma_vs_Serrated,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 1,
    missing = "no"
  ) %>%
  add_p() %>%
  bold_labels() %>%
  bold_p()

table1_by_type %>% as_tibble() %>%
  write.csv(file.path(RES, "Table1_adenoma_vs_serrated.csv"), row.names = FALSE)


# ------------------------------------------------------------------------------
# 6. TOP 10 MUTATED GENES (Figure 4: full cohort / adenoma / serrated)
# ------------------------------------------------------------------------------

make_top10_bar <- function(mut_df, ids, title, fill_col = "steelblue",
                           drop_genes = NULL, n_top = 10) {
  m <- mut_df[ids, sapply(mut_df, is.numeric), drop = FALSE]
  df <- data.frame(Gene = colnames(m),
                   n = colSums(m == 1, na.rm = TRUE)) %>%
    mutate(Percentage = round(100 * n / length(ids), 1)) %>%
    filter(!Gene %in% drop_genes) %>%
    arrange(desc(Percentage)) %>%
    head(n_top)
  
  ggplot(df, aes(x = reorder(Gene, Percentage), y = Percentage)) +
    geom_col(fill = fill_col, width = 0.75) +
    geom_text(aes(label = paste0(Percentage, "%")), hjust = -0.15, size = 4.5, fontface = "bold") +
    coord_flip() +
    ylim(0, max(df$Percentage) * 1.18) +
    labs(title = paste0(title, " (n=", length(ids), ")"),
         x = "Gene", y = "Mutation frequency (%)") +
    theme_thesis(base_size = 13)
}

p_top10_full     <- make_top10_bar(mut_clean, rownames(mut_clean), "Full cohort", "purple4")
p_top10_adenoma  <- make_top10_bar(mut_adenoma, rownames(mut_adenoma), "Adenoma subset", "steelblue")
p_top10_serrated <- make_top10_bar(mut_serrated, rownames(mut_serrated), "Serrated subset", "darkgreen")

fig4 <- ggarrange(p_top10_full, p_top10_adenoma, p_top10_serrated,
                  ncol = 3, labels = c("A", "B", "C"))
ggsave(file.path(PLOTS, "Figure4_top10_mutated_genes.png"), fig4,
       width = 18, height = 7, dpi = 300)


# ------------------------------------------------------------------------------
# 7. TOP 10 MUTATED GENES BY LOCATION AND MUTATION SUBGROUP (Figure 5)
# ------------------------------------------------------------------------------
# Left/right colon subset established here is reused for every subsequent
# analysis in this script (PCA restricted comparisons, ssGSEA, dualGSEA,
# decoupleR) 

meta_lr <- meta_adenoma[meta_adenoma$Location %in% c("Left_colon", "Right_colon"), ]
expr_lr <- expr_adenoma[, meta_lr$Study.ID]
mut_lr  <- mut_adenoma[rownames(mut_adenoma) %in% meta_lr$Study.ID, ]

left_ids  <- meta_lr$Study.ID[meta_lr$Location == "Left_colon"]   # n=526
right_ids <- meta_lr$Study.ID[meta_lr$Location == "Right_colon"]  # n=130

# APC-only and APC+KRAS defined strictly (BRAF wild-type in both), matching
# Methods §Mutation Analysis. These four ID vectors are the basis of every
# mutation-progression comparison in Sections 10 onward.
a_all  <- rownames(mut_lr)[mut_lr$APC == 1 & mut_lr$KRAS == 0 & mut_lr$BRAF == 0]  # n=372
ak_all <- rownames(mut_lr)[mut_lr$APC == 1 & mut_lr$KRAS == 1 & mut_lr$BRAF == 0]  # n=162

length(ak_all)

a_left   <- intersect(a_all,  left_ids)   # n=320
ak_left  <- intersect(ak_all, left_ids)   # n=116
a_right  <- intersect(a_all,  right_ids)  # n=52
ak_right <- intersect(ak_all, right_ids)  # n=46

cat("Left:", length(left_ids), " Right:", length(right_ids), "\n")
cat("APC-only combined:", length(a_all), " APC+KRAS combined:", length(ak_all), "\n")
cat("APC-only left:", length(a_left), " right:", length(a_right), "\n")
cat("APC+KRAS left:", length(ak_left), " right:", length(ak_right), "\n")

p5_left  <- make_top10_bar(mut_lr, left_ids,  "Left colon",  "steelblue")
p5_right <- make_top10_bar(mut_lr, right_ids, "Right colon", "palevioletred")
# APC and KRAS excluded from the mutation-subgroup panels: they define the
# groups being compared, so including them would be circular.
p5_a     <- make_top10_bar(mut_lr, a_all,  "APC-only", unname(COL_MUT["APC_only"]), drop_genes = c("APC", "KRAS"))
p5_ak    <- make_top10_bar(mut_lr, ak_all, "APC+KRAS", unname(COL_MUT["APC_KRAS"]), drop_genes = c("APC", "KRAS"))

fig5 <- ggarrange(p5_left, p5_right, p5_a, p5_ak,
                  ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(file.path(PLOTS, "Figure5_top10_by_group.png"), fig5,
       width = 14, height = 12, dpi = 300)


# ------------------------------------------------------------------------------
# 8. MUTATION FREQUENCY PIE CHARTS (Figure 6 — adenoma subset)
# ------------------------------------------------------------------------------

# Shared helper: mutant/wild-type counts and percentages for one gene vector.
# Used both by the pie charts below and by the Figure 9B stacked bars.
make_mut_pie_df <- function(gene_vec) {
  df <- as.data.frame(table(ifelse(gene_vec == 1, "Mutant", "Wild-type")))
  colnames(df) <- c("Status", "Count")
  df$Percentage <- round(df$Count / sum(df$Count) * 100, 1)
  df
}

make_mut_pie <- function(gene_vec, title) {
  df <- make_mut_pie_df(gene_vec)
  
  ggplot(df, aes(x = "", y = Count, fill = Status)) +
    geom_bar(stat = "identity", width = 1, colour = "white") +
    coord_polar("y") +
    geom_text(aes(label = paste0(Percentage, "%")),
              position = position_stack(vjust = 0.5), size = 6, fontface = "bold", colour = "white") +
    scale_fill_manual(values = c("Mutant" = "palevioletred", "Wild-type" = "steelblue")) +
    labs(title = title) +
    theme_void(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          legend.text = element_text(size = 13), legend.title = element_text(size = 14))
}

p6_apc  <- make_mut_pie(mut_adenoma$APC,  "APC")
p6_kras <- make_mut_pie(mut_adenoma$KRAS, "KRAS")
p6_braf <- make_mut_pie(mut_adenoma$BRAF, "BRAF")

fig6 <- ggarrange(p6_apc, p6_kras, p6_braf, ncol = 3,
                  labels = c("A", "B", "C"), common.legend = TRUE, legend = "bottom")
ggsave(file.path(PLOTS, "Figure6_mutation_pies.png"), fig6,
       width = 13, height = 6, dpi = 300)


# ------------------------------------------------------------------------------
# 9. MUTATION COMBINATION SUBGROUPS (Figures 7 & 8)
# ------------------------------------------------------------------------------

# Figure 7 — overall adenoma subset
combo_df <- as.data.frame(table(mut_adenoma$combo))
colnames(combo_df) <- c("Combination", "Count")
combo_df$Percentage <- round(combo_df$Count / sum(combo_df$Count) * 100, 1)

fig7 <- ggplot(combo_df, aes(x = reorder(Combination, -Count), y = Percentage)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            vjust = -0.5, size = 4.5) +
  labs(x = "Mutation combination", y = "Percentage (%)") +
  theme_thesis() +
  ylim(0, max(combo_df$Percentage) * 1.15)
ggsave(file.path(PLOTS, "Figure7_mutation_combinations.png"), fig7,
       width = 10, height = 7, dpi = 300)

# Figure 8 — stratified by location (left / rectum / right)
meta_adenoma$Location <- factor(meta_adenoma$Location,
                                levels = c("Left_colon", "Rectum", "Right_colon"))

combo_loc <- mut_adenoma %>%
  rownames_to_column("Study.ID") %>%
  left_join(meta_adenoma %>% select(Study.ID, Location), by = "Study.ID") %>%
  filter(!is.na(Location)) %>%
  count(Location, combo) %>%
  group_by(Location) %>%
  mutate(Percentage = round(n / sum(n) * 100, 1)) %>%
  ungroup()

combo_order <- combo_df %>% arrange(desc(Count)) %>% pull(Combination) %>% as.character()
combo_loc$combo <- factor(combo_loc$combo, levels = combo_order)

fig8 <- ggplot(combo_loc, aes(x = combo, y = Percentage, fill = Location)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = paste0(Percentage, "%")),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3.2) +
  scale_fill_manual(values = c("Left_colon" = "steelblue", "Rectum" = "mediumpurple",
                               "Right_colon" = "palevioletred")) +
  labs(x = "Mutation combination", y = "Percentage (%)", fill = "Location") +
  theme_thesis(base_size = 13)
ggsave(file.path(PLOTS, "Figure8_mutation_by_location.png"), fig8,
       width = 13, height = 7, dpi = 300)


# ------------------------------------------------------------------------------
# 10. CLINICOPATHOLOGICAL COMPARISON — LEFT vs RIGHT (Figure 9, Table 1 text)
# ------------------------------------------------------------------------------

meta_lr_tbl <- meta_lr %>%
  mutate(
    APC  = factor(mut_lr$APC[match(Study.ID, rownames(mut_lr))],  levels = c(0,1), labels = c("Wild-type","Mutant")),
    KRAS = factor(mut_lr$KRAS[match(Study.ID, rownames(mut_lr))], levels = c(0,1), labels = c("Wild-type","Mutant")),
    BRAF = factor(mut_lr$BRAF[match(Study.ID, rownames(mut_lr))], levels = c(0,1), labels = c("Wild-type","Mutant")),
    Location = factor(Location, levels = c("Left_colon", "Right_colon"))
  )

table_lr <- meta_lr_tbl %>%
  select(Sex, Age, APC, KRAS, BRAF, Location) %>%
  tbl_summary(
    by = Location,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 1, missing = "no"
  ) %>%
  add_p(test = list(BRAF ~ "fisher.test")) %>%
  bold_labels() %>% bold_p()

table_lr %>% as_tibble() %>%
  write.csv(file.path(RES, "Table_LeftRight_clinicopathological.csv"), row.names = FALSE)

# Figure 9B — stacked bars, mutant vs wild-type proportion, Left vs Right
mut_stack_df <- bind_rows(
  data.frame(Gene = "APC",  Location = "Left_colon",  make_mut_pie_df(mut_lr$APC[match(left_ids, rownames(mut_lr))])),
  data.frame(Gene = "APC",  Location = "Right_colon", make_mut_pie_df(mut_lr$APC[match(right_ids, rownames(mut_lr))])),
  data.frame(Gene = "KRAS", Location = "Left_colon",  make_mut_pie_df(mut_lr$KRAS[match(left_ids, rownames(mut_lr))])),
  data.frame(Gene = "KRAS", Location = "Right_colon", make_mut_pie_df(mut_lr$KRAS[match(right_ids, rownames(mut_lr))])),
  data.frame(Gene = "BRAF", Location = "Left_colon",  make_mut_pie_df(mut_lr$BRAF[match(left_ids, rownames(mut_lr))])),
  data.frame(Gene = "BRAF", Location = "Right_colon", make_mut_pie_df(mut_lr$BRAF[match(right_ids, rownames(mut_lr))]))
)

fig9B <- ggplot(mut_stack_df, aes(x = Location, y = Percentage, fill = Status)) +
  geom_col() +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_stack(vjust = 0.5), size = 4.5, colour = "white") +
  facet_wrap(~Gene) +
  scale_fill_manual(values = c("Mutant" = "palevioletred", "Wild-type" = "steelblue")) +
  labs(x = "Polyp location", y = "Percentage (%)") +
  theme_thesis(base_size = 13)
ggsave(file.path(PLOTS, "Figure9B_mutation_stacked_bars.png"), fig9B,
       width = 11, height = 6, dpi = 300)


# ------------------------------------------------------------------------------
# 11. MUTATION FREQUENCY ACROSS THREE LOCATIONS (Bonferroni-corrected)
# ------------------------------------------------------------------------------

meta_adenoma_mut <- meta_adenoma %>%
  mutate(
    APC  = mut_adenoma$APC[match(Study.ID, rownames(mut_adenoma))],
    KRAS = mut_adenoma$KRAS[match(Study.ID, rownames(mut_adenoma))],
    BRAF = mut_adenoma$BRAF[match(Study.ID, rownames(mut_adenoma))]
  ) %>%
  filter(!is.na(Location))

p_apc_3loc  <- chisq.test(table(meta_adenoma_mut$Location, meta_adenoma_mut$APC))$p.value
p_kras_3loc <- chisq.test(table(meta_adenoma_mut$Location, meta_adenoma_mut$KRAS))$p.value
p_braf_3loc <- fisher.test(table(meta_adenoma_mut$Location, meta_adenoma_mut$BRAF),
                           simulate.p.value = TRUE)$p.value

p.adjust(c(APC = p_apc_3loc, KRAS = p_kras_3loc, BRAF = p_braf_3loc), method = "bonferroni")
# Reported in text: KRAS adjusted p = 5.6e-05, BRAF adjusted p = 0.022, APC adjusted p = 1.000


# ------------------------------------------------------------------------------
# 12. PRINCIPAL COMPONENT ANALYSIS (Figure — Supplementary)
# ------------------------------------------------------------------------------

gene_vars   <- apply(expr_adenoma, 1, var)
top3000     <- names(sort(gene_vars, decreasing = TRUE))[1:3000]
expr_top3000 <- expr_adenoma[top3000, ]

pca_meta <- data.frame(
  row.names = colnames(expr_top3000),
  Location  = meta_adenoma$Location[match(colnames(expr_top3000), meta_adenoma$Study.ID)],
  Combo     = mut_adenoma$combo[match(colnames(expr_top3000), rownames(mut_adenoma))]
)

p_obj <- pca(mat = expr_top3000, metadata = pca_meta, removeVar = 0.1)

# % variance explained by PC1 + PC2 — quoted in text as 7.19%
sum(p_obj$variance[1:2])

fig_pca <- biplot(p_obj, colby = "Location", shape = "Combo",
                  pointSize = 2.5, legendPosition = "right",
                  title = "PCA of Adenoma Gene Expression",
                  caption = paste0("n = ", ncol(expr_top3000), " adenoma samples"))
ggsave(file.path(PLOTS, "FigureS1_PCA.png"), fig_pca, width = 11, height = 8, dpi = 300)


# ------------------------------------------------------------------------------
# 13. GENE SETS FOR PATHWAY ANALYSIS
# ------------------------------------------------------------------------------

# MSigDB Hallmark collection (50 gene sets), via msigdbr.
hallmark_df   <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
length(hallmark_list)  # 50

# Cancer stem cell signatures (CBC, RSC, proCSC, revCSC): supervisor-provided
# gene lists (Gil Vazquez et al., 2022; Qin et al., 2023), supplied as
# stem_signatures.xlsx. Not derivable from a public database — the file
# itself must accompany this script for the stemness analyses to run.
stem_sigs <- readxl::read_excel(file.path(D, "stem_signatures.xlsx"))
stem_list <- list(
  CBC    = na.omit(stem_sigs$CBC),
  RSC    = na.omit(stem_sigs$RSC),
  proCSC = na.omit(stem_sigs$proCSC),
  revCSC = na.omit(stem_sigs$revCSC)
)
sapply(stem_list, length)

# Strips "HALLMARK_" prefix and converts underscores to spaces for display;
# stemness signature names (CBC, RSC, proCSC, revCSC) pass through unchanged.
clean_pathway_name <- function(pathway) {
  ifelse(grepl("^HALLMARK_", pathway),
         gsub("_", " ", sub("^HALLMARK_", "", pathway)),
         pathway)
}


# ------------------------------------------------------------------------------
# 14. ssGSEA ACROSS THREE LOCATIONS (Figures 10 & 11)
# ------------------------------------------------------------------------------
# Performed on all adenoma samples with recorded location (n=741) prior to
# the pairwise dualGSEA restriction to left/right (Section 16), following
# the rationale that dualGSEA cannot accommodate a three-level comparison.

meta_3loc <- meta_adenoma %>% filter(!is.na(Location))
expr_3loc <- expr_adenoma[, meta_3loc$Study.ID]

ssgsea_hallmark_param <- ssgseaParam(exprData = as.matrix(expr_3loc),
                                     geneSets = hallmark_list, normalize = TRUE)
ssgsea_hallmark <- gsva(ssgsea_hallmark_param, verbose = TRUE)

ssgsea_stem_param <- ssgseaParam(exprData = as.matrix(expr_3loc),
                                 geneSets = stem_list, normalize = TRUE)
ssgsea_stem <- gsva(ssgsea_stem_param, verbose = TRUE)

saveRDS(ssgsea_hallmark, file.path(D, "ssgsea_hallmark_3loc.rds"))
saveRDS(ssgsea_stem,     file.path(D, "ssgsea_stem_3loc.rds"))

# --- Kruskal-Wallis across all pathways within a gene set collection, with
#     BH correction across pathways, eta-squared effect size (rstatix), and
#     ranking to select the largest-effect pathways for display (Methods:
#     "the six with the largest effect sizes were selected"). ---

test_ssgsea_by_location <- function(ssgsea_mat, meta) {
  pathways <- rownames(ssgsea_mat)
  purrr::map_dfr(pathways, function(pw) {
    df <- data.frame(Score = ssgsea_mat[pw, ],
                     Location = meta$Location[match(colnames(ssgsea_mat), meta$Study.ID)])
    kw <- kruskal.test(Score ~ Location, data = df)
    eff <- rstatix::kruskal_effsize(df, Score ~ Location)
    data.frame(Pathway = pw, KW_stat = kw$statistic, p_value = kw$p.value,
               eta2 = eff$effsize)
  }) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(desc(eta2))
}

kw_hallmark <- test_ssgsea_by_location(ssgsea_hallmark, meta_3loc)
sum(kw_hallmark$p_adj < 0.05)  # reported in text as 23

kw_hallmark_sig <- kw_hallmark %>% filter(p_adj < 0.05)
write.csv(kw_hallmark, file.path(RES, "Table_ssGSEA_Hallmark_KruskalWallis_full.csv"), row.names = FALSE)

top6_pathways <- kw_hallmark_sig %>% arrange(desc(eta2)) %>% head(6) %>% pull(Pathway)
top6_pathways

kw_stem <- test_ssgsea_by_location(ssgsea_stem, meta_3loc)
write.csv(kw_stem, file.path(RES, "Table_ssGSEA_Stemness_KruskalWallis.csv"), row.names = FALSE)


# --- Reusable boxplot builder: KW result (top strip) + pairwise Wilcoxon
#     brackets (BH-adjusted within pathway), used for both Figure 10
#     (top 6 Hallmark pathways) and Figure 11 (four stemness signatures). ---

plot_ssgsea_boxplots <- function(ssgsea_mat, pathways, meta, kw_table,
                                 ncol = 3, base_size = 14) {
  long_df <- as.data.frame(t(ssgsea_mat[pathways, , drop = FALSE])) %>%
    rownames_to_column("Study.ID") %>%
    pivot_longer(-Study.ID, names_to = "Pathway", values_to = "Score") %>%
    left_join(meta %>% select(Study.ID, Location), by = "Study.ID") %>%
    mutate(Location = factor(Location, levels = c("Left_colon", "Rectum", "Right_colon")))
  
  # Build a combined facet label — "PATHWAY NAME\nKW p = X, eta2 = X" — so
  # both statistics appear in the strip title, matching the report figure.
  strip_labels <- kw_table %>%
    filter(Pathway %in% pathways) %>%
    mutate(clean_name = clean_pathway_name(Pathway),
           strip = paste0(toupper(clean_name), "\nKW p = ",
                          formatC(p_value, format = "e", digits = 2),
                          ", eta2 = ", round(eta2, 3))) %>%
    select(Pathway, strip) %>%
    tibble::deframe()
  
  long_df$Pathway <- factor(strip_labels[long_df$Pathway], levels = unname(strip_labels))
  
  ggplot(long_df, aes(x = Location, y = Score, fill = Location)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 0.8, colour = "black") +
    facet_wrap(~Pathway, scales = "free_y", ncol = ncol) +
    stat_compare_means(comparisons = list(c("Left_colon", "Rectum"),
                                          c("Left_colon", "Right_colon"),
                                          c("Rectum", "Right_colon")),
                       method = "wilcox.test", p.adjust.method = "BH",
                       label = "p.format", size = 4, tip.length = 0.01) +
    labs(x = NULL, y = "ssGSEA enrichment score") +
    theme_thesis(base_size = base_size) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

fig10 <- plot_ssgsea_boxplots(ssgsea_hallmark, top6_pathways, meta_3loc, kw_hallmark, ncol = 3)
ggsave(file.path(PLOTS, "Figure10_ssGSEA_top6_Hallmark.png"), fig10,
       width = 16, height = 10, dpi = 300)

fig11 <- plot_ssgsea_boxplots(ssgsea_stem, rownames(ssgsea_stem), meta_3loc, kw_stem, ncol = 4)
ggsave(file.path(PLOTS, "Figure11_ssGSEA_stemness.png"), fig11,
       width = 18, height = 6, dpi = 300)


# ==============================================================================
# 15. dualGSEA — SIX COMPARISONS (Figures 12, 14, 16–20)
# ==============================================================================
# Restricted to left/right colon adenomas throughout (Methods §dualGSEA
# pipeline): rectal adenomas excluded, and dualGSEA itself is
# limited to two-group comparisons.
#
# All six comparisons share the same publication-format barplot and boxplot
# functions defined immediately below, so a single change to font size,
# colour, or significance filtering applies to every one of Figures 12 and
# 14–21 at once.
# ------------------------------------------------------------------------------
# 15a. Publication-format plotting functions
# ------------------------------------------------------------------------------

# Barplot: NON-SIGNIFICANT PATHWAYS ARE DROPPED ENTIRELY (not greyed out),
# and the two enrichment directions are given distinct, high-contrast
# colours rather than a single "significant" colour 
plot_dualgsea_barplot_clean <- function(result_table, group_pos, group_neg,
                                        col_pos, col_neg,
                                        title = NULL, base_size = 15) {
  df <- result_table %>%
    filter(padj < 0.05) %>%
    mutate(
      # Stemness pathway names (CBC/RSC/proCSC/revCSC) pass through
      # unchanged; only the "HALLMARK_GENE_SET_NAME" style gets prefix
      # stripped and underscores converted to spaces.
      Pathway   = ifelse(grepl("^HALLMARK_", pathway),
                         gsub("_", " ", sub("^HALLMARK_", "", pathway)),
                         pathway),
      Direction = ifelse(NES > 0, group_pos, group_neg)
    ) %>%
    arrange(NES)
  df$Pathway <- factor(df$Pathway, levels = df$Pathway)
  
  ggplot(df, aes(x = Pathway, y = NES, fill = Direction)) +
    geom_col(width = 0.75) +
    geom_hline(yintercept = 0, colour = "grey40") +
    coord_flip() +
    scale_fill_manual(values = setNames(c(col_pos, col_neg), c(group_pos, group_neg)),
                      labels = DISPLAY_LABELS[c(group_pos, group_neg)]) +
    labs(title = title, x = NULL, y = "Normalised Enrichment Score (NES)", fill = NULL) +
    theme_thesis(base_size = base_size) +
    theme(legend.position = "top")
}

# Boxplot for a Hallmark comparison: the single most positive- and single
# most negative-NES SIGNIFICANT pathway, with ssGSEA scores drawn from the
# dualGSEA object's own single-sample result table (so the boxplot is
# guaranteed to reflect the same samples/scores as the barplot beside it).
plot_dualgsea_top_boxplots <- function(dualgsea_result, group_data, group_col,
                                       group_pos, group_neg, col_pos, col_neg,
                                       base_size = 15) {
  res_tbl <- dualgsea_result$Pairwise_ResultTable %>% filter(padj < 0.05)
  top_pos <- res_tbl %>% arrange(desc(NES)) %>% slice(1) %>% pull(pathway)
  top_neg <- res_tbl %>% arrange(NES)       %>% slice(1) %>% pull(pathway)
  
  ss_df <- dualgsea_result$SingleSample_ResultTable
  ss_df$Group <- group_data[[group_col]][match(ss_df$ID, rownames(group_data))]
  
  plot_df <- ss_df %>%
    select(Group, all_of(c(top_pos, top_neg))) %>%
    pivot_longer(-Group, names_to = "Pathway", values_to = "Score") %>%
    mutate(Pathway = gsub("HALLMARK_", "", Pathway), Pathway = gsub("_", " ", Pathway),
           Pathway = dplyr::recode(Pathway,
             "EPITHELIAL MESENCHYMAL TRANSITION" = "EMT",
             "TNFA SIGNALING VIA NFKB" = "TNFa"),
           Group = factor(Group, levels = c(group_pos, group_neg)))
  
  ggplot(plot_df, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(outlier.size = 0.8) +
    facet_wrap(~Pathway, scales = "free_y") +
    scale_fill_manual(values = setNames(c(col_pos, col_neg), c(group_pos, group_neg)),
                      labels = DISPLAY_LABELS[c(group_pos, group_neg)]) +
    scale_x_discrete(labels = DISPLAY_LABELS) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 5) +
    labs(x = NULL, y = "ssGSEA score") +
    theme_thesis(base_size = base_size) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
}

# Boxplot for a stemness comparison: ALL FOUR signatures shown (not just
# the top hit in each direction), matching Figures 13/15/19/21.
plot_dualgsea_stem_boxplots <- function(dualgsea_result, group_data, group_col,
                                        group_pos, group_neg, col_pos, col_neg,
                                        base_size = 15) {
  ss_df <- dualgsea_result$SingleSample_ResultTable
  ss_df$Group <- group_data[[group_col]][match(ss_df$ID, rownames(group_data))]
  
  plot_df <- ss_df %>%
    select(Group, CBC, RSC, proCSC, revCSC) %>%
    pivot_longer(-Group, names_to = "Signature", values_to = "Score") %>%
    mutate(Signature = factor(Signature, levels = c("CBC", "RSC", "proCSC", "revCSC")),
           Group = factor(Group, levels = c(group_pos, group_neg)))
  
  ggplot(plot_df, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(outlier.size = 0.8) +
    facet_wrap(~Signature, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = setNames(c(col_pos, col_neg), c(group_pos, group_neg)),
                      labels = DISPLAY_LABELS[c(group_pos, group_neg)]) +
    scale_x_discrete(labels = DISPLAY_LABELS) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 5) +
    labs(x = NULL, y = "ssGSEA score") +
    theme_thesis(base_size = base_size) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
}


# ------------------------------------------------------------------------------
# 15b. Comparison definitions
# ------------------------------------------------------------------------------
# Each entry supplies the sample IDs, the group data.frame dualGSEA expects,
# and the two group labels/colours used consistently across that
# comparison's barplot and boxplot(s).

build_group_df <- function(ids_pos, ids_neg, label_pos, label_neg) {
  data.frame(
    Group = factor(c(rep(label_pos, length(ids_pos)), rep(label_neg, length(ids_neg))),
                   levels = c(label_pos, label_neg)),
    row.names = c(ids_pos, ids_neg)
  )
}

comparisons <- list(
  AllAdenomas_LR = list(
    ids_pos = left_ids, ids_neg = right_ids,
    label_pos = "Left_colon", label_neg = "Right_colon",
    col_pos = unname(COL_LEFT_RIGHT["Left_colon"]), col_neg = unname(COL_LEFT_RIGHT["Right_colon"]),
    has_stem = TRUE,  fig_hallmark = "Figure12", fig_stem = "Figure13"
  ),
  APC_vs_APCKRAS_Combined = list(
    ids_pos = ak_all, ids_neg = a_all,
    label_pos = "APC_KRAS", label_neg = "APC_only",
    col_pos = unname(COL_MUT["APC_KRAS"]), col_neg = unname(COL_MUT["APC_only"]),
    has_stem = TRUE,  fig_hallmark = "Figure14", fig_stem = "Figure15"
  ),
  APConly_LR = list(
    ids_pos = a_left, ids_neg = a_right,
    label_pos = "Left_colon", label_neg = "Right_colon",
    col_pos = unname(COL_LEFT_RIGHT["Left_colon"]), col_neg = unname(COL_LEFT_RIGHT["Right_colon"]),
    has_stem = FALSE, fig_hallmark = "Figure16", fig_stem = NA
  ),
  APCKRAS_LR = list(
    ids_pos = ak_left, ids_neg = ak_right,
    label_pos = "Left_colon", label_neg = "Right_colon",
    col_pos = unname(COL_LEFT_RIGHT["Left_colon"]), col_neg = unname(COL_LEFT_RIGHT["Right_colon"]),
    has_stem = FALSE, fig_hallmark = "Figure17", fig_stem = NA
  ),
  APC_vs_APCKRAS_Left = list(
    ids_pos = ak_left, ids_neg = a_left,
    label_pos = "APC_KRAS", label_neg = "APC_only",
    col_pos = unname(COL_MUT["APC_KRAS"]), col_neg = unname(COL_MUT["APC_only"]),
    has_stem = TRUE,  fig_hallmark = "Figure18", fig_stem = "Figure19"
  ),
  APC_vs_APCKRAS_Right = list(
    ids_pos = ak_right, ids_neg = a_right,
    label_pos = "APC_KRAS", label_neg = "APC_only",
    col_pos = unname(COL_MUT["APC_KRAS"]), col_neg = unname(COL_MUT["APC_only"]),
    has_stem = TRUE,  fig_hallmark = "Figure20", fig_stem = "Figure21"
  )
)


# ------------------------------------------------------------------------------
# 15c. Run dualGSEA for every comparison (Hallmark, then stemness)
# ------------------------------------------------------------------------------

dualgsea_hallmark_results <- list()
dualgsea_stem_results     <- list()

for (cmp_name in names(comparisons)) {
  cmp <- comparisons[[cmp_name]]
  ids <- c(cmp$ids_pos, cmp$ids_neg)
  group_df <- build_group_df(cmp$ids_pos, cmp$ids_neg, cmp$label_pos, cmp$label_neg)
  
  cat("\n=== Running dualGSEA (Hallmark):", cmp_name, "===\n")
  res_h <- dualGSEA(
    data = as.data.frame(expr_lr[, ids]),
    group_data = group_df,
    group_colname = "Group",
    geneset_list = hallmark_list
  )
  dualgsea_hallmark_results[[cmp_name]] <- res_h
  saveRDS(res_h, file.path(D, paste0("dualgsea_hallmark_", cmp_name, ".rds")))
  
  bar_h <- plot_dualgsea_barplot_clean(
    res_h$Pairwise_ResultTable, cmp$label_pos, cmp$label_neg, cmp$col_pos, cmp$col_neg)
  box_h <- plot_dualgsea_top_boxplots(
    res_h, group_df, "Group", cmp$label_pos, cmp$label_neg, cmp$col_pos, cmp$col_neg)
  
  fig_h <- ggarrange(bar_h, box_h, ncol = 2, widths = c(1.3, 1), labels = c("A", "B"))
  ggsave(file.path(PLOTS, paste0(cmp$fig_hallmark, "_", cmp_name, "_Hallmark.png")),
         fig_h, width = 16, height = 8, dpi = 300)
  
  if (isTRUE(cmp$has_stem)) {
    cat("=== Running dualGSEA (Stemness):", cmp_name, "===\n")
    res_s <- dualGSEA(
      data = as.data.frame(expr_lr[, ids]),
      group_data = group_df,
      group_colname = "Group",
      geneset_list = stem_list
    )
    dualgsea_stem_results[[cmp_name]] <- res_s
    saveRDS(res_s, file.path(D, paste0("dualgsea_stem_", cmp_name, ".rds")))
    
    bar_s <- plot_dualgsea_barplot_clean(
      res_s$Pairwise_ResultTable, cmp$label_pos, cmp$label_neg, cmp$col_pos, cmp$col_neg)
    box_s <- plot_dualgsea_stem_boxplots(
      res_s, group_df, "Group", cmp$label_pos, cmp$label_neg, cmp$col_pos, cmp$col_neg)
    
    fig_s <- ggarrange(bar_s, box_s, ncol = 2, widths = c(1, 1.3), labels = c("A", "B"))
    ggsave(file.path(PLOTS, paste0(cmp$fig_stem, "_", cmp_name, "_Stemness.png")),
           fig_s, width = 16, height = 6, dpi = 300)
  }
}

# Convenience handles matching the thesis narrative directly (used in
# Section 16's Wilcoxon validation and Section 17's Venn diagrams).
result_location <- dualgsea_hallmark_results$AllAdenomas_LR
result_a         <- dualgsea_hallmark_results$APConly_LR
result_ak        <- dualgsea_hallmark_results$APCKRAS_LR
result_left_AvsAK  <- dualgsea_hallmark_results$APC_vs_APCKRAS_Left
result_right_AvsAK <- dualgsea_hallmark_results$APC_vs_APCKRAS_Right
result_AvsAK_LR    <- dualgsea_hallmark_results$APC_vs_APCKRAS_Combined

result_stem_location <- dualgsea_stem_results$AllAdenomas_LR
result_stem_left     <- dualgsea_stem_results$APC_vs_APCKRAS_Left
result_stem_right    <- dualgsea_stem_results$APC_vs_APCKRAS_Right
result_stem_AvsAK    <- dualgsea_stem_results$APC_vs_APCKRAS_Combined


# ------------------------------------------------------------------------------
# 16. WILCOXON VALIDATION OF ssGSEA SCORES (Table 2)
# ------------------------------------------------------------------------------
# Individual-sample-level validation of the four pathways highlighted in the
# "All adenomas Left vs Right" dualGSEA result (Section 15, result_location).

validate_pathway <- function(ss_result_table, group_data, group_col, pathway) {
  df <- data.frame(
    Score = ss_result_table[[pathway]],
    Group = group_data[[group_col]][match(ss_result_table$ID, rownames(group_data))]
  )
  wt <- wilcox.test(Score ~ Group, data = df)
  data.frame(Pathway = pathway, W = unname(wt$statistic), p_value = wt$p.value)
}

key_pathways <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                  "HALLMARK_E2F_TARGETS",
                  "HALLMARK_INTERFERON_ALPHA_RESPONSE")

group_df_lr <- build_group_df(left_ids, right_ids, "Left_colon", "Right_colon")

table2_wilcoxon <- bind_rows(lapply(key_pathways, function(pw)
  validate_pathway(result_location$SingleSample_ResultTable, group_df_lr, "Group", pw)
))
table2_wilcoxon
write.csv(table2_wilcoxon, file.path(RES, "Table2_Wilcoxon_validation.csv"), row.names = FALSE)


# ------------------------------------------------------------------------------
# 17. VENN DIAGRAMS (Figure 22, Table 3)
# ------------------------------------------------------------------------------
# Panels A/B: APC+KRAS-enriched and APC-only-enriched pathways, compared
# across the three mutation-progression analyses (combined, left, right).
# Panels C/D: Left-enriched and Right-enriched pathways, compared across the
# three location analyses (all adenomas, APC-only, APC+KRAS).

get_sig_by_direction <- function(result_table, direction = c("pos", "neg")) {
  direction <- match.arg(direction)
  result_table %>%
    filter(padj < 0.05, if (direction == "pos") NES > 0 else NES < 0) %>%
    pull(pathway)
}

# --- Mutation-progression sets (Panels A/B) ---
ak_combined_up <- get_sig_by_direction(result_AvsAK_LR$Pairwise_ResultTable, "pos")   # APC+KRAS-enriched
ak_left_up     <- get_sig_by_direction(result_left_AvsAK$Pairwise_ResultTable, "pos")
ak_right_up    <- get_sig_by_direction(result_right_AvsAK$Pairwise_ResultTable, "pos")

a_combined_up  <- get_sig_by_direction(result_AvsAK_LR$Pairwise_ResultTable, "neg")   # APC-only-enriched
a_left_up      <- get_sig_by_direction(result_left_AvsAK$Pairwise_ResultTable, "neg")
a_right_up     <- get_sig_by_direction(result_right_AvsAK$Pairwise_ResultTable, "neg")

venn_A <- list("Combined" = ak_combined_up, "Left" = ak_left_up, "Right" = ak_right_up)
venn_B <- list("Combined" = a_combined_up,  "Left" = a_left_up,  "Right" = a_right_up)

# --- Location sets (Panels C/D) ---
left_up_all    <- get_sig_by_direction(result_location$Pairwise_ResultTable, "pos")
left_up_a      <- get_sig_by_direction(result_a$Pairwise_ResultTable, "pos")
left_up_ak     <- get_sig_by_direction(result_ak$Pairwise_ResultTable, "pos")

right_up_all   <- get_sig_by_direction(result_location$Pairwise_ResultTable, "neg")
right_up_a     <- get_sig_by_direction(result_a$Pairwise_ResultTable, "neg")
right_up_ak    <- get_sig_by_direction(result_ak$Pairwise_ResultTable, "neg")

venn_C <- list("All adenomas" = left_up_all, "APC-only" = left_up_a, "APC+KRAS" = left_up_ak)
venn_D <- list("All adenomas" = right_up_all, "APC-only" = right_up_a, "APC+KRAS" = right_up_ak)

make_venn <- function(venn_list, title, fill_cols) {
  ggvenn(venn_list, fill_color = fill_cols, stroke_size = 0.5,
         set_name_size = 4.5, text_size = 4.5) +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
}

p_venn_A <- make_venn(venn_A, "APC+KRAS-enriched pathways", c("steelblue", "palevioletred", "mediumpurple"))
p_venn_B <- make_venn(venn_B, "APC-only-enriched pathways",  c("steelblue", "palevioletred", "mediumpurple"))
p_venn_C <- make_venn(venn_C, "Left colon-enriched pathways",  c("steelblue", "palevioletred", "mediumpurple"))
p_venn_D <- make_venn(venn_D, "Right colon-enriched pathways", c("steelblue", "palevioletred", "mediumpurple"))

fig22 <- ggarrange(p_venn_A, p_venn_B, p_venn_C, p_venn_D,
                   ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(file.path(PLOTS, "Figure22_Venn_diagrams.png"), fig22, width = 14, height = 12, dpi = 300)

# Table 3 — pathway identity per Venn region, for every panel.
venn_region_table <- function(venn_list, panel_label) {
  all_pw <- unique(unlist(venn_list))
  purrr::map_dfr(all_pw, function(pw) {
    membership <- names(venn_list)[sapply(venn_list, function(x) pw %in% x)]
    data.frame(Panel = panel_label, Pathway = pw,
               Present_in = paste(membership, collapse = "; "),
               N_analyses = length(membership))
  })
}

table3 <- bind_rows(
  venn_region_table(venn_A, "A (APC+KRAS-enriched)"),
  venn_region_table(venn_B, "B (APC-only-enriched)"),
  venn_region_table(venn_C, "C (Left-enriched)"),
  venn_region_table(venn_D, "D (Right-enriched)")
) %>% arrange(Panel, desc(N_analyses))

write.csv(table3, file.path(RES, "Table3_Venn_pathway_overlap.csv"), row.names = FALSE)


# ==============================================================================
# 18. TRANSCRIPTION FACTOR ACTIVITY INFERENCE — decoupleR (Figures 23–25)
# ==============================================================================
# CollecTRI regulon network (Müller-Dott et al., 2023), univariate linear
# model (run_ulm), minimum regulon size 5, following the reference workflow
# supplied for this analysis (GSE117606 script). Restricted to the same
# left/right colon adenoma subset used throughout dualGSEA.

net <- readRDS(file.path(D, "collectri_database.rds"))  # supervisor-provided network file

set.seed(121)
sample_acts <- decoupleR::run_ulm(
  mat = as.matrix(expr_lr), net = net,
  .source = "source", .target = "target", .mor = "mor", minsize = 5
)

sample_acts_mat <- sample_acts %>%
  filter(statistic == "ulm") %>%
  pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
  column_to_rownames("condition") %>%
  as.matrix()

dim(sample_acts_mat)  # expect 656 samples x 695 TFs
saveRDS(sample_acts_mat, file.path(D, "collectri_activity_matrix.rds"))

# Long format with mutation and location labels attached, for all subsequent
# TF-level testing (combined, left-only, right-only).
mut_labels_tf <- data.frame(
  Study.ID = c(a_all, ak_all),
  Mutation = factor(c(rep("APC-only", length(a_all)), rep("APC+KRAS", length(ak_all))),
                    levels = c("APC-only", "APC+KRAS"))
)

acts_long <- as.data.frame(sample_acts_mat) %>%
  rownames_to_column("Study.ID") %>%
  pivot_longer(-Study.ID, names_to = "TF", values_to = "Activity") %>%
  left_join(meta_lr %>% select(Study.ID, Location), by = "Study.ID") %>%
  left_join(mut_labels_tf, by = "Study.ID")

# --- Wilcoxon + BH testing, run identically for combined / left / right ---
test_tf_activity <- function(acts_long_df, restrict_location = NULL) {
  df <- acts_long_df %>% filter(!is.na(Mutation))
  if (!is.null(restrict_location)) df <- df %>% filter(Location == restrict_location)
  
  df %>%
    group_by(TF) %>%
    summarise(
      med_a  = median(Activity[Mutation == "APC-only"], na.rm = TRUE),
      med_ak = median(Activity[Mutation == "APC+KRAS"], na.rm = TRUE),
      p_value = tryCatch(wilcox.test(Activity ~ Mutation, exact = FALSE)$p.value,
                         error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(diff = med_ak - med_a, p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_adj)
}

res_mut       <- test_tf_activity(acts_long)                              # combined, n=534
res_mut_left  <- test_tf_activity(acts_long, restrict_location = "Left_colon")
res_mut_right <- test_tf_activity(acts_long, restrict_location = "Right_colon")

write.csv(res_mut,       file.path(RES, "collectri_TF_results_combined.csv"), row.names = FALSE)
write.csv(res_mut_left,  file.path(RES, "collectri_TF_results_left.csv"),     row.names = FALSE)
write.csv(res_mut_right, file.path(RES, "collectri_TF_results_right.csv"),    row.names = FALSE)

cat("TFs tested:", nrow(res_mut),
    "| significant (BH<0.05):", sum(res_mut$p_adj < 0.05, na.rm = TRUE), "\n")
cat("Higher in APC-only:", sum(res_mut$p_adj < 0.05 & res_mut$diff < 0, na.rm = TRUE),
    "| Higher in APC+KRAS:", sum(res_mut$p_adj < 0.05 & res_mut$diff > 0, na.rm = TRUE), "\n")

# --- Figure 23: volcano, top 5 labelled each direction ---
# Three-way colour: significant+higher in APC+KRAS, significant+higher in
# APC-only, and non-significant — matching Figure 23's description exactly
# ("blue... APC-only, red... APC+KRAS; grey... not reaching significance").
volc <- res_mut %>%
  mutate(Direction = case_when(
    p_adj < 0.05 & diff > 0 ~ "Higher in APC+KRAS",
    p_adj < 0.05 & diff < 0 ~ "Higher in APC-only",
    TRUE                    ~ "Not significant"
  ))

top_ak <- volc %>% filter(Direction == "Higher in APC+KRAS") %>% arrange(desc(diff)) %>% head(5) %>% pull(TF)
top_a  <- volc %>% filter(Direction == "Higher in APC-only")  %>% arrange(diff)       %>% head(5) %>% pull(TF)
top10_tfs <- c(top_ak, top_a)

fig23 <- ggplot(volc, aes(x = diff, y = -log10(p_adj), colour = Direction)) +
  geom_point(alpha = 0.7, size = 2.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = volc %>% filter(TF %in% top10_tfs),
                  aes(label = TF), size = 5, max.overlaps = 20, colour = "black") +
  scale_colour_manual(values = c("Higher in APC+KRAS" = unname(COL_MUT["APC_KRAS"]),
                                 "Higher in APC-only"  = unname(COL_MUT["APC_only"]),
                                 "Not significant"     = "grey80")) +
  labs(x = "Difference in median activity (APC+KRAS − APC-only)",
       y = expression(-log[10]("adjusted p-value")), colour = NULL) +
  theme_thesis(base_size = 15) +
  theme(legend.position = "bottom")
ggsave(file.path(PLOTS, "Figure23_TF_volcano.png"), fig23, width = 11, height = 8, dpi = 300)



#-----------Figure 24---------------
top10_tfs <- c(top_ak, top_a)

tf_order <- c(top_ak, top_a)  # FOXO3, LRRFIP1, REST, TGIF1, TBPL2, LMO2, GATA2, SP1, SPI1, ETV5
row_grp  <- factor(c(rep("Higher in APC+KRAS", length(top_ak)),
                     rep("Higher in APC-only", length(top_a))),
                   levels = c("Higher in APC-only", "Higher in APC+KRAS"))

heat_ids <- c(a_all, ak_all)  # restrict to the 534 samples used in TF testing

hm_mat   <- t(sample_acts_mat[heat_ids, tf_order])
hm_mat_z <- t(scale(t(hm_mat)))

col_grp <- factor(
  paste(mut_labels_tf$Mutation[match(heat_ids, mut_labels_tf$Study.ID)],
        meta_lr$Location[match(heat_ids, meta_lr$Study.ID)]),
  levels = c("APC-only Left_colon", "APC-only Right_colon",
             "APC+KRAS Left_colon", "APC+KRAS Right_colon"),
  labels = c("APC only Left", "APC only Right", "APC+KRAS Left", "APC+KRAS Right")
)

col_annotation <- HeatmapAnnotation(
  Location = meta_lr$Location[match(heat_ids, meta_lr$Study.ID)],
  Mutation = mut_labels_tf$Mutation[match(heat_ids, mut_labels_tf$Study.ID)],
  col = list(Mutation = c("APC-only" = unname(COL_MUT["APC_only"]), "APC+KRAS" = unname(COL_MUT["APC_KRAS"])),
             Location = c("Left_colon" = "seagreen3", "Right_colon" = "slateblue3"))
)

png(file.path(PLOTS, "Figure24_TF_heatmap.png"), width = 14, height = 6, units = "in", res = 300)
Heatmap(hm_mat_z, name = "TF activity\n(z-score)",
        top_annotation = col_annotation,
        column_split = col_grp,
        row_split = row_grp,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        cluster_column_slices = FALSE,
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        show_column_names = FALSE,
        row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        row_names_gp = grid::gpar(fontsize = 13),
        column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
        col = colorRamp2(c(-2, 0, 2), c(unname(COL_MUT["APC_only"]), "white", unname(COL_MUT["APC_KRAS"]))))
dev.off()


#----------Figure 25----------------

# Two separate TF sets
tf_apc_only <- top_a   # LMO2, GATA2, SP1, SPI1, ETV5
tf_apc_kras <- top_ak  # FOXO3, LRRFIP1, REST, TGIF1, TBPL2

# Shared p-value label builder (same as before, just called twice)
build_pval_labels <- function(tf_set) {
  bind_rows(
    res_mut_left  %>% filter(TF %in% tf_set) %>%
      transmute(TF, label = paste0("Left: ",  formatC(p_adj, format = "e", digits = 2))),
    res_mut_right %>% filter(TF %in% tf_set) %>%
      transmute(TF, label = paste0("Right: ", formatC(p_adj, format = "e", digits = 2)))
  ) %>%
    group_by(TF) %>%
    summarise(label = paste(label, collapse = "   "), .groups = "drop") %>%
    mutate(TF = factor(TF, levels = tf_set))
}

build_tf_boxplot <- function(tf_set, filename) {
  pval_labels <- build_pval_labels(tf_set)
  
  plot_df <- acts_long %>%
    filter(TF %in% tf_set, !is.na(Mutation)) %>%
    mutate(TF = factor(TF, levels = tf_set),
           Location = factor(Location, levels = c("Left_colon", "Right_colon"),
                             labels = c("Left", "Right")))
  
  p <- ggplot(plot_df, aes(x = Location, y = Activity, fill = Mutation)) +
    geom_boxplot(outlier.size = 1, width = 0.6) +
    geom_text(data = pval_labels, aes(x = 1.5, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.5, size = 4.2) +
    facet_wrap(~TF, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("APC-only" = unname(COL_MUT["APC_only"]), "APC+KRAS" = unname(COL_MUT["APC_KRAS"]))) +
    labs(x = NULL, y = "Inferred TF activity") +
    theme_thesis(base_size = 15) +
    theme(legend.position = "top",
          strip.text = element_text(size = 15, face = "bold"),
          axis.text.x = element_text(size = 13))
  
  ggsave(file.path(PLOTS, filename), p, width = 13, height = 11, dpi = 300)
  p
}

fig25a <- build_tf_boxplot(tf_apc_only, "Figure25a_TF_boxplots_APConly.png")
fig25b <- build_tf_boxplot(tf_apc_kras, "Figure25b_TF_boxplots_APCKRAS.png")

# ------------------------------------------------------------------------------
# 18b. TF ACTIVITY BY LOCATION — independent of mutation status
# ------------------------------------------------------------------------------

test_tf_activity_location <- function(acts_long_df) {
  acts_long_df %>%
    group_by(TF) %>%
    summarise(
      med_left  = median(Activity[Location == "Left_colon"],  na.rm = TRUE),
      med_right = median(Activity[Location == "Right_colon"], na.rm = TRUE),
      p_value = tryCatch(wilcox.test(Activity ~ Location, exact = FALSE)$p.value,
                         error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(diff = med_right - med_left, p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_adj)
}

res_location_tf <- test_tf_activity_location(acts_long)
write.csv(res_location_tf, file.path(RES, "collectri_TF_results_location.csv"), row.names = FALSE)

cat("TFs tested:", nrow(res_location_tf),
    "| significant (BH<0.05):", sum(res_location_tf$p_adj < 0.05, na.rm = TRUE), "\n")

# --- Volcano plot, same style as Figure 23, coloured by location instead of mutation ---
volc_loc <- res_location_tf %>%
  mutate(Direction = case_when(
    p_adj < 0.05 & diff > 0 ~ "Higher in Right colon",
    p_adj < 0.05 & diff < 0 ~ "Higher in Left colon",
    TRUE                    ~ "Not significant"
  ))

top_right <- volc_loc %>% filter(Direction == "Higher in Right colon") %>% arrange(desc(diff)) %>% head(5) %>% pull(TF)
top_left  <- volc_loc %>% filter(Direction == "Higher in Left colon")  %>% arrange(diff)       %>% head(5) %>% pull(TF)
top10_tfs_loc <- c(top_right, top_left)

fig23b <- ggplot(volc_loc, aes(x = diff, y = -log10(p_adj), colour = Direction)) +
  geom_point(alpha = 0.7, size = 2.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = volc_loc %>% filter(TF %in% top10_tfs_loc),
                  aes(label = TF), size = 5, max.overlaps = 20, colour = "black") +
  scale_colour_manual(values = c("Higher in Right colon" = unname(COL_LEFT_RIGHT["Right_colon"]),
                                 "Higher in Left colon"  = unname(COL_LEFT_RIGHT["Left_colon"]),
                                 "Not significant"       = "grey80")) +
  labs(x = "Difference in median activity (Right − Left)",
       y = expression(-log[10]("adjusted p-value")), colour = NULL) +
  theme_thesis(base_size = 15) +
  theme(legend.position = "bottom")
ggsave(file.path(PLOTS, "Figure23b_TF_volcano_location.png"), fig23b, width = 11, height = 8, dpi = 300)


# ------------------------------------------------------------------------------
# Master NES table — all Hallmark + stemness results, all six comparisons
# ------------------------------------------------------------------------------

comparison_names <- c("AllAdenomas_LR", "APC_vs_APCKRAS_Combined", "APConly_LR",
                      "APCKRAS_LR", "APC_vs_APCKRAS_Left", "APC_vs_APCKRAS_Right")
stem_comparisons <- c("AllAdenomas_LR", "APC_vs_APCKRAS_Combined",
                      "APC_vs_APCKRAS_Left", "APC_vs_APCKRAS_Right")


dualgsea_hallmark_results <- setNames(
  lapply(comparison_names, function(n) readRDS(file.path(D, paste0("dualgsea_hallmark_", n, ".rds")))),
  comparison_names
)
dualgsea_stem_results <- setNames(
  lapply(stem_comparisons, function(n) readRDS(file.path(D, paste0("dualgsea_stem_", n, ".rds")))),
  stem_comparisons
)

# Extract just the columns you actually need — pathway, NES, padj (no leadingEdge list)
extract_nes_table <- function(results_list, geneset_type) {
  bind_rows(lapply(names(results_list), function(cmp_name) {
    results_list[[cmp_name]]$Pairwise_ResultTable %>%
      transmute(Comparison = cmp_name, GeneSetType = geneset_type,
                Pathway = pathway, NES = round(NES, 3), padj = padj)
  }))
}

nes_full_master <- bind_rows(
  extract_nes_table(dualgsea_hallmark_results, "Hallmark"),
  extract_nes_table(dualgsea_stem_results, "Stemness")
) %>%
  arrange(GeneSetType, Comparison, padj)

write.csv(nes_full_master, file.path(RES, "Master_NES_all_comparisons_FULL.csv"), row.names = FALSE)
write.csv(nes_full_master %>% filter(padj < 0.05),
          file.path(RES, "Master_NES_all_comparisons_SIGNIFICANT.csv"), row.names = FALSE)

cat("Total rows:", nrow(nes_full_master), "| Significant (padj<0.05):", sum(nes_full_master$padj < 0.05), "\n")

# ------------------------------------------------------------------------------
# Master ssGSEA (3-location) table — Hallmark + stemness, Kruskal-Wallis results
# ------------------------------------------------------------------------------


kw_hallmark <- read.csv(file.path(RES, "Table_ssGSEA_Hallmark_KruskalWallis_full.csv"))
kw_stem     <- read.csv(file.path(RES, "Table_ssGSEA_Stemness_KruskalWallis.csv"))

ssgsea_3loc_master <- bind_rows(
  kw_hallmark %>% mutate(GeneSetType = "Hallmark"),
  kw_stem     %>% mutate(GeneSetType = "Stemness")
) %>%
  select(GeneSetType, Pathway, KW_stat, p_value, eta2, p_adj) %>%
  mutate(across(c(KW_stat, eta2), ~ round(., 3))) %>%
  arrange(GeneSetType, p_adj)

write.csv(ssgsea_3loc_master, file.path(RES, "Master_ssGSEA_3location_FULL.csv"), row.names = FALSE)
write.csv(ssgsea_3loc_master %>% filter(p_adj < 0.05),
          file.path(RES, "Master_ssGSEA_3location_SIGNIFICANT.csv"), row.names = FALSE)

