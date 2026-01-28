################################################################################
#
# Figure 3: Dynamic changes in adipocytes during cold adaptation
#
# Behrens et al. (2025) - Molecular Metabolism
# "Single-nucleus mRNA-sequencing reveals dynamics of lipogenic and thermogenic 
#  adipocyte populations in murine brown adipose tissue in response to cold exposure"
#
# This script generates the snRNA-seq panels of Figure 3:
#   - Panel A: Relative frequencies of adipocyte clusters by treatment (WT only)
#   - Panel B: Pathway scores heatmap across adipocyte subtypes
#   - Panel C: Apoptosis score violin plots for white-like and lipogenic adipocytes
#   - Panel D: Fold change of aggregated expression of Acaca and Fasn
#   - Panel F: Aggregated expression of Acaca and Fasn by adipocyte subtype
#   - Panels G-H: RNA velocity analysis (Acute Cold and Chronic Cold)
#
# Note: Panel E (qPCR validation) is generated from separate experiments 
#       and not included in this script
#
################################################################################

# ==============================================================================
# REQUIRED INPUT FILES
# ==============================================================================
# 
# 1. Main Seurat object (pre-processed, integrated, annotated):
#    - ChREBP_SeuratObject_final.rds
#    - Available from: doi: 10.25592/uhhfdm.18248
#
# 2. Supporting files (provided in the repository):
#    - data/color_assignment.rds
#    - data/rna_velocity_export_with_CellID_WT_AC.csv
#    - data/rna_velocity_export_WT_CC.csv
#    - AdlungLab_HelperFunctions.R
#
# ==============================================================================

# ==============================================================================
# SETUP: Load required packages
# ==============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(pheatmap)
library(data.table)
library(stringr)

# Load custom theme functions
# This file should contain theme_adlunglab() function and ALcols color palette
source("AdlungLab_HelperFunctions.R")

# Set working directory to the repository root
# setwd("path/to/ChREBP/")

# Define figure dimensions for publication
fontlist <- list("max" = 10, "min" = 7)
widthlist <- list("min" = 30, "singlecolumn" = 90, "oneandahalf" = 140, "fullwidth" = 190)

# ==============================================================================
# SECTION 1: LOAD DATA AND DEFINE CONSTANTS
# ==============================================================================
#
# Load the main Seurat object and define cell type mappings
#
# ==============================================================================

# --- Load main Seurat object ---
# Download from: doi: 10.25592/uhhfdm.18248
seurat_obj <- readRDS("ChREBP_SeuratObject_final.rds")

# Load color assignment for consistent plotting
color_assignment <- readRDS("data/color_assignment.rds")

# --- Define adipocyte subtypes ---
# These are the cell type labels used in the Seurat object
adipocyte_celltypes <- c(
  "1 Brown adipocyte (Ucp1 low)",
  "2 Brown adipocyte (Ucp1 high)",
  "4 White adipocyte",
  "8 Brown adipocyte (lipogenic)",
  "9 Brown adipocyte (from endothelial cell)",
  "10 Brown adipocyte (from stromal cell)",
  "16 Brown adipocyte (contractile)"
)

# Publication-ready adipocyte names (short form)
adipocyte_names <- c(
  "basal brown adipocytes",
  "OXPHOS-high adipocytes",
  "white-like adipocytes",
  "lipogenic adipocytes",
  "endothelium-derived adipocytes",
  "stroma-derived adipocytes",
  "contractile adipocytes"
)

# --- Define cell type colors for RNA velocity plots ---
celltype_colors <- c(
  "basal brown adipocytes"           = "#af6729",
  "OXPHOS-high brown adipocytes"     = "#876346",
  "stromal cells type 1"             = "#ff990a",
  "white-like adipocytes"            = "#dfc2a9",
  "fast muscle fibers"               = "#bf3335",
  "capillary endothelial cells"      = "#A0626A",
  "arterial endothelial cells"       = "#9e64a7",
  "lipogenic brown adipocytes"       = "#c79469",
  "endothelium-derived adipocytes"   = "#d4bad8",
  "stroma-derived adipocytes"        = "#dda8aa",
  "MonDC"                            = "#54A453",
  "smooth muscle cells"              = "#d8b642",
  "stromal cells type 2"             = "#ffad3a",
  "venous endothelial cells"         = "#94539e",
  "pericytes"                        = "#a98b24",
  "contractile brown adipocytes"     = "#af494b",
  "stromal cells type 3"             = "#ffcc84",
  "satellite glia cells"             = "#526E9F",
  "slow muscle fibers"               = "#af497e",
  "lymphocytes"                      = "#003200",
  "schwann cells"                    = "#3C8A9B"
)

# ==============================================================================
# SECTION 2: FIGURE 3A - ADIPOCYTE CLUSTER FREQUENCIES (WT ONLY)
# ==============================================================================
#
# Relative frequencies of adipocyte clusters across treatment conditions
# in Wild-Type (Cre-) mice only
#
# ==============================================================================

# --- Subset to Wild-Type (Cre-) mice ---
Idents(seurat_obj) <- "Genotype"
wt <- subset(seurat_obj, idents = "Wild-Type")

# --- Subset to adipocytes only ---
Idents(wt) <- "celltype"
adipo.behrens <- subset(wt, idents = adipocyte_celltypes)

message(paste0("Figure 3A: Wild-Type adipocyte subset contains ", 
               ncol(adipo.behrens), " nuclei"))

# --- Prepare data for plotting ---
cplot <- as.data.frame(adipo.behrens@meta.data)

# Set factor levels for Treatment
cplot$Treatment <- factor(
  cplot$Treatment, 
  levels = c("Room Temp.", "Acute Cold", "Chronic Cold"),
  labels = c("RT", "AC", "CC")
)

# --- Generate stacked bar plot ---
fig3a <- ggplot(cplot) +
  geom_bar(aes(x = Treatment, fill = celltype), 
           position = "fill", 
           color = "#000000",
           size = 0.1) +
  theme_adlunglab(base_size = 7) +
  scale_fill_manual(name = NULL, 
                    values = color_assignment, 
                    labels = adipocyte_names) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = NULL, 
       x = NULL, 
       y = "Relative frequency") +
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  guides(fill = guide_legend(ncol = 1))

print(fig3a)

# --- Save Panel A ---
fw <- widthlist$singlecolumn
ggsave("figures/Figure3A_AllAdipocytes.pdf",
       plot = fig3a,
       width = fw * 0.85, 
       height = fw * 4/3/2,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

# ==============================================================================
# SECTION 3: FIGURE 3B - PATHWAY SCORES HEATMAP
# ==============================================================================
#
# Heatmap showing metabolic pathway enrichment scores across adipocyte subtypes
# and treatment conditions. Scores are calculated using Seurat's AddModuleScore
# and normalized relative to basal adipocytes at room temperature.
#
# ==============================================================================

# --- Reload and prepare Seurat object ---
snuc <- readRDS("ChREBP_SeuratObject_final.rds")
seurat_obj <- JoinLayers(snuc)

# --- Subset to Wild-Type adipocytes ---
Idents(seurat_obj) <- "celltype"
seurat_obj <- subset(seurat_obj, idents = adipocyte_celltypes)

Idents(seurat_obj) <- "Genotype"
seurat_obj <- subset(seurat_obj, idents = "Wild-Type")

message(paste0("Figure 3B: Wild-Type adipocyte subset contains ", 
               ncol(seurat_obj), " nuclei, ", nrow(seurat_obj), " features"))

# --- Define metabolic gene sets ---
# Curated gene sets representing key metabolic pathways in brown adipose tissue
gene_sets <- list(
  `Oxphos regulators` = c("Ppargc1a", "Tfam"),
  `Ucp1` = c("Ucp1"),
  `Creatine futile cycling` = c("Alpl", "Ckb"),
  `Calcium futile cycling` = c("Atp2a1", "Atp2a2", "Ryr1", "Ryr2"),
  `Fatty acid transporters` = c("Cd36", "Slc27a1", "Slc27a2"),
  `Fatty acid binding proteins` = c("Fabp3", "Fabp4", "Fabp5"),
  `Carnitine plasma membrane transporters` = c("Slc22a4", "Slc22a5", "Slc22a21"),
  `Carnitine shuttle` = c("Slc25a20", "Cpt1b", "Cpt1a", "Cpt2"),
  `Fatty acid beta oxidation` = c("Acads", "Acadm", "Acadl", "Acadvl", "Echs1", 
                                   "Ehhadh", "Hadh", "Hadha", "Hadhb", "Acaa1a", "Acaa1b"),
  `OXPHOS, nucleus-encoded` = c("Cox4i1", "Cox8b", "Cox6c", "Cox7c", "Cox5a", 
                                 "Cox6b1", "Cox8a", "Cox7a1", "Cox6a1", "Cox5b", "Cox7a2"),
  `OXPHOS, mitochondria-encoded` = c("mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2", "mt-Atp8", 
                                      "mt-Atp6", "mt-Co3", "mt-Nd3", "mt-Nd4l", "mt-Nd4", 
                                      "mt-Nd5", "mt-Nd6", "mt-Cytb"),
  `DNL transcription factors` = c("Srebf1", "Mlxipl"),
  `Glucose transporters` = c("Slc2a1", "Slc2a3", "Slc2a4"),
  `Upper glycolysis` = c("Hk1", "Hk2", "Gpi1", "Pfkp", "Pfkl"),
  `Lower glycolysis` = c("Aldoa", "Tpi1", "Gapdh", "Pgk1", "Pgam1", "Pgam2", 
                          "Eno1", "Eno2", "Eno3", "Pkm", "Ldha", "Ldhb"),
  `Pentose phosphate pathway` = c("G6pdx", "Pgls", "Pgd", "Rpia", "Rpe", "Tkt", "Taldo1"),
  `De novo lipogenesis` = c("Me1", "Slc25a1", "Acly", "Acss2", "Aacs", "Acat2", 
                             "Acaca", "Acacb", "Elovl6", "Scd1", "Scd2", "Fasn"),
  `Triglyceride synthesis` = c("Gpd1", "Gpd2", "Gpam", "Gpat3", "Gpat4", "Agpat1", 
                                "Agpat2", "Agpat3", "Agpat4", "Agpat5", "Lpin1", 
                                "Lpin2", "Lpin3", "Dgat1", "Dgat2"),
  `Acyl CoA synthesis` = c("Acsl1", "Acsl3", "Acsl4", "Acsl5"),
  `Lipolysis` = c("Lipe", "Pnpla2", "Mgll"),
  `Pyruvate dehydrogenase` = c("Mpc1", "Mpc2", "Pdha1", "Pdhb", "Dlat", "Dld"),
  `Citric acid cycle` = c("Cs", "Aco2", "Idh3a", "Idh3g", "Ogdh", "Sucla2", 
                           "Suclg1", "Suclg2", "Sdha", "Sdhb", "Fh1", "Mdh1", "Mdh2")
)

# --- Calculate module scores for each pathway ---
for (pathway in names(gene_sets)) {
  seurat_obj <- AddModuleScore(
    seurat_obj, 
    features = list(gene_sets[[pathway]]), 
    name = pathway
  )
}

# --- Recode treatment labels for proper ordering ---
seurat_obj@meta.data$Treatment <- recode(
  seurat_obj@meta.data$Treatment,
  "Room Temp." = "01RT",
  "Acute Cold" = "02AC",
  "Chronic Cold" = "03CC"
)

# --- Extract scores and metadata ---
# Select the module score columns (last N columns, where N = number of gene sets)
scores <- seurat_obj@meta.data %>%
  dplyr::select(tail(seq_along(.), length(gene_sets))) %>%
  # Remove "_1" suffix added by Seurat's AddModuleScore
  rename_with(~ gsub("1$", "", .))

scores$CellType <- seurat_obj@meta.data$celltype
scores$Treatment <- seurat_obj@meta.data$Treatment

# --- Aggregate scores by cell type and treatment ---
agg_scores <- scores %>%
  group_by(CellType, Treatment) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(
    -c(CellType, Treatment), 
    names_to = "GeneSet", 
    values_to = "Score"
  )

# --- Create matrix for heatmap ---
heatmap_matrix <- agg_scores %>%
  pivot_wider(names_from = c(Treatment, CellType), values_from = Score) %>%
  as.data.frame()

# Set row names and remove the GeneSet column
rownames(heatmap_matrix) <- heatmap_matrix$GeneSet
heatmap_matrix <- heatmap_matrix[, -1]

# --- Create annotation data frame ---
col_names <- colnames(heatmap_matrix)

annotations <- data.frame(
  Treatment = sapply(strsplit(col_names, "_"), `[`, 1),
  CellType = sapply(strsplit(col_names, "_"), `[`, 2)
)
rownames(annotations) <- col_names

# Set factor levels for proper ordering
annotations$Treatment <- factor(
  annotations$Treatment, 
  levels = c("01RT", "02AC", "03CC"), 
  labels = c("Room Temp.", "Acute Cold", "Chronic Cold")
)

# Define publication labels for adipocyte types
labels.adipocytes <- c(
  "basal brown adipocytes", 
  "oxphos-high brown adipocytes",
  "white-like adipocytes", 
  "lipogenic adipocytes", 
  "endothelium-derived adipocytes", 
  "stromal-derived adipocytes", 
  "contractile brown adipocytes"
)

annotations$CellType <- factor(
  annotations$CellType, 
  levels = adipocyte_celltypes,
  labels = labels.adipocytes
)

# --- Normalize to first column (RT basal adipocytes) ---
reference_column <- heatmap_matrix[, 1]
heatmap_matrix_scaled <- heatmap_matrix - reference_column

# --- Create symmetric color breaks centered at 0 ---
max_val <- max(abs(heatmap_matrix_scaled), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = 51)

# --- Generate heatmap ---
# Save heatmap to PDF
pdf("figures/Figure3B_PathwayScoresHeatmap.pdf", 
    width = 380 / 25.4, 
    height = (380 * 3 / 4) / 25.4)

pheatmap(
  heatmap_matrix_scaled,
  annotation_col = annotations,
  annotation_colors = list(
    CellType = c(
      "basal brown adipocytes" = "#af6729",
      "oxphos-high brown adipocytes" = "#876346",
      "white-like adipocytes" = "#dfc2a9",
      "lipogenic adipocytes" = "#c79469",
      "endothelium-derived adipocytes" = "#d4bad8", 
      "stromal-derived adipocytes" = "#dda8aa",
      "contractile brown adipocytes" = "#af494b"
    ),
    Treatment = c(
      "Room Temp." = ALcols[1],
      "Acute Cold" = ALcols[3],
      "Chronic Cold" = ALcols[5]
    )
  ),
  cluster_rows = TRUE,
  cluster_cols = FALSE, 
  show_colnames = FALSE,
  scale = "none",
  color = colorRampPalette(c("blue", "white", "red"))(50), 
  breaks = breaks
)

dev.off()

message("Figure 3B heatmap saved to figures/Figure3B_PathwayScoresHeatmap.pdf")

# ==============================================================================
# SECTION 4: FIGURE 3C - APOPTOSIS SCORE
# ==============================================================================
#
# Violin plots showing apoptosis pathway scores in white-like and lipogenic
# adipocytes across treatment conditions (Wild-Type only)
#
# ==============================================================================

# --- Reload Seurat object ---
snuc <- readRDS("ChREBP_SeuratObject_final.rds")
snuc <- JoinLayers(snuc)

# --- Define apoptosis gene set ---
# Gene set from WikiPathways WP_APOPTOSIS (https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/WP_APOPTOSIS.html)
apoptosis_genes <- list(
  c("Akt1", "Apaf1", "Birc3", "Birc2", "Xiap", "Birc5", "Bad", "Bax", "Bcl2", "Bcl2l1", 
    "Bcl2l2", "Bid", "Hrk", "Bcl2l11", "Bnip3l", "Casp1", "Casp4", "Casp2", "Casp3", 
    "Casp6", "Casp7", "Casp8", "Casp9", "Cflar", "Chuk", "Cradd", "Dffa", "Dffb", "Fas", 
    "Fasl", "Gzmc", "Hells", "Igf1", "Igf1r", "Igf2", "Ikbkb", "Ikbkg", "Irf1", "Irf2", 
    "Irf4", "Jun", "Lta", "Mcl1", "Mdm2", "Myc", "Nfkb1", "Nfkbia", "Nfkbib", "Nfkbie", 
    "Prf1", "Pik3r1", "Rela", "Ripk1", "Tnf", "Tnfrsf10b", "Tnfrsf1a", "Tnfrsf1b", 
    "Traf1", "Traf2", "Traf3", "Tnfsf10", "Trp53", "Trp63", "Trp73", "Map2k4", "Map3k1", 
    "Mapk10", "Irf5", "Bok", "Irf7", "Irf3", "Irf6", "Pmaip1", "Diablo", "Tradd", 
    "Scaf11", "Tnfrsf25", "Tnfrsf21")
)

# --- Calculate apoptosis module score ---
snuc <- AddModuleScore(snuc, features = apoptosis_genes, name = "Death")

# --- Prepare data for plotting ---
cplot <- as.data.frame(snuc@meta.data)

# Set factor levels
cplot$Genotype <- factor(
  cplot$Genotype, 
  levels = c("Wild-Type", "Knock-Out"),
  labels = c("Cre-", "Cre+")
)

cplot$Treatment <- factor(
  cplot$Treatment, 
  levels = c("Room Temp.", "Acute Cold", "Chronic Cold"),
  labels = c("RT", "AC", "CC")
)

# --- Subset to white-like and lipogenic adipocytes (WT only) ---
cplot3 <- cplot %>%
  dplyr::filter(
    celltype %in% c("4 White adipocyte", "8 Brown adipocyte (lipogenic)"),
    Genotype == "Cre-"
  )

cplot3$celltype <- factor(
  cplot3$celltype, 
  levels = c("4 White adipocyte", "8 Brown adipocyte (lipogenic)"),
  labels = c("white-like adipocytes", "lipogenic adipocytes")
)

# --- Generate violin plot ---
fig3c <- ggplot(cplot3, aes(x = Treatment, y = Death1, fill = Treatment)) + 
  geom_violin(color = "#000000", 
              draw_quantiles = c(0.5), 
              trim = TRUE, 
              scale = "width",
              width = 0.66) + 
  labs(title = NULL,
       x = NULL, 
       y = "Apoptosis score\n(normalized expression)") +
  scale_fill_manual(values = c("RT" = ALcols[1],
                                "AC" = ALcols[3],
                                "CC" = ALcols[5]), 
                    guide = "none") +
  theme_adlunglab(base_size = fontlist$min) +
  facet_wrap(~celltype, ncol = 6) +
  theme(strip.text = element_text(size = fontlist$min))

# --- Add statistical comparisons ---
# Perform Wilcoxon test with Bonferroni correction
stat_test <- compare_means(
  Death1 ~ Treatment, 
  data = cplot3, 
  group.by = c("celltype"),
  method = "wilcox.test",
  p.adjust.method = "bonferroni"
)

print(stat_test)

# Add significance annotations to plot
fig3c_stats <- fig3c + stat_compare_means(
  comparisons = list(c("RT", "AC"), c("RT", "CC"), c("AC", "CC")), 
  label = "p.signif", 
  method = "wilcox.test", 
  p.adjust.method = "bonferroni",
  hide.ns = FALSE, 
  size = 2
)

print(fig3c_stats)

# --- Save Panel C ---
fw <- widthlist$singlecolumn
ggsave("figures/Figure3C_Apoptosis_Score_WhiteLipo.pdf",
       plot = fig3c_stats,
       width = fw, 
       height = fw * 3/4,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

# ==============================================================================
# SECTION 5: FIGURE 3D - AGGREGATED DNL GENE EXPRESSION (FOLD CHANGE)
# ==============================================================================
#
# Fold change of aggregated expression of key DNL genes (Acaca, Fasn)
# across samples relative to Wild-Type Room Temperature
#
# ==============================================================================

# --- Reload and subset Seurat object ---
seurat_obj <- readRDS("ChREBP_SeuratObject_final.rds")
seurat_obj <- JoinLayers(seurat_obj)

Idents(seurat_obj) <- "celltype"
adipo.behrens <- subset(seurat_obj, idents = adipocyte_celltypes)

Idents(adipo.behrens) <- "Genotype"
adipo.behrens <- subset(adipo.behrens, idents = "Wild-Type")

message(paste0("Figure 3D: Wild-Type adipocyte subset contains ", 
               ncol(adipo.behrens), " nuclei"))

# --- Define genes of interest ---
genes <- c("Acaca", "Fasn")

# --- Aggregate expression by sample ---
agg_list <- AggregateExpression(
  adipo.behrens, 
  features = genes, 
  group.by = "sample", 
  assays = "RNA", 
  slot = "counts"
)

# Extract and transform data
agg_expr <- agg_list$RNA
df_agg <- as.data.frame(t(agg_expr))
df_agg$sample <- rownames(df_agg)

# --- Convert to long format ---
df_long <- df_agg %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "gene",
    values_to = "agg_expression"
  )

# --- Calculate fold change relative to WT Room Temperature ---
df_long <- df_long %>%
  group_by(gene) %>%
  mutate(
    ref_expression = agg_expression[sample == "g01-Wild-Type-Room Temp."],
    fold_change = agg_expression / ref_expression
  ) %>%
  ungroup()

# --- Generate bar plot ---
fig3d <- ggplot(df_long, aes(x = gene, group = sample, y = fold_change, 
                              color = sample, fill = sample)) +
  geom_col(position = "dodge") +
  theme_adlunglab(base_size = 7) +
  labs(title = NULL,
       x = NULL,
       y = "Fold expression",
       fill = "Sample",
       color = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(fig3d)

# --- Save Panel D ---
ggsave("figures/Figure3D_DNL_FoldChange.pdf",
       plot = fig3d,
       width = widthlist$singlecolumn, 
       height = widthlist$singlecolumn * 3/4,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

# ==============================================================================
# SECTION 6: FIGURE 3F - DNL GENES BY ADIPOCYTE SUBTYPE
# ==============================================================================
#
# Aggregated expression of Fasn and Acaca across adipocyte subtypes
# and treatment conditions (Wild-Type only)
#
# ==============================================================================

# --- Reload Seurat object ---
seurat_obj <- readRDS("ChREBP_SeuratObject_final.rds")
seurat_obj <- JoinLayers(seurat_obj)

# --- Subset to Wild-Type adipocytes ---
Idents(seurat_obj) <- "celltype"
adipo.behrens <- subset(seurat_obj, idents = adipocyte_celltypes)

Idents(adipo.behrens) <- "Genotype"
adipo.behrens <- subset(adipo.behrens, idents = "Wild-Type")

# --- Create combined grouping variable ---
adipo.behrens$celltype_Treatment <- paste(
  adipo.behrens$celltype, 
  adipo.behrens$Treatment, 
  sep = "_"
)

# --- Aggregate expression by celltype and treatment ---
genes <- c("Fasn", "Acaca")

agg_list <- AggregateExpression(
  adipo.behrens, 
  features = genes, 
  group.by = "celltype_Treatment", 
  assays = "RNA", 
  slot = "counts"
)

agg_expr <- agg_list$RNA
df_agg <- as.data.frame(t(agg_expr))
df_agg$celltype_Treatment <- rownames(df_agg)

# --- Convert to long format ---
df_long <- df_agg %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "gene",
    values_to = "agg_expression"
  )

# --- Split combined column back into celltype and Treatment ---
df_long <- df_long %>%
  separate(celltype_Treatment, into = c("celltype", "Treatment"), sep = "-")

# --- Fix celltype factor levels ---
# Note: The split creates "g1", "g2", etc. prefixes that need to be handled
df_long$celltype <- factor(
  df_long$celltype, 
  levels = c("g1 Brown adipocyte (Ucp1 low)",
             "g2 Brown adipocyte (Ucp1 high)",
             "g4 White adipocyte",
             "g8 Brown adipocyte (lipogenic)",
             "g9 Brown adipocyte (from endothelial cell)",
             "g10 Brown adipocyte (from stromal cell)",
             "g16 Brown adipocyte (contractile)"),
  labels = c("1 Brown adipocyte (Ucp1 low)",
             "2 Brown adipocyte (Ucp1 high)",
             "4 White adipocyte",
             "8 Brown adipocyte (lipogenic)",
             "9 Brown adipocyte (from endothelial cell)",
             "10 Brown adipocyte (from stromal cell)",
             "16 Brown adipocyte (contractile)")
)

df_long$Treatment <- factor(
  df_long$Treatment, 
  levels = c("Room Temp.", "Acute Cold", "Chronic Cold"),
  labels = c("RT", "AC", "CC")
)

# --- Generate faceted bar plot ---
fig3f <- ggplot(df_long, aes(x = Treatment, y = agg_expression, fill = celltype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(
    values = color_assignment, 
    labels = c("basal brown adipocytes",
               "oxphos-high brown adipocytes",
               "white-like adipocytes",
               "lipogenic adipocytes",
               "endothelium-derived adipocytes",
               "stroma-derived adipocytes",
               "contractile adipocytes"),
    name = NULL
  ) +
  theme_adlunglab(base_size = 10) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = NULL,
       y = "Aggregated expression",
       x = NULL) +
  theme(legend.position = "right")

print(fig3f)

# --- Save Panel F ---
fw <- widthlist$fullwidth
ggsave("figures/Figure3F_DNL_SingleGenes_AggregatedExpression.pdf",
       plot = fig3f,
       width = fw, 
       height = fw * 9/16/1.5,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

# ==============================================================================
# SECTION 7: FIGURE 3G - RNA VELOCITY (ACUTE COLD)
# ==============================================================================
#
# Violin plot showing RNA velocity magnitudes across cell types
# for Wild-Type mice under Acute Cold (AC) conditions
#
# RNA velocity data was pre-computed using scVelo and exported to CSV
#
# ==============================================================================

# --- Reload Seurat object and subset to WT Acute Cold ---
seurat_obj <- readRDS("ChREBP_SeuratObject_final.rds")

Idents(seurat_obj) <- "id"
seurat_obj <- subset(seurat_obj, idents = c("02"))  # WT Acute Cold sample ID

message(paste0("Figure 3G: WT Acute Cold subset contains ", 
               ncol(seurat_obj), " nuclei"))

# --- Load pre-computed RNA velocity data ---
velocity_data <- fread("data/rna_velocity_export_with_CellID_WT_AC.csv")

# --- Align CellIDs between velocity data and Seurat object ---
# Fix CellIDs: Remove prefix and replace 'x' with '-1' to match Seurat format
velocity_data[, CellID := str_replace(CellID, "^.*:", "")]
velocity_data[, CellID := str_replace(CellID, "x$", "-1")]

# Extract the sample prefix from Seurat object column names
seurat_prefix <- unique(str_extract(colnames(seurat_obj), "^.*_"))

# Append prefix to velocity CellIDs
velocity_data[, CellID := paste0(seurat_prefix, CellID)]

# Check alignment
matching_cells <- sum(velocity_data$CellID %in% colnames(seurat_obj))
total_cells <- nrow(velocity_data)
message(sprintf("Matching CellIDs: %d/%d", matching_cells, total_cells))

# --- Keep only matching cells ---
matching_cells <- intersect(colnames(seurat_obj), velocity_data$CellID)
velocity_data <- velocity_data[CellID %in% matching_cells]
seurat_obj <- subset(seurat_obj, cells = matching_cells)

# Ensure ordering matches
velocity_data <- velocity_data[match(colnames(seurat_obj), velocity_data$CellID), ]

# --- Merge with UMAP coordinates ---
umap_data <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
umap_data$CellID <- rownames(umap_data)

velocity_data <- merge(velocity_data, umap_data, by = "CellID")

# --- Add cell type information ---
velocity_data[, CellType := seurat_obj$celltype[match(velocity_data$CellID, colnames(seurat_obj))]]

# Verify alignment
if (!all(velocity_data$CellID == colnames(seurat_obj))) {
  stop("Error: CellIDs are not aligned between velocity_data and Seurat object!")
}

# --- Add velocity data to Seurat object ---
rownames(velocity_data) <- velocity_data$CellID
velocity_data$CellID <- NULL
seurat_obj <- AddMetaData(seurat_obj, metadata = velocity_data)

# --- Calculate velocity magnitude ---
seurat_obj$velocity_magnitude <- sqrt(seurat_obj$Velocity1^2 + seurat_obj$Velocity2^2)

# --- Prepare data for plotting ---
df_vel <- seurat_obj@meta.data %>%
  rename(
    Velocity = velocity_magnitude,
    Celltype = cellname
  )

# --- Generate violin plot ---
fig3g <- ggplot(df_vel, aes(x = Velocity, y = Celltype, fill = Celltype)) +
  geom_violin(trim = TRUE, scale = "width", show.legend = FALSE) +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  theme_adlunglab(base_size = 7) +
  labs(x = "RNA velocity absolute",
       y = "Cell type") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8)) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 1.5,
    colour = "white",
    stroke = 0.5
  )

print(fig3g)

# ==============================================================================
# SECTION 8: FIGURE 3H - RNA VELOCITY (CHRONIC COLD)
# ==============================================================================
#
# Violin plot showing RNA velocity magnitudes across cell types
# for Wild-Type mice under Chronic Cold (CC) conditions
#
# ==============================================================================

# --- Reload Seurat object and subset to WT Chronic Cold ---
seurat_obj <- readRDS("ChREBP_SeuratObject_final.rds")

Idents(seurat_obj) <- "id"
seurat_obj <- subset(seurat_obj, idents = c("03"))  # WT Chronic Cold sample ID

message(paste0("Figure 3H: WT Chronic Cold subset contains ", 
               ncol(seurat_obj), " nuclei"))

# --- Load pre-computed RNA velocity data ---
velocity_data <- fread("data/rna_velocity_export_WT_CC.csv")

# --- Align CellIDs (same process as Section 7) ---
velocity_data[, CellID := str_replace(CellID, "^.*:", "")]
velocity_data[, CellID := str_replace(CellID, "x$", "-1")]

seurat_prefix <- unique(str_extract(colnames(seurat_obj), "^.*_"))
velocity_data[, CellID := paste0(seurat_prefix, CellID)]

matching_cells <- sum(velocity_data$CellID %in% colnames(seurat_obj))
total_cells <- nrow(velocity_data)
message(sprintf("Matching CellIDs: %d/%d", matching_cells, total_cells))

# --- Keep only matching cells ---
matching_cells <- intersect(colnames(seurat_obj), velocity_data$CellID)
velocity_data <- velocity_data[CellID %in% matching_cells]
seurat_obj <- subset(seurat_obj, cells = matching_cells)

velocity_data <- velocity_data[match(colnames(seurat_obj), velocity_data$CellID), ]

# --- Merge with UMAP and cell type info ---
umap_data <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
umap_data$CellID <- rownames(umap_data)

velocity_data <- merge(velocity_data, umap_data, by = "CellID")
velocity_data[, CellType := seurat_obj$celltype[match(velocity_data$CellID, colnames(seurat_obj))]]

if (!all(velocity_data$CellID == colnames(seurat_obj))) {
  stop("Error: CellIDs are not aligned between velocity_data and Seurat object!")
}

rownames(velocity_data) <- velocity_data$CellID
velocity_data$CellID <- NULL
seurat_obj <- AddMetaData(seurat_obj, metadata = velocity_data)

# --- Calculate velocity magnitude ---
seurat_obj$velocity_magnitude <- sqrt(seurat_obj$Velocity1^2 + seurat_obj$Velocity2^2)

# --- Prepare data for plotting ---
df_vel <- seurat_obj@meta.data %>%
  rename(
    Velocity = velocity_magnitude,
    Celltype = cellname
  )

# --- Generate violin plot ---
fig3h <- ggplot(df_vel, aes(x = Velocity, y = Celltype, fill = Celltype)) +
  geom_violin(trim = TRUE, scale = "width", show.legend = FALSE) +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  theme_adlunglab(base_size = 7) +
  labs(x = "RNA velocity absolute",
       y = "Cell type") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8)) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 1.5,
    colour = "white",
    stroke = 0.5
  )

print(fig3h)

# ==============================================================================
# SECTION 9: COMBINE AND SAVE RNA VELOCITY PANELS
# ==============================================================================

# --- Combine Panels G and H ---
scV <- ggarrange(
  fig3g, 
  fig3h, 
  nrow = 1, 
  labels = c("G", "H"),
  widths = c(1/2, 1/2)
)

print(scV)

# --- Save combined velocity figure ---
fw <- widthlist$fullwidth
ggsave("figures/Figure3GH_scVelo_Violin.pdf",
       plot = scV,
       width = fw, 
       height = fw * 9/16/1.5,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

ggsave("figures/Figure3GH_scVelo_Violin.png",
       plot = scV,
       width = fw, 
       height = fw * 9/16/1.5,
       units = "mm", 
       dpi = 300)

message("All Figure 3 panels saved to figures/ directory")

# ==============================================================================
# SESSION INFO
# ==============================================================================

sessionInfo()

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
