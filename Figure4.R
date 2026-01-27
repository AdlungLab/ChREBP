################################################################################
#
# Figure 4: Lipogenic adipocyte identity depends on ChREBP and is linked to 
#           Wnt signaling
#
# Behrens et al. (2025) - Molecular Metabolism
# "Single-nucleus mRNA-sequencing reveals dynamics of lipogenic and thermogenic 
#  adipocyte populations in murine brown adipose tissue in response to cold exposure"
#
# This script generates the snRNA-seq panels of Figure 4:
#   - Panel C: Fold change of aggregated expression of Mlxipl, Acaca, and Fasn
#   - Panel E: Relative frequencies of adipocyte clusters by genotype/treatment
#   - Panel F: PROGENy pathway activity heatmap across adipocyte subtypes
#   - Panel G: WNT PROGENy scores in lipogenic adipocytes (Cre- vs Cre+)
#   - Panel H: Aggregated expression of Wnt ligands across treatments
#   - Panel I: LIANA cell-cell interaction analysis (Wnt ligand-receptor)
#   - Panel J: Chord diagram of cell-cell interactions
#
# Note: Panels A (qPCR), B (Western blot), and D (fatty acid quantification) 
#       are generated from separate experiments and not included in this script
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
#    - data/LIANAout.rds (pre-computed LIANA results)
#    - AdlungLab_HelperFunctions.R
#
# ==============================================================================

# ==============================================================================
# SETUP: Load required packages
# ==============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(progeny)
library(circlize)
library(scales)

# Load custom theme functions
# This file should contain theme_adlunglab() function
source("AdlungLab_HelperFunctions.R")

# Set working directory to the repository root
# setwd("path/to/ChREBP/")

# Define figure dimensions for publication
fontlist <- list("max" = 10, "min" = 7)
widthlist <- list("min" = 30, "singlecolumn" = 90, "oneandahalf" = 140, "fullwidth" = 190)

# ==============================================================================
# SECTION 1: LOAD DATA
# ==============================================================================
#
# Load the main Seurat object and supporting files
#
# ==============================================================================

# --- Load main Seurat object ---
# Download from: doi: 10.25592/uhhfdm.18248
snuc <- readRDS("ChREBP_SeuratObject_final.rds")

# Join layers (required for Seurat v5)
seurat_obj <- JoinLayers(snuc)

# Load color assignment for consistent plotting
color_assignment <- readRDS("data/color_assignment.rds")

# Define cluster name mappings for publication-ready labels
cluster_names <- c(
  "basal brown adipocytes", 
  "OXPHOS-high brown adipocytes", 
  "stromal cells type 1", 
  "white-like adipocytes", 
  "fast muscle fibers", 
  "capillary endothelial cells", 
  "arterial endothelial cells", 
  "lipogenic brown adipocytes", 
  "endothelial cell-derived brown adipocytes", 
  "stromal cell-derived brown adipocytes", 
  "MonDC", 
  "smooth muscle cells", 
  "stromal cells type 2", 
  "venous endothelial cells", 
  "pericytes", 
  "contractile brown adipocytes", 
  "stromal cells type 3", 
  "satellite glia cells", 
  "slow muscle fibers", 
  "lymphocytes", 
  "schwann cells"
)

# Define adipocyte subtypes for filtering
adipocyte_celltypes <- c(
  "1 Brown adipocyte (Ucp1 low)",
  "2 Brown adipocyte (Ucp1 high)",
  "4 White adipocyte",
  "8 Brown adipocyte (lipogenic)",
  "9 Brown adipocyte (from endothelial cell)",
  "10 Brown adipocyte (from stromal cell)",
  "16 Brown adipocyte (contractile)"
)

# Publication-ready adipocyte names
adipocyte_names <- c(
  "basal brown adipocytes",
  "OXPHOS-high brown adipocytes",
  "white-like adipocytes",
  "lipogenic brown adipocytes",
  "endothelial cell-derived brown adipocytes",
  "stromal cell-derived brown adipocytes",
  "contractile brown adipocytes"
)

# ==============================================================================
# SECTION 2: FIGURE 4C - AGGREGATED EXPRESSION OF DNL GENES
# ==============================================================================
#
# Calculate fold change of aggregated expression of Mlxipl, Acaca, and Fasn
# across experimental conditions (genotype and treatment)
#
# ==============================================================================

# --- Subset to adipocytes only ---
Idents(seurat_obj) <- "celltype"

adipo.behrens <- subset(seurat_obj, idents = adipocyte_celltypes)

message(paste0("Adipocyte subset: ", ncol(adipo.behrens), " nuclei, ", 
               nrow(adipo.behrens), " features"))

# --- Define genes of interest ---
genes <- c("Mlxipl", "Acaca", "Fasn")

# --- Aggregate expression by sample ---
# Sum counts across all adipocytes within each sample
agg_list <- AggregateExpression(
  adipo.behrens, 
  features = genes, 
  group.by = "sample", 
  assays = "RNA", 
  slot = "counts"
)

# Extract RNA assay aggregated data
agg_expr <- agg_list$RNA

# Convert to data frame with samples as rows
df_agg <- as.data.frame(t(agg_expr))
df_agg$sample <- rownames(df_agg)

# --- Convert to long format for plotting ---
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
    # Reference: Wild-Type at Room Temperature
    ref_expression = agg_expression[sample == "g01-Wild-Type-Room Temp."],
    fold_change = agg_expression / ref_expression
  ) %>%
  ungroup()

# --- Generate bar plot ---
fig4c <- ggplot(df_long, aes(x = gene, group = sample, y = fold_change, 
                              color = sample, fill = sample)) +
  geom_col(position = "dodge") +
  theme_adlunglab(base_size = 7) +
  labs(title = NULL,
       x = NULL,
       y = "Fold expression",
       fill = "Sample",
       color = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(fig4c)

# ==============================================================================
# SECTION 3: FIGURE 4E - ADIPOCYTE CLUSTER FREQUENCIES
# ==============================================================================
#
# Relative frequencies of adipocyte clusters by genotype across treatments
#
# ==============================================================================

# --- Prepare metadata for plotting ---
cplot <- as.data.frame(adipo.behrens@meta.data)

# Set factor levels for Treatment
cplot$Treatment <- factor(
  cplot$Treatment, 
  levels = c("Room Temp.", "Acute Cold", "Chronic Cold"),
  labels = c("RT", "AC", "CC")
)

# Set factor levels for Genotype
cplot$Genotype <- factor(
  cplot$Genotype, 
  levels = c("Wild-Type", "Knock-Out"),
  labels = c("Cre-", "Cre+")
)

# --- Generate stacked bar plot ---
fig4e <- ggplot(cplot) +
  geom_bar(aes(x = Genotype, fill = celltype), 
           position = "fill", 
           color = "#000000",
           size = 0.1) +
  theme_adlunglab(base_size = 7) +
  facet_wrap(~Treatment) +
  scale_fill_manual(name = NULL, 
                    values = color_assignment, 
                    labels = adipocyte_names) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = NULL, 
       x = NULL, 
       y = "Relative frequency") +
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "right") +
  guides(fill = guide_legend(override.aes = list(size = 5), ncol = 2))

print(fig4e)

# --- Save Panel E ---
fw <- widthlist$singlecolumn
ggsave("figures/Figure4E_FreqBand_AllAdi_AllConditions.pdf",
       plot = fig4e,
       width = fw, 
       height = fw * 3/4,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

# ==============================================================================
# SECTION 4: FIGURE 4F - PROGENy PATHWAY ANALYSIS
# ==============================================================================
#
# Pathway activity inference using PROGENy across adipocyte subtypes
#
# ==============================================================================

# --- Update celltype labels to publication names ---
seurat_obj_labeled <- seurat_obj

seurat_obj_labeled$celltype <- factor(
  seurat_obj_labeled$celltype, 
  levels = sort(unique(seurat_obj_labeled$celltype)),
  labels = cluster_names
)

# --- Subset to adipocytes ---
Idents(seurat_obj_labeled) <- "celltype"

seurat_adipocytes <- subset(seurat_obj_labeled, idents = c(
  "basal brown adipocytes",
  "OXPHOS-high brown adipocytes",
  "white-like adipocytes",
  "lipogenic brown adipocytes",
  "endothelial cell-derived brown adipocytes",
  "stromal cell-derived brown adipocytes",
  "contractile brown adipocytes"
))

message(paste0("Adipocyte subset for PROGENy: ", ncol(seurat_adipocytes), " nuclei"))

# --- Compute PROGENy activity scores ---
# Uses top 500 most responsive genes per pathway
seurat_adipocytes <- progeny(
  seurat_adipocytes, 
  scale = FALSE, 
  organism = "Mouse", 
  top = 500, 
  perm = 1, 
  return_assay = TRUE
)

# --- Scale pathway activity scores ---
seurat_adipocytes <- Seurat::ScaleData(seurat_adipocytes, assay = "progeny")

# --- Extract PROGENy scores to data frame ---
progeny_scores_df <- as.data.frame(
  t(GetAssayData(seurat_adipocytes, slot = "scale.data", assay = "progeny"))
) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

# --- Match with cell clusters ---
CellsClusters <- data.frame(
  Cell = names(Idents(seurat_adipocytes)), 
  CellType = as.character(Idents(seurat_adipocytes)),
  stringsAsFactors = FALSE
)

progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters, by = "Cell")

# --- Summarize PROGENy scores by cell population ---
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity), .groups = "drop")

# --- Prepare matrix for Z-score normalization ---
progeny_matrix <- summarized_progeny_scores %>%
  dplyr::select(CellType, Pathway, avg) %>%
  spread(CellType, avg, fill = 0) %>%
  column_to_rownames("Pathway")

# --- Perform Z-score normalization per pathway (row-wise) ---
zscore_matrix <- t(scale(t(progeny_matrix)))

# --- Convert back to tidy format for plotting ---
zscore_data <- as.data.frame(zscore_matrix) %>%
  rownames_to_column("Pathway") %>%
  gather(CellType, Zscore, -Pathway)

# --- Apply k-means clustering for ordering ---
set.seed(42)  # For reproducibility
celltype_clusters <- kmeans(t(zscore_matrix), centers = 4)$cluster
pathway_clusters <- kmeans(zscore_matrix, centers = 4)$cluster

# Order by clustering
ordered_celltypes <- names(sort(celltype_clusters))
ordered_pathways <- names(sort(pathway_clusters))

# Update factor levels
zscore_data$CellType <- factor(zscore_data$CellType, levels = ordered_celltypes)
zscore_data$Pathway <- factor(zscore_data$Pathway, levels = ordered_pathways)

# --- Generate heatmap ---
fig4f <- ggplot(zscore_data, aes(x = CellType, y = Pathway, fill = Zscore)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#004488", high = "#BB5566", mid = "white", midpoint = 0) +
  theme_adlunglab(base_size = 7) +
  theme(axis.text.x = element_text(angle = 33, hjust = 1)) +
  labs(title = NULL,
       x = NULL, 
       y = "Pathway", 
       fill = "PROGENy\nZ-score")

print(fig4f)

# ==============================================================================
# SECTION 5: FIGURE 4G - WNT PROGENY SCORES IN LIPOGENIC ADIPOCYTES
# ==============================================================================
#
# Compare WNT pathway activity in lipogenic adipocytes between genotypes
#
# ==============================================================================

# --- Add metadata to PROGENy scores ---
metadata <- seurat_adipocytes@meta.data %>%
  rownames_to_column("Cell")

progeny_scores_df <- progeny_scores_df %>%
  inner_join(metadata, by = "Cell")

# --- Filter for WNT pathway in lipogenic adipocytes ---
wnt_data <- progeny_scores_df %>%
  filter(Pathway == "WNT", CellType == "lipogenic brown adipocytes")

# --- Count cells per genotype ---
cell_counts <- wnt_data %>%
  group_by(Genotype) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Label = paste0("n = ", Count))

# --- Update factor levels for plotting ---
wnt_data$Genotype <- factor(
  wnt_data$Genotype, 
  levels = c("Wild-Type", "Knock-Out"),
  labels = c("Cre-", "Cre+")
)

cell_counts$Genotype <- factor(
  cell_counts$Genotype, 
  levels = c("Wild-Type", "Knock-Out"),
  labels = c("Cre-", "Cre+")
)

# --- Statistical testing ---
# Wilcoxon rank-sum test
stat_test <- wilcox.test(Activity ~ Genotype, data = wnt_data)
message(paste0("Wilcoxon test p-value: ", signif(stat_test$p.value, 3)))

# Linear model controlling for treatment
stat_test2 <- lm(Activity ~ Genotype + Treatment, data = wnt_data)
message("Linear model summary:")
print(summary(stat_test2))

# --- Generate boxplot ---
library(ggsignif)
library(ggbeeswarm)

fig4g <- ggplot(wnt_data, aes(x = Genotype, y = Activity)) +
  geom_boxplot(outlier.shape = NA, fill = "#c79469", size = 0.25) +
  theme_adlunglab(base_size = 7) +
  xlab(NULL) +
  ylab("WNT PROGENy score") +
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()) +
  # Add significance annotation
  geom_signif(
    comparisons = list(c("Cre-", "Cre+")),
    annotations = paste("p =", signif(stat_test$p.value, digits = 2)),
    y_position = max(wnt_data$Activity) * 1.1,
    tip_length = 0.02,
    size = 0.25,
    textsize = 1.5
  ) +
  # Add sample size labels below x-axis
  geom_text(
    data = cell_counts, 
    aes(x = Genotype, y = min(wnt_data$Activity) - 0.1, label = Label),
    size = 2, 
    vjust = 1, 
    color = "black"
  )

print(fig4g)

# ==============================================================================
# SECTION 6: FIGURE 4H - WNT LIGAND EXPRESSION ACROSS TREATMENTS
# ==============================================================================
#
# Aggregated expression of Wnt ligands in Wild-Type mice across treatments
#
# ==============================================================================

# --- Define Wnt ligands of interest ---
wnt_ligands <- c("Wnt2", "Wnt4", "Wnt5a", "Wnt5b", "Wnt6", "Wnt8b", 
                 "Wnt9a", "Wnt9b", "Wnt10b", "Wnt11", "Wnt16")

# --- Extract expression data for Wnt ligands ---
wnt_expr <- Seurat::FetchData(seurat_obj, vars = c(wnt_ligands, "Genotype", "Treatment"))

# --- Filter for Wild-Type and calculate aggregated expression ---
wnt_expr_aggregated <- wnt_expr %>%
  filter(Genotype == "Wild-Type") %>%
  pivot_longer(
    cols = all_of(wnt_ligands), 
    names_to = "Ligand", 
    values_to = "Expression"
  ) %>%
  group_by(Ligand, Treatment) %>%
  summarise(aggregated_expression = mean(Expression, na.rm = TRUE), .groups = "drop")

# --- Set factor levels ---
wnt_expr_aggregated$Ligand <- factor(wnt_expr_aggregated$Ligand, levels = wnt_ligands)
wnt_expr_aggregated$Treatment <- factor(
  wnt_expr_aggregated$Treatment, 
  levels = c("Room Temp.", "Acute Cold", "Chronic Cold"),
  labels = c("RT", "AC", "CC")
)

# --- Generate stacked bar plot ---
fig4h <- ggplot(wnt_expr_aggregated, 
                aes(x = Treatment, y = aggregated_expression, fill = Ligand)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.25) +
  theme_adlunglab(base_size = 7) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(NULL) +
  ylab("Aggregated expression (norm.)") +
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.2, "cm")) +
  scale_fill_brewer(palette = "Paired")

print(fig4h)

# --- Combine Panels F, G, H ---
fig4fgh <- ggarrange(fig4f, fig4g, fig4h, 
                     nrow = 1, 
                     labels = c("F", "G", "H"), 
                     widths = c(1/3, 1/3, 1/3))

print(fig4fgh)

# --- Save combined figure ---
fw <- widthlist$fullwidth
ggsave("figures/Figure4F-H_PROGENy.pdf",
       plot = fig4fgh,
       width = fw, 
       height = fw * 9/16/1.5,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

# ==============================================================================
# OPTIONAL SECTION 7: LIANA CELL-CELL INTERACTION ANALYSIS
# ==============================================================================
#
# This section performs LIANA cell-cell interaction analysis.
# If using pre-computed results, skip to SECTION 8.
#
# Note: LIANA analysis can be computationally intensive.
# Pre-computed results are provided in data/LIANAout.rds
#
# ==============================================================================

# --- Uncomment to run LIANA analysis from scratch ---
# 
# library(liana)
# library(magrittr)
# 
# # Update celltype labels
# seurat_obj_labeled <- seurat_obj
# seurat_obj_labeled$celltype <- factor(
#   seurat_obj_labeled$celltype, 
#   levels = sort(unique(seurat_obj_labeled$celltype)),
#   labels = cluster_names
# )
# 
# Idents(seurat_obj_labeled) <- "celltype"
# 
# # Convert LIANA's Consensus resource to murine symbols
# op_resource <- select_resource("Consensus")[[1]]
# 
# # Generate orthologous resource for mouse (NCBI taxonomy ID: 10090)
# ortholog_resource <- generate_homologs(
#   op_resource = op_resource,
#   target_organism = 10090
# )
# 
# # Run LIANA with orthologous resource
# liana_res <- liana_wrap(
#   seurat_obj_labeled,
#   resource = 'custom',
#   external_resource = ortholog_resource,
#   method = c('sca', 'natmi')  # Run with SCA and NATMI methods
# )
# 
# # Aggregate results
# chrebp_agg <- liana_res %>%
#   liana_aggregate()
# 
# # Save results
# saveRDS(chrebp_agg, "data/LIANAout.rds")

# ==============================================================================
# SECTION 8: FIGURE 4I - LIANA WNT INTERACTION DOT PLOT
# ==============================================================================
#
# Visualize Wnt ligand-receptor interactions with lipogenic brown adipocytes
# as targets using pre-computed LIANA results
#
# ==============================================================================

library(liana)

# --- Load pre-computed LIANA results ---
chrebp_agg <- readRDS("data/LIANAout.rds")

# --- Filter for Wnt ligand interactions ---
# Target: lipogenic brown adipocytes
fig4i <- chrebp_agg %>%
  dplyr::filter(ligand.complex %in% wnt_ligands) %>%
  liana_dotplot(target_groups = "lipogenic brown adipocytes") +
  theme_adlunglab(base_size = 7) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(size = 0.25, colour = "grey")) +
  scale_color_gradient(low = "#F0E634", high = "#220589") +
  labs(x = "lipogenic brown adipocytes\nTarget", 
       y = NULL)

print(fig4i)

# --- Save Panel I ---
fw <- widthlist$fullwidth
ggsave("figures/Figure4I_LR_analysis_LIANA.pdf",
       plot = fig4i,
       width = fw/2, 
       height = fw/2,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

# ==============================================================================
# SECTION 9: FIGURE 4J - CHORD DIAGRAM OF CELL-CELL INTERACTIONS
# ==============================================================================
#
# Chord diagram showing interactions between lipogenic adipocytes,
# stroma-derived adipocytes, and muscle fibers
#
# ==============================================================================

# --- Define custom colors for cell types ---
custom_colors <- c(
  "slow muscle fibers" = "#af497e",
  "fast muscle fibers" = "#bf3335",
  "stromal cell-derived brown adipocytes" = "#dda8aa",
  "lipogenic brown adipocytes" = "#c79469"
)

# --- Prepare interaction data ---
# Filter for interactions between specified cell types
cell_types_of_interest <- c(
  "slow muscle fibers", 
  "fast muscle fibers",
  "stromal cell-derived brown adipocytes", 
  "lipogenic brown adipocytes"
)

freqs <- chrebp_agg %>%
  dplyr::filter(
    source %in% cell_types_of_interest,
    target %in% cell_types_of_interest
  ) %>%
  dplyr::select(source, target, aggregate_rank) %>%
  # Convert rank to interaction score (higher score = stronger interaction)
  dplyr::mutate(
    interaction_score = scales::rescale(1 / aggregate_rank, to = c(1, 10))
  ) %>%
  dplyr::select(source, target, interaction_score)

message(paste0("Number of interactions: ", nrow(freqs)))

# --- Generate chord diagram ---
# Save to PDF
pdf("figures/Figure4J_Chord.pdf", width = 75/25.4, height = 75/25.4)

# Clear previous circos plots
circos.clear()

# Set gap between sectors
circos.par(gap.degree = 5)

# Draw chord diagram
chordDiagram(
  freqs,
  directional = 1, 
  direction.type = c("diffHeight", "arrows"), 
  annotationTrack = "grid",
  preAllocateTracks = 1,
  transparency = 0.1,
  grid.col = custom_colors,
  link.arr.col = NA
)

dev.off()

# Also display in R session
circos.clear()
circos.par(gap.degree = 5)
chordDiagram(
  freqs,
  directional = 1, 
  direction.type = c("diffHeight", "arrows"), 
  annotationTrack = "grid",
  preAllocateTracks = 1,
  transparency = 0.1,
  grid.col = custom_colors,
  link.arr.col = NA
)

message("Chord diagram saved to figures/Figure4J_Chord.pdf")

# ==============================================================================
# SECTION 10: SAVE ALL PANELS
# ==============================================================================
#
# Save individual panels and combined figures
#
# ==============================================================================

# --- Save Panel C ---
ggsave("figures/Figure4C_AggregatedExpression_DNL.pdf",
       plot = fig4c,
       width = widthlist$singlecolumn, 
       height = widthlist$singlecolumn * 3/4,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

# --- Save Panel F ---
ggsave("figures/Figure4F_PROGENy_Heatmap.pdf",
       plot = fig4f,
       width = widthlist$singlecolumn, 
       height = widthlist$singlecolumn * 3/4,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

# --- Save Panel G ---
ggsave("figures/Figure4G_WNT_Boxplot.pdf",
       plot = fig4g,
       width = widthlist$min * 2, 
       height = widthlist$min * 2.5,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

# --- Save Panel H ---
ggsave("figures/Figure4H_WntLigands_StackedBar.pdf",
       plot = fig4h,
       width = widthlist$singlecolumn, 
       height = widthlist$singlecolumn * 3/4,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

# ==============================================================================
# SESSION INFO
# ==============================================================================

sessionInfo()

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
