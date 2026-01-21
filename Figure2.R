################################################################################
#
# Figure 2: Ttc25 is a specific marker of lipogenic adipocytes strongly 
#           associated with DNL and thermogenic metabolism
#
# Behrens et al. (2025) - Molecular Metabolism
# "Single-nucleus mRNA-sequencing reveals dynamics of lipogenic and thermogenic 
#  adipocyte populations in murine brown adipose tissue in response to cold exposure"
#
# This script generates the snRNA-seq panels of Figure 2:
#   - Panel A: Scatter plot comparing lipogenic adipocytes (this study vs Lundgren et al.)
#   - Panel B: Volcano plot of differential expression (lipogenic vs basal adipocytes)
#   - Panel C: UMAP showing Ttc25 expression
#   - Panel D: UMAP showing Fasn expression
#   - Panel G: Waterfall plot of genes correlating with Ttc25
#   - Panel H: GO term analysis of Ttc25-correlated genes
#
# Note: Panels E (qPCR ChREBP KO) and F (qPCR ChREBPβ OE) are generated from 
#       separate experiments and not included in this script
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
# 2. External data from Lundgren et al. (2023) Nature Metabolism:
#    - Download from GEO: GSE218711
#    - Required folders: GSE218711_RAW/TN/ and GSE218711_RAW/CYC/
#    - Each folder contains: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
#
# 3. Supporting files (provided in the repository):
#    - data/color_assignment.rds
#    - AdlungLab_HelperFunctions.R
#
# 4. Pre-computed intermediate files (optional, provided in repository):
#    - data/mergedLundgrenBehrensLipogenic.rds  (Panel A data)
#    - data/Lipo_vs_basal_WTRT.csv              (Panel B data)
#
# ==============================================================================

# ==============================================================================
# SETUP: Load required packages
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# Load custom theme functions
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
snuc <- JoinLayers(snuc)

# Load color assignment for consistent plotting
color_assignment <- readRDS("data/color_assignment.rds")

# ==============================================================================
# OPTIONAL SECTION 2: PROCESS LUNDGREN ET AL. DATA FOR COMPARISON
# ==============================================================================
#
# This section processes publicly available snRNA-seq data from:
# Lundgren et al. (2023) Nature Metabolism
# "A subpopulation of lipogenic brown adipocytes drives thermogenic memory"
# GEO: GSE218711
#
# If using pre-computed data, skip to SECTION 3
#
# ==============================================================================

# --- Load Lundgren et al. data from GEO ---
# Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218711
# Place in: OtherData/Lundgren_etal_NatMetab_2023/GSE218711_RAW/

# Thermoneutrality (TN) sample
data <- Read10X(data.dir = "OtherData/Lundgren_etal_NatMetab_2023/GSE218711_RAW/TN/")
TN <- CreateSeuratObject(counts = data, 
                         project = "Lundgren",
                         min.cells = 3)
TN$Treatment <- "Thermoneutrality"

# Cold cycling (CYC) sample  
data <- Read10X(data.dir = "OtherData/Lundgren_etal_NatMetab_2023/GSE218711_RAW/CYC/")
CYC <- CreateSeuratObject(counts = data, 
                          project = "Lundgren",
                          min.cells = 3)
CYC$Treatment <- "Cold"

rm(data)

# --- Quality control: TN sample ---
TN[["percent.mt"]] <- PercentageFeatureSet(TN, pattern = "^mt-")
# VlnPlot(TN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
TN <- subset(TN, subset = percent.mt < 25 & nCount_RNA < 5000 & nFeature_RNA < 2000)
TN <- NormalizeData(TN)

# --- Quality control: CYC sample ---
CYC[["percent.mt"]] <- PercentageFeatureSet(CYC, pattern = "^mt-")
# VlnPlot(CYC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CYC <- subset(CYC, subset = percent.mt < 50 & nCount_RNA < 2500 & nFeature_RNA < 1500)
CYC <- NormalizeData(CYC)

# --- Merge samples ---
lundgren.combined <- merge(TN, y = c(CYC),
                           add.cell.ids = c("TN", "CYC"), 
                           project = "lundgren", 
                           merge.data = TRUE)

message(paste0("Lundgren merged: ", ncol(lundgren.combined), " nuclei"))

rm(TN, CYC)

# --- Integration with Harmony ---
library(harmony)

lundgren.combined$sample <- paste0(lundgren.combined$Treatment)

lundgren.combined.integrated <- NormalizeData(lundgren.combined) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE)

lundgren.combined.integrated <- RunHarmony(lundgren.combined.integrated, 
                                            group.by.vars = "sample")

lundgren.combined.integrated <- RunUMAP(lundgren.combined.integrated, 
                                         reduction = "harmony", 
                                         dims = 1:30)

lundgren.combined.integrated <- FindNeighbors(lundgren.combined.integrated, 
                                               reduction = "harmony", 
                                               dims = 1:30) %>% 
  FindClusters(resolution = 0.5)

# DimPlot(lundgren.combined.integrated, group.by = c("sample", "seurat_clusters"), ncol = 2)

# --- Identify adipocyte clusters and subset to lipogenic adipocytes ---
lundgren <- lundgren.combined.integrated
lundgren2 <- JoinLayers(lundgren)

# Identify adipocyte clusters based on marker gene expression
# VlnPlot(lundgren, features = c("Acaca", "Fasn", "Acly", "Scd1", "Ppara", "Cidea", "Lipe"), sort = TRUE)
# DimPlot(lundgren, label = TRUE)

# Subset to adipocyte clusters (identified as 0, 1, 2, 4, 10)
adipocytes <- subset(lundgren2, idents = c(0, 1, 2, 4, 10))

# Recluster adipocytes
adipocytes <- NormalizeData(adipocytes)
adipocytes <- FindVariableFeatures(adipocytes, selection.method = "vst", nfeatures = 2000)
adipocytes <- ScaleData(adipocytes)
adipocytes <- RunPCA(adipocytes, npcs = 30)
adipocytes <- FindNeighbors(adipocytes, dims = 1:30)
adipocytes <- FindClusters(adipocytes, resolution = 0.5)

# Identify lipogenic cluster (cluster 2 based on DNL gene expression)
# VlnPlot(adipocytes, features = c("Acaca", "Fasn", "Acly", "Scd1", "Ppara", "Cidea", "Lipe"), sort = TRUE)

lipo.lundgren <- subset(adipocytes, idents = c(2))
Idents(lipo.lundgren) <- "Treatment"
lipocold.lundgren <- subset(lipo.lundgren, idents = c("Cold"))

# Calculate average expression for Lundgren lipogenic adipocytes (cold)
pb.lundgren <- AverageExpression(lipocold.lundgren)
pb.exp.lundgren <- as.data.frame(pb.lundgren$RNA)

# --- Extract lipogenic adipocytes from our data (Behrens et al.) ---
Idents(snuc) <- "sample"
wtrt <- subset(snuc, idents = "03_Wild-Type_Chronic Cold")
Idents(wtrt) <- "celltype"

lipocold.behrens <- subset(wtrt, idents = "8 Brown adipocyte (lipogenic)")

# Calculate average expression for Behrens lipogenic adipocytes (cold)
pb.behrens <- AverageExpression(lipocold.behrens)
pb.exp.behrens <- as.data.frame(pb.behrens$RNA)

# --- Merge expression data from both studies ---
merged <- merge(pb.exp.lundgren, pb.exp.behrens, by = "row.names", all = TRUE)
rownames(merged) <- merged$Row.names
colnames(merged) <- c("gene_name", "Lundgren", "Behrens")

# Define signature genes to highlight
signaturegenes <- c("Fasn", "Mlxipl", "Acly", "Acaca", "Ttc25")

# Remove NAs and add signature annotation
merged.plot <- na.omit(merged)
merged.plot$signature <- ifelse(merged.plot$gene_name %in% signaturegenes, 
                                 "signature", "none")

# --- Save intermediate file ---
saveRDS(merged.plot, "data/mergedLundgrenBehrensLipogenic.rds")

# Clean up
rm(lundgren.combined, lundgren.combined.integrated, lundgren, lundgren2, 
   adipocytes, lipo.lundgren, lipocold.lundgren, lipocold.behrens,
   pb.lundgren, pb.behrens, pb.exp.lundgren, pb.exp.behrens, merged)

# ==============================================================================
# SECTION 3: FIGURE 2 PANEL A - SCATTER PLOT COMPARISON
# ==============================================================================
#
# Compare gene expression in lipogenic brown adipocytes between this study
# and Lundgren et al. (2023), both under cold exposure conditions
#
# ==============================================================================

# --- Load pre-computed data (or use data from SECTION 2) ---
merged.plot <- readRDS("data/mergedLundgrenBehrensLipogenic.rds")

# --- Calculate correlation statistics ---
# Overall correlation
ct <- cor.test(merged.plot$Lundgren, merged.plot$Behrens, 
               exact = FALSE, method = "pearson")
message(paste0("Overall Pearson correlation: r = ", round(ct$estimate, 2)))

# Correlation for signature genes only
test <- merged.plot[which(merged.plot$signature %in% c("signature")), ]
ct_sig <- cor.test(test$Lundgren, test$Behrens)
message(paste0("Signature genes correlation: r = ", round(ct_sig$estimate, 2)))

# --- Generate scatter plot ---
p2a <- ggplot() + 
  labs(title = bquote("Lipogenic brown adipocytes, cold," ~ rho ~ "=" ~ .(round(ct$estimate, 2))),
       y = "log10(norm. mRNA), this work",
       x = "log10(norm. mRNA), Lundgren") +
  # Diagonal reference line (y = x)
  geom_abline(slope = 1, intercept = 0, colour = "#000000", size = 0.5) +
  # Non-signature genes (grey)
  geom_point(data = merged.plot[-which(merged.plot$signature %in% c("signature")), ], 
             aes(x = log10(Lundgren), y = log10(Behrens)), 
             size = 0.66, color = "grey") +
  # Signature genes (highlighted)
  geom_point(data = merged.plot[which(merged.plot$signature %in% c("signature")), ], 
             aes(x = log10(Lundgren), y = log10(Behrens)), 
             size = 0.66, color = "#c79469") +
  # Gene labels for signature genes
  geom_text_repel(data = merged.plot[which(merged.plot$signature %in% c("signature")), ], 
                  aes(x = log10(Lundgren), y = log10(Behrens), label = gene_name), 
                  color = "#000000", 
                  max.overlaps = 23,
                  segment.size = 0.1, 
                  size = 2.5,
                  min.segment.length = 0, 
                  direction = "both") +
  scale_colour_manual(values = rev(c("#000000", "#FFFFFF")), guide = "none") + 
  scale_fill_manual(values = rev(c("#220589", "grey")), guide = "none") + 
  theme_adlunglab(base_size = 7) +
  theme(plot.title = element_text(hjust = 0.5, size = 7), 
        legend.position = "none")

print(p2a)

# ==============================================================================
# OPTIONAL SECTION 4: DIFFERENTIAL EXPRESSION ANALYSIS
# ==============================================================================
#
# Calculate differential expression between lipogenic and basal brown adipocytes
# in Wild-Type mice at room temperature
#
# If using pre-computed data, skip to SECTION 5
#
# ==============================================================================

# --- Subset to WT Room Temperature ---
Idents(snuc) <- "id"
snuc2 <- subset(snuc, idents = "01")
Idents(snuc2) <- "celltype"

# --- Find differentially expressed genes ---
# Compare lipogenic (cluster 8) vs basal (cluster 1) brown adipocytes
diff_markers <- FindMarkers(snuc2, 
                            ident.1 = "8 Brown adipocyte (lipogenic)", 
                            ident.2 = "1 Brown adipocyte (Ucp1 low)", 
                            only.pos = FALSE, 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25)
diff_markers$gene <- rownames(diff_markers)

# --- Save results ---
write.csv(diff_markers, "data/Lipo_vs_basal_WTRT.csv", row.names = FALSE)

# ==============================================================================
# SECTION 5: FIGURE 2 PANEL B - VOLCANO PLOT
# ==============================================================================
#
# Volcano plot showing differential gene expression between lipogenic and
# basal brown adipocytes (Wild-Type, room temperature)
#
# ==============================================================================

# --- Load differential expression results ---
vulc <- read.csv("data/Lipo_vs_basal_WTRT.csv")

# --- Calculate significance thresholds ---
# Use 1.96 standard deviations of log2FC as threshold
fc_threshold <- 1.96 * sd(vulc$avg_log2FC)

# Determine highlighted genes
vulc$highlight <- vulc$p_val_adj < 0.05 & abs(vulc$avg_log2FC) > fc_threshold

# Assign colors based on direction and significance
vulc$color <- ifelse(vulc$highlight & vulc$avg_log2FC > 0, "#c79469",
                     ifelse(vulc$highlight & vulc$avg_log2FC < 0, "#af6729", "grey"))

# Label key DNL genes
vulc$label <- ifelse(vulc$gene %in% c("Acly", "Ttc25", "Acaca", "Fasn"), 
                     vulc$gene, NA)

# --- Generate volcano plot ---
library(grid)  # For arrow() function

gvl <- ggplot(vulc, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_vline(xintercept = 0, colour = "#000000", size = 0.5) +
  geom_point(aes(color = color), size = 0.66) +
  scale_color_identity() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) +
  geom_text_repel(aes(label = label), 
                  segment.size = 0.1, 
                  size = 2.5,
                  min.segment.length = 0,
                  direction = "both") +
  theme_adlunglab(base_size = 7) +
  labs(title = NULL,
       x = "log2(fold-change)",
       y = "-log10(p adj.)",
       color = NULL) +
  # Arrow annotation: lipogenic direction
  annotate("segment",
           x = 0, xend = 1.25,
           y = 140, yend = 140,
           arrow = arrow(length = unit(2, "mm"), type = "closed"),
           size = 0.5) +
  annotate("text",
           x = 0.22, y = 145,
           label = "lipogenic",
           hjust = 0,
           size = 2.5) +
  # Arrow annotation: basal direction
  annotate("segment",
           x = 0, xend = -1.25,
           y = 140, yend = 140,
           arrow = arrow(length = unit(2, "mm"), type = "closed"),
           size = 0.5) +
  annotate("text",
           x = -0.22, y = 145,
           label = "basal",
           hjust = 1,
           size = 2.5)

print(gvl)

# ==============================================================================
# SECTION 6: FIGURE 2 PANELS C & D - UMAP FEATURE PLOTS
# ==============================================================================
#
# UMAP plots showing expression of Ttc25 (Panel C) and Fasn (Panel D)
# in Wild-Type mice at room temperature
#
# ==============================================================================

# --- Subset to WT Room Temperature ---
Idents(snuc) <- "sample"
wtrt <- subset(snuc, idents = "01_Wild-Type_Room Temp.")
Idents(wtrt) <- "celltype"

# --- Prepare data for plotting ---
cplot <- as.data.frame(wtrt@meta.data)
cplot$x <- wtrt[["umap"]]@cell.embeddings[, 1]
cplot$y <- wtrt[["umap"]]@cell.embeddings[, 2]

# Extract gene expression values
cplot$Ttc25 <- as.numeric(FetchData(wtrt, vars = "Ttc25", slot = "data")$Ttc25)
cplot$Fasn <- as.numeric(FetchData(wtrt, vars = "Fasn", slot = "data")$Fasn)

# --- Panel C: Ttc25 expression UMAP ---
c2 <- ggplot(cplot, aes(x = x, y = y, color = Ttc25)) + 
  geom_point(shape = 16, size = 0.25) + 
  scale_color_gradient(low = "grey89", high = "#220589", name = "Ttc25") +
  labs(x = NULL, y = NULL) +
  theme_adlunglab(base_size = 7) +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  guides(color = guide_colorbar(barwidth = 0.3))

print(c2)

# --- Panel D: Fasn expression UMAP ---
c1 <- ggplot(cplot, aes(x = x, y = y, color = Fasn)) + 
  geom_point(shape = 16, size = 0.25) + 
  scale_color_gradient(low = "grey89", high = "#220589", name = "Fasn") +
  labs(x = NULL, y = NULL) +
  theme_adlunglab(base_size = 7) +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  guides(color = guide_colorbar(barwidth = 0.3))

print(c1)

# --- Combine UMAP panels ---
fig2cd <- ggarrange(c2, c1, nrow = 2, heights = c(0.5, 0.5), labels = c("C", "D"))

# ==============================================================================
# SECTION 7: FIGURE 2 PANELS G & H - TTC25 CORRELATION ANALYSIS
# ==============================================================================
#
# Identify genes correlated with Ttc25 in lipogenic brown adipocytes
# and perform Gene Ontology enrichment analysis
#
# ==============================================================================

# --- Subset to lipogenic adipocytes in Wild-Type mice ---
Idents(snuc) <- "celltype"
adi.behrens <- subset(snuc, idents = c("8 Brown adipocyte (lipogenic)"))
Idents(adi.behrens) <- "Genotype"
adi.behrens <- subset(adi.behrens, idents = "Wild-Type")

message(paste0("Lipogenic adipocytes (WT): ", ncol(adi.behrens), " cells"))

# --- Calculate Ttc25 correlation with all genes ---
# Extract expression matrix
expr_matrix <- GetAssayData(adi.behrens, layer = "data")

# Get Ttc25 expression vector
ttc25_expr <- expr_matrix["Ttc25", ]

# Calculate Pearson correlation for each gene
gene_correlations <- apply(expr_matrix, 1, function(x) {
  cor(x, ttc25_expr, method = "pearson", use = "all.obs")
})

# Create correlation dataframe
cor_df <- data.frame(
  Gene = names(gene_correlations),
  Correlation = gene_correlations
)

# Sort by correlation and remove NAs and Ttc25 itself
cor_df <- cor_df %>% 
  arrange(desc(Correlation)) %>%
  filter(!is.na(Correlation)) %>%
  filter(Gene != "Ttc25")

# Add labels for key genes
cor_df$text <- ifelse(cor_df$Gene %in% c("Fasn", "Acaca", "Acly", "Mlxipl"), 
                      cor_df$Gene, NA)

# --- Panel G: Waterfall plot of gene correlations ---
wf <- ggplot(cor_df, aes(x = reorder(Gene, Correlation), y = Correlation)) +
  geom_bar(stat = "identity", color = "#c79469", fill = "#c79469", size = 0) +
  geom_text_repel(aes(label = text), 
                  min.segment.length = 0, 
                  segment.size = 0.1, 
                  size = 2.5) + 
  theme_adlunglab(base_size = 7) +
  geom_hline(yintercept = 0, size = 0.5, colour = "#000000") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  coord_flip() +
  labs(title = NULL, 
       x = NULL,
       y = "Pearson's correlation")

print(wf)

# --- Calculate correlation statistics with p-values ---
cor_results <- apply(expr_matrix, 1, function(x) {
  test <- cor.test(x, ttc25_expr, method = "pearson", use = "complete.obs")
  return(c(r = test$estimate, pval = test$p.value))
})

# Convert to dataframe
cor_df_stats <- as.data.frame(t(cor_results))
cor_df_stats$Gene <- rownames(cor_df_stats)

# Adjust p-values for multiple testing (Benjamini-Hochberg)
cor_df_stats$adj_pval <- p.adjust(cor_df_stats$pval, method = "BH")

# Filter for significantly positively correlated genes
sig_genes <- cor_df_stats %>% 
  filter(r.cor > 0 & adj_pval < 0.001) %>% 
  arrange(desc(r.cor))

message(paste0("Significantly correlated genes (r > 0, adj.p < 0.001): ", nrow(sig_genes)))

# --- Panel H: Gene Ontology enrichment analysis ---
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse annotation database
library(enrichplot)

# Get top correlated genes
top_genes <- sig_genes %>% 
  filter(adj_pval < 0.001) %>% 
  pull(Gene)

# Convert gene symbols to Entrez IDs
entrez_ids <- mapIds(org.Mm.eg.db, 
                     keys = top_genes, 
                     column = "ENTREZID", 
                     keytype = "SYMBOL", 
                     multiVals = "first")

# Remove NA values
entrez_ids <- na.omit(entrez_ids)

message(paste0("Genes with Entrez IDs: ", length(entrez_ids)))

# Perform GO enrichment analysis (Biological Process)
go_results <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Mm.eg.db, 
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Generate dot plot
dp <- dotplot(go_results, showCategory = 15) + 
  theme_adlunglab(base_size = 7) +
  geom_segment(aes(y = reorder(Description, Count), 
                   yend = reorder(Description, Count), 
                   x = 0, xend = GeneRatio), 
               color = "grey", linetype = "dashed", size = 0.1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.115))

print(dp)

# ==============================================================================
# SECTION 8: COMBINE AND SAVE FIGURE 2
# ==============================================================================
#
# Arrange all panels and save as PDF/PNG
#
# ==============================================================================

# --- Combine panels A, B, C, D ---
fig2abcd <- ggarrange(p2a, gvl, fig2cd, 
                      nrow = 1, 
                      widths = c(0.35, 0.35, 0.3), 
                      labels = c("A", "B", ""))

print(fig2abcd)

# --- Save panels A-D ---
fw <- widthlist$fullwidth
ggsave("figures/Figure2A-D.pdf",
       plot = fig2abcd,
       width = fw, 
       height = fw * 9/16/1.5, 
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

ggsave("figures/Figure2A-D.png",
       plot = fig2abcd,
       width = fw, 
       height = fw * 9/16/1.5, 
       units = "mm", 
       dpi = 300)

# --- Save waterfall plot (Panel G) ---
fw_small <- widthlist$fullwidth * 0.25
ggsave("figures/Figure2G_WaterfallPlot.pdf",
       plot = wf,
       width = fw_small, 
       height = fw_small * 16/9, 
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

# --- Combine waterfall and GO analysis (Panels G & H) ---
fig2gh <- ggarrange(wf, dp, nrow = 1, widths = c(0.25, 0.75), labels = c("G", "H"))

ggsave("figures/Figure2G-H_Ttc25Correlation.pdf", 
       plot = fig2gh,
       width = 190, 
       height = 90,
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

ggsave("figures/Figure2G-H_Ttc25Correlation.png", 
       plot = fig2gh,
       width = 190, 
       height = 90,
       units = "mm", 
       dpi = 300)

# ==============================================================================
# SESSION INFO
# ==============================================================================

sessionInfo()

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
