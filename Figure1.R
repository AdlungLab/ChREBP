################################################################################
#
# Figure 1: Single-nucleus RNA-seq analysis of murine brown adipose tissue
#
# Behrens et al. (2025) - Molecular Metabolism
# "Single-nucleus mRNA-sequencing reveals dynamics of lipogenic and thermogenic 
#  adipocyte populations in murine brown adipose tissue in response to cold exposure"
#
# This script generates the snRNA-seq panels of Figure 1:
#   - Panel C: UMAP visualization of 21 cell clusters
#   - Panel D: Brown and white adipocyte scores
#   - Panel E: Dot plot of marker gene expression
#
# Note: Panels A (qPCR) and B (Western blot) are generated from separate experiments
#
################################################################################

# ==============================================================================
# REQUIRED INPUT FILES
# ==============================================================================
# 
# 1. Raw 10X Genomics data folders (each containing barcodes.tsv.gz, 
#    features.tsv.gz, matrix.mtx.gz):
#    - AlignedFiltered/KO_AC-SCI7T058_HGVMFDSX7/
#    - AlignedFiltered/KO_CC-SCI7T070_HGW3CDSX7/
#    - AlignedFiltered/KO_RT-SCI7T046_HGVMFDSX7/
#    - AlignedFiltered/WT_RT-SCI7T010_HGWFKDSX7/
#    - AlignedFiltered/WT_AC-SCI7T022_HGVMFDSX7/
#    - AlignedFiltered/WT_CC-SCI7T034_HGVMFDSX7/
#   https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15819 
#
# 2. Supporting files (provided in the repository):
#    - data/doublet_barcodes.csv          (DoubletFinder output)
#    - data/CellTypeAnnotation_ChREBP_final.csv  (Cell type annotations)
#    - data/SupplementalTable1_GeneScores.csv    (Gene lists for scoring)
#
# 3. Alternative: Final Seurat object:
#    - doi: 10.25592/uhhfdm.18248
#
# ==============================================================================

# ==============================================================================
# SETUP: Load required packages
# ==============================================================================

library(Seurat)
library(dplyr)
library(patchwork)
library(harmony)
library(ggplot2)

# Set your working directory to the repository root
# setwd("path/to/ChREBP/")

# ==============================================================================
# OPTIONAL SECTION 1: DATA LOADING
# ==============================================================================
# 
# This section loads raw 10X Genomics data from six samples:
# - 3 wild-type (WT) samples: Room temp, Acute cold, Chronic cold
# - 3 knock-out (KO) samples: Room temp, Acute cold, Chronic cold
#
# If you have the pre-processed Seurat object, skip to SECTION 5
# ==============================================================================

# --- Sample 1: KO Acute Cold ---
data <- Read10X(data.dir = "AlignedFiltered/KO_AC-SCI7T058_HGVMFDSX7/")
KO_AC <- CreateSeuratObject(counts = data, 
                            project = "ChREBP",
                            min.cells = 3)
KO_AC$id <- "05"
KO_AC$Genotype <- "Knock-Out"
KO_AC$Treatment <- "Acute Cold"

# --- Sample 2: KO Chronic Cold ---
data <- Read10X(data.dir = "AlignedFiltered/KO_CC-SCI7T070_HGW3CDSX7/")
KO_CC <- CreateSeuratObject(counts = data, 
                            project = "ChREBP",
                            min.cells = 3)
KO_CC$id <- "06"
KO_CC$Genotype <- "Knock-Out"
KO_CC$Treatment <- "Chronic Cold"

# --- Sample 3: KO Room Temperature ---
data <- Read10X(data.dir = "AlignedFiltered/KO_RT-SCI7T046_HGVMFDSX7/")
KO_RT <- CreateSeuratObject(counts = data, 
                            project = "ChREBP",
                            min.cells = 3)
KO_RT$id <- "04"
KO_RT$Genotype <- "Knock-Out"
KO_RT$Treatment <- "Room Temp."

# --- Sample 4: WT Room Temperature ---
data <- Read10X(data.dir = "AlignedFiltered/WT_RT-SCI7T010_HGWFKDSX7/")
WT_RT <- CreateSeuratObject(counts = data, 
                            project = "ChREBP",
                            min.cells = 3)
WT_RT$id <- "01"
WT_RT$Genotype <- "Wild-Type"
WT_RT$Treatment <- "Room Temp."

# --- Sample 5: WT Acute Cold ---
data <- Read10X(data.dir = "AlignedFiltered/WT_AC-SCI7T022_HGVMFDSX7/")
WT_AC <- CreateSeuratObject(counts = data, 
                            project = "ChREBP",
                            min.cells = 3)
WT_AC$id <- "02"
WT_AC$Genotype <- "Wild-Type"
WT_AC$Treatment <- "Acute Cold"

# --- Sample 6: WT Chronic Cold ---
data <- Read10X(data.dir = "AlignedFiltered/WT_CC-SCI7T034_HGVMFDSX7/")
WT_CC <- CreateSeuratObject(counts = data, 
                            project = "ChREBP",
                            min.cells = 3)
WT_CC$id <- "03"
WT_CC$Genotype <- "Wild-Type"
WT_CC$Treatment <- "Chronic Cold"

# Clean up temporary variable
rm(data)

# ==============================================================================
# OPTIONAL SECTION 2: QUALITY CONTROL AND NORMALIZATION
# ==============================================================================
#
# QC filtering criteria (applied to all samples):
# - Mitochondrial content < 15%
# - Total UMI counts < 25,000
# - Number of detected genes < 5,000
#
# ==============================================================================

# --- QC function to apply consistent filtering ---
qc_and_normalize <- function(seurat_obj, sample_name) {
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  # Visualize QC metrics (optional - uncomment to view)
  # print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  
  # Apply QC filters
  seurat_obj <- subset(seurat_obj, 
                       subset = percent.mt < 15 & 
                                nCount_RNA < 25000 & 
                                nFeature_RNA < 5000)
  
  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj)
  
  message(paste0("QC complete for ", sample_name, ": ", 
                 ncol(seurat_obj), " cells retained"))
  return(seurat_obj)
}

# Apply QC and normalization to each sample
WT_RT <- qc_and_normalize(WT_RT, "WT_RT")
WT_AC <- qc_and_normalize(WT_AC, "WT_AC")
WT_CC <- qc_and_normalize(WT_CC, "WT_CC")
KO_RT <- qc_and_normalize(KO_RT, "KO_RT")
KO_AC <- qc_and_normalize(KO_AC, "KO_AC")
KO_CC <- qc_and_normalize(KO_CC, "KO_CC")

# ==============================================================================
# OPTIONAL SECTION 3: DATA INTEGRATION
# ==============================================================================
#
# Merge all samples and integrate using Harmony to remove batch effects
# while preserving biological variation
#
# ==============================================================================

# --- Merge all samples into a single Seurat object ---
chrebp.combined <- merge(WT_RT, 
                         y = c(WT_AC, WT_CC, KO_RT, KO_AC, KO_CC),
                         add.cell.ids = c("WT_RT", "WT_AC", "WT_CC",
                                          "KO_RT", "KO_AC", "KO_CC"), 
                         project = "ChREBP", 
                         merge.data = TRUE)

# Check merged object
message(paste0("Merged object: ", ncol(chrebp.combined), " nuclei, ",
               nrow(chrebp.combined), " features"))

# Optional: Save intermediate merged object
# saveRDS(chrebp.combined, "Output/merged_ChREBP_redone.rds")

# Clean up individual sample objects to free memory
rm(WT_RT, WT_AC, WT_CC, KO_RT, KO_AC, KO_CC)
gc()

# --- Create sample identifier for Harmony integration ---
chrebp.combined$sample <- paste0(chrebp.combined$id, "_", 
                                  chrebp.combined$Genotype, "_", 
                                  chrebp.combined$Treatment)

# --- Standard Seurat processing pipeline ---
# Normalize, find variable features, scale, and run PCA
chrebp.combined.integrated <- NormalizeData(chrebp.combined) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE)

# --- Harmony integration ---
# Corrects for batch effects based on sample identity
chrebp.combined.integrated <- RunHarmony(chrebp.combined.integrated, 
                                          group.by.vars = "sample")

# --- UMAP and clustering ---
# Using Harmony-corrected dimensions (1:30)
chrebp.combined.integrated <- RunUMAP(chrebp.combined.integrated, 
                                       reduction = "harmony", 
                                       dims = 1:30)

chrebp.combined.integrated <- FindNeighbors(chrebp.combined.integrated, 
                                             reduction = "harmony", 
                                             dims = 1:30) %>% 
  FindClusters()

# Visualize integration results
DimPlot(chrebp.combined.integrated, 
        group.by = c("sample", "seurat_clusters"), 
        ncol = 2)

# Check final object dimensions
message(paste0("Integrated object: ", ncol(chrebp.combined.integrated), " nuclei"))

# Optional: Save integrated object
# saveRDS(chrebp.combined.integrated, "Output/harmony_ChREBP_redone.rds")

# ==============================================================================
# OPTIONAL SECTION 4: DOUBLET REMOVAL AND CELL TYPE ANNOTATION
# ==============================================================================
#
# Remove RBC cluster and annotate doublets using DoubletFinder results
#
# ==============================================================================

# --- Remove RBC cluster ---
# Cluster 0 was identified as red blood cells based on marker gene expression
chrebp <- chrebp.combined.integrated
snuc <- subset(x = chrebp, idents = c("0 RBC"), invert = TRUE)

# Optional: Save object without RBCs
# saveRDS(snuc, "ChREBP_SeuratObject.rds")

# --- Add doublet information ---
# Load doublet calls from DoubletFinder analysis
dbl <- read.csv("data/doublet_barcodes.csv", stringsAsFactors = FALSE)

# Create full barcode including sample prefix
dbl <- dbl %>%
  mutate(FullBarcode = paste0(Sample, "_", Barcode))

# Mark doublets in the Seurat object
cells <- colnames(snuc)
is_dbl <- cells %in% dbl$FullBarcode
snuc$DoubletCall <- ifelse(is_dbl, "Doublet", "Singlet")
snuc$DoubletCall <- factor(snuc$DoubletCall, levels = c("Singlet", "Doublet"))

# View doublet distribution
table(snuc$DoubletCall)

# ==============================================================================
# OPTIONAL SECTION 5: CELL TYPE ANNOTATION
# ==============================================================================
#
# Apply cell type annotations based on marker gene analysis
# Mapping from original cluster IDs to descriptive cell type names
#
# ==============================================================================

# --- Load cell type annotations ---
celltypes <- read.csv("data/CellTypeAnnotation_ChREBP_final.csv")
celltypes$showName <- paste0(celltypes$Number, " ", celltypes$CellType)

# --- Apply cluster names ---
clusterID <- celltypes$showName
names(clusterID) <- levels(snuc)
snuc <- RenameIdents(snuc, clusterID)
snuc$celltype <- Idents(snuc)
rownames(celltypes) <- celltypes$ClusterID

# --- Define color scheme for cell types ---
color_assignment <- setNames(celltypes$celltypecolor, celltypes$showName)

# --- Define mapping to new descriptive names ---
# This mapping converts numbered cluster names to descriptive cell type names
mapping <- c(
  "1 Brown adipocyte (Ucp1 low)"               = "basal brown adipocytes",
  "2 Brown adipocyte (Ucp1 high)"              = "OXPHOS-high brown adipocytes",
  "3 Stromal cell-1"                            = "stromal cells type 1",
  "4 White adipocyte"                           = "white-like adipocytes",
  "5 Muscle fiber (fast)"                       = "fast muscle fibers",
  "6 Endothelial cell (capillary)"              = "capillary endothelial cells",
  "7 Endothelial cell (arterial)"               = "arterial endothelial cells",
  "8 Brown adipocyte (lipogenic)"               = "lipogenic brown adipocytes",
  "9 Brown adipocyte (from endothelial cell)"   = "endothelium-derived adipocytes",
  "10 Brown adipocyte (from stromal cell)"      = "stroma-derived adipocytes",
  "11 MonDC"                                    = "MonDC",
  "12 Smooth muscle"                            = "smooth muscle cells",
  "13 Stromal cell-2"                           = "stromal cells type 2",
  "14 Endothelial cell (venous)"                = "venous endothelial cells",
  "15 Pericyte"                                 = "pericytes",
  "16 Brown adipocyte (contractile)"            = "contractile brown adipocytes",
  "17 Stromal cell-3"                           = "stromal cells type 3",
  "18 Satelite Glia cell"                       = "satellite glia cells",
  "19 Slow muscle fiber"                        = "slow muscle fibers",
  "20 Lymphocytes"                              = "lymphocytes",
  "21 Schwann cell"                             = "schwann cells",
  "0 RBC"                                       = "RBC"
)

# --- Apply new cell type names ---
so2 <- snuc
so2 <- AddMetaData(
  so2,
  metadata = data.frame(
    cellname = recode(as.character(so2$celltype), !!!mapping),
    row.names = colnames(so2)
  )
)

# --- Define cluster order for plotting ---
cluster_order <- c(
  "basal brown adipocytes",
  "OXPHOS-high brown adipocytes",
  "lipogenic brown adipocytes",
  "contractile brown adipocytes",
  "endothelium-derived adipocytes",
  "stroma-derived adipocytes",
  "white-like adipocytes",
  "stromal cells type 1",
  "stromal cells type 2",
  "stromal cells type 3",
  "capillary endothelial cells",
  "arterial endothelial cells",
  "venous endothelial cells",
  "pericytes",
  "fast muscle fibers",
  "slow muscle fibers",
  "MonDC",
  "smooth muscle cells",
  "satellite glia cells",
  "schwann cells",
  "lymphocytes",
  "RBC"
)

# Set factor levels for consistent ordering
so2$cellname <- factor(so2$cellname, levels = cluster_order)

# --- Define colors for new cell type names ---
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

# Add cell color as metadata
so2$cellcolor <- as.character(celltype_colors[as.character(so2$cellname)])

# Check doublet distribution across cell types
table(so2$DoubletCall, so2$cellname)

# --- Save final annotated object ---
# saveRDS(so2, "data/ChREBP_SeuratObject_final.rds")

# ==============================================================================
# SECTION 6: START HERE FOR FIGURE 1 PANEL C - UMAP VISUALIZATION
# ==============================================================================
#
# Generate the main UMAP plot showing all 21 cell clusters
#
# ==============================================================================

# --- Load final object (if starting from pre-processed data) ---
snuc <- readRDS("ChREBP_SeuratObject_final.rds") # from doi: 

# --- Load cell type annotations ---
celltypes <- read.csv("data/CellTypeAnnotation_ChREBP_final.csv")
celltypes$showName <- paste0(celltypes$Number, " ", celltypes$CellType)

# --- Apply cluster names ---
clusterID <- celltypes$showName
names(clusterID) <- levels(snuc)
snuc <- RenameIdents(snuc, clusterID)
snuc$celltype <- Idents(snuc)
rownames(celltypes) <- celltypes$ClusterID

# --- Define color scheme for cell types ---
color_assignment <- setNames(celltypes$celltypecolor, celltypes$showName)

# --- Figure 1C: UMAP with cell type labels ---
# This is the main visualization showing all 21 cell clusters
figure1c <- DimPlot(snuc, 
                     cols = color_assignment, 
                     label = TRUE,
                     label.box = TRUE, 
                     repel = TRUE, 
                     label.size = 4.25) + 
  NoAxes() + 
  NoLegend() +
  ggtitle("Figure 1C: UMAP of BAT cell populations")

# Display the plot
print(figure1c)

# Save the plot
fw <- 150
ggsave("figures/Figure1C_UMAP.pdf", 
       plot = figure1c, 
       width = fw, 
       height = fw*9/16,units = "mm", dpi=300, useDingbats=F)

ggsave("figures/Figure1C_UMAP.png", 
       plot = figure1c, 
       width = fw, 
       height = fw*9/16, 
       dpi = 300)


# ==============================================================================
# SECTION 7: FIGURE 1 PANEL D - BROWN AND WHITE ADIPOCYTE SCORES
# ==============================================================================
#
# Calculate module scores for brown and white adipocyte gene signatures
# Gene lists are from BATLAS (Perdikari et al., 2018 Cell Reports)
# Human genes were converted to mouse symbols using nichenetr
#
# ==============================================================================

# --- Load color assignment ---
# This should match the colors used in Figure 1C
color_assignment <- readRDS("data/color_assignment.rds")

# --- Optional: Source custom theme function ---
# Uncomment if you have the AdlungLab helper functions
# source("path/to/AdlungLab_HelperFunctions.R")

# --- Define WATLAS (White Adipocyte) gene signature ---
# These genes distinguish white from brown adipocytes
mousewatlas <- c("Gatb", "Col4a2", "Acvr1c", "Eef2k", "Nupr1l", "Qsox1", 
                 "Nrip1", "Igf1", "Ccnd2", "Lrp1", "Col3a1", "Eepd1", "Nnat", 
                 "Cavin3", "Dmrt2", "Lpgat1", "Ccdc80", "Pik3r1", "Lep", 
                 "Pygb", "Gadd45a", "Ndrg1")

# --- Define BATLAS (Brown Adipocyte) gene signature ---
# These genes are characteristic of brown adipocytes
mousebatlas <- c("Mrps22", "Ppif", "Mrpl15", "Mrps18b", "Phb2", "C1qbp", 
                 "March5", "Timm50", "Ndufaf5", "mt-Nd2", "mt-Co3", "mt-Nd3", 
                 "Hccs", "Bckdhb", "Mrpl34", "Cyc1", "Ndufs3", "Cox5a", 
                 "Ndufs8", "Bsg", "Etfa", "Ech1", "Crls1", "Aurkaip1", 
                 "Ndufb2", "Cox7a2", "Uqcr10", "Uqcrb", "Ndufb7", "Ndufa13", 
                 "Ndufab1", "Ndufb9", "Acsf2", "Poldip2", "Slc25a11", "Flad1", 
                 "Sdc4", "Adcy3", "Ehhadh", "Glrx5", "Mrps7", "Them4", 
                 "Mrps5", "Des", "Poln", "Sdhc", "Pank1", "Ucp1", "Suclg1", 
                 "Slc25a39", "Sod2", "Idh3a", "Gk", "Chchd3", "Kcnk3", 
                 "Cpt1b", "Letmd1", "Uqcrc1", "Acsl5", "Got1", "Amacr", 
                 "Ciapin1", "Cox10", "Ndufv1", "Ndufs2", "Hadha", "Acadvl", 
                 "Idh3b", "Coq6", "Ecsit", "Uqcc1", "Ptcd3", "Timm44", 
                 "Dnajc11", "Agpat3", "Dlst", "Aco2", "Dnaja3", "Dld", 
                 "Tbrg4", "Hspa9", "Mtif2", "Acads", "Oxnad1", "Echs1", 
                 "Pck1", "Sgpl1", "Ppargc1b", "Akap1", "Acat1", "Vwa8", 
                 "Immt", "Ogdh", "Hspd1", "Sdha", "Ndufs1", "Pdha1")

# --- Join layers and calculate module scores ---
snuc <- JoinLayers(snuc)
snuc <- AddModuleScore(snuc, features = list(mousewatlas), name = "WATLAS")
snuc <- AddModuleScore(snuc, features = list(mousebatlas), name = "BATLAS")

# --- Extract metadata for plotting ---
cplot <- as.data.frame(snuc@meta.data)

# Subset to adipocyte clusters only
adipocyte_clusters <- c("1 Brown adipocyte (Ucp1 low)",
                        "2 Brown adipocyte (Ucp1 high)",
                        "4 White adipocyte",
                        "8 Brown adipocyte (lipogenic)",
                        "9 Brown adipocyte (from endothelial cell)",
                        "10 Brown adipocyte (from stromal cell)",
                        "16 Brown adipocyte (contractile)")

cplot2 <- cplot[which(cplot$celltype %in% adipocyte_clusters), ]

# --- Define labels for legend ---
adipocyte_labels <- c("basal brown adipocytes",
                      "OXPHOS-high brown adipocytes",
                      "white-like adipocytes",
                      "lipogenic adipocytes",
                      "endothelial-derived adipocytes",
                      "stromal-derived adipocytes",
                      "contractile brown adipocytes")

# --- Figure 1D: Brown adipocyte score (BATLAS) ---
pb <- ggplot(cplot2, aes(x = celltype, y = BATLAS1, fill = celltype)) + 
  geom_violin(color = "#000000", 
              draw_quantiles = c(0.5), 
              trim = TRUE, 
              scale = "width",
              width = 0.66) + 
  labs(title = NULL, 
       x = NULL, 
       y = "Brown-adipocyte score\n(normalized expression)") +
  scale_fill_manual(values = color_assignment, 
                    name = NULL,
                    labels = adipocyte_labels) +
  theme_adlunglab(base_size = 8) + 
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "right")

# --- Figure 1D: White adipocyte score (WATLAS) ---
pa <- ggplot(cplot2, aes(x = celltype, y = WATLAS1, fill = celltype)) + 
  geom_violin(color = "#000000", 
              draw_quantiles = c(0.5), 
              trim = TRUE, 
              scale = "width",
              width = 0.66) + 
  labs(title = NULL, 
       x = NULL, 
       y = "White-adipocyte score\n(normalized expression)") +
  scale_fill_manual(values = color_assignment, 
                    name = NULL,
                    labels = adipocyte_labels) +
  theme_adlunglab(base_size = 8) + 
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "right")

# --- Combine plots using ggpubr ---
library(ggpubr)

figure1d <- ggarrange(pb, pa, 
                      ncol = 2, 
                      labels = c("D", ""), 
                      common.legend = TRUE, 
                      legend = "right",
                      widths = c(1/2, 1/2))

# Display
print(figure1d)

# Save
fontlist <- list("max"=c(10), "min"=c(7))
widthlist <- list("min"=c(30), "singlecolumn"=c(90),"oneandahalf"=c(140), "fullwidth"=c(190))

fw <- widthlist$fullwidth
ggsave(paste0("figures/Figure1D_Scores.pdf"),
       width=fw,height=fw*9/16/1.5,units = "mm", dpi=300, useDingbats=FALSE
)

# ==============================================================================
# SECTION 8: FIGURE 1 PANEL E - MARKER GENE DOT PLOT
# ==============================================================================
#
# Generate dot plot showing average expression and percent expressed
# for key marker genes across adipocyte clusters
# Data shown is from Cre- (Wild-Type) mice at room temperature only
#
# ==============================================================================

# --- Ensure layers are joined (required for Seurat v5) ---
snuc <- JoinLayers(snuc)

# --- Load color assignment ---
color_assignment <- readRDS("data/color_assignment.rds")

# --- Define marker genes for dot plot ---
# Genes are grouped by cell type characteristics:
# - Basal brown adipocytes: Ucp1, Prdm16, Esrrg, Ppargc1a, Ppargc1b, Cidea, Dio2
# - OXPHOS-high brown adipocytes: Fabp4, Scd1
# - White-like adipocytes: Retn, Lep, Adipoq, Cyp2e1, Slc7a10, Aldh1a1, Slc2a3, Tshr
# - Lipogenic adipocytes: Mlxipl, Pparg, Srebf1, Nr1h2, Nr1h3, Fasn, Acaca, Bhlhe40
# - Contractile brown adipocytes: Myh2, Myh7, Tnnc2, Tnni2, Ttn

adipocyte.bubblegenes <- c("Ucp1", "Prdm16", "Esrrg", "Ppargc1a", "Ppargc1b", 
                           "Cidea", "Dio2",
                           "Fabp4", "Scd1",
                           "Retn", "Lep", "Adipoq",
                           "Cyp2e1", "Slc7a10", "Aldh1a1", "Slc2a3", "Tshr",
                           "Mlxipl", "Pparg", "Srebf1", "Nr1h2", "Nr1h3", 
                           "Fasn", "Acaca", "Bhlhe40",
                           "Myh2", "Myh7", "Tnnc2", "Tnni2", "Ttn")

# --- Subset to Wild-Type Room Temperature sample only ---
Idents(snuc) <- "sample"
wtrt <- subset(snuc, idents = c("01_Wild-Type_Room Temp."))

# Check dimensions
message(paste0("WT RT subset: ", ncol(wtrt), " cells"))

# --- Subset to adipocyte clusters only ---
Idents(wtrt) <- "celltype"

cells.adipocytes <- c("1 Brown adipocyte (Ucp1 low)",
                      "2 Brown adipocyte (Ucp1 high)",
                      "4 White adipocyte",
                      "8 Brown adipocyte (lipogenic)",
                      "9 Brown adipocyte (from endothelial cell)",
                      "10 Brown adipocyte (from stromal cell)",
                      "16 Brown adipocyte (contractile)")

labels.adipocytes <- c("basal brown\nadipocytes",
                       "OXPHOS-high brown\nadipocytes",
                       "white-like\nadipocytes",
                       "lipogenic\nadipocytes",
                       "endothelial-derived\nadipocytes",
                       "stromal-derived\nadipocytes",
                       "contractile brown\nadipocytes")

wtrt.adi <- subset(wtrt, idents = cells.adipocytes)

# Check dimensions
message(paste0("WT RT adipocytes: ", ncol(wtrt.adi), " cells"))

# --- Relevel cell types with custom labels (reversed for y-axis ordering) ---
wtrt.adi$celltype <- factor(wtrt.adi$celltype, 
                            levels = rev(cells.adipocytes),
                            labels = rev(labels.adipocytes))
Idents(wtrt.adi) <- "celltype"

# --- Define figure dimensions ---
# Using publication-ready dimensions (full width = 190mm)
widthlist <- list("min" = 30, "singlecolumn" = 90, "fullwidth" = 190)
fw <- widthlist$fullwidth

# --- Figure 1E: Dot plot of marker genes ---
library(ggplot2)

s2c <- DotPlot(wtrt.adi, 
               features = adipocyte.bubblegenes,
               cols = c("grey89", "#220589")) +
  theme_classic(base_size = 7) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 33, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7)) +
  labs(title = NULL, x = NULL, y = NULL)

# Display
print(s2c)

# --- Combine with label using ggpubr ---
library(ggpubr)

figure1e <- ggarrange(s2c, 
                      ncol = 1, 
                      labels = c("E"), 
                      legend = "bottom",
                      widths = c(1))

# Display
print(figure1e)

# --- Save with exact dimensions from original ---
ggsave("figures/Figure1E_DotPlot_Adipocytes.pdf",
       plot = figure1e,
       width = fw, 
       height = fw * 9/16, 
       units = "mm", 
       dpi = 300, 
       useDingbats = FALSE)

ggsave("figures/Figure1E_DotPlot_Adipocytes.png",
       plot = figure1e,
       width = fw, 
       height = fw * 9/16, 
       units = "mm", 
       dpi = 300)
# ==============================================================================
# SESSION INFO
# ==============================================================================

sessionInfo()

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
