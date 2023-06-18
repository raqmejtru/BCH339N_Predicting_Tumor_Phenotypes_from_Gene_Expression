# //////////////////////////////////////////////////////////////////////////////
#
# Bioinformatics Project
# Fall 2022
# Raquel Mejia-Trujillo | rm57578
#
# Objective:
# Which biological processes characterize breast cancer, skin cancer, low grade
# glioma, and mesothelioma tumors based on differences in gene expression?
#
# //////////////////////////////////////////////////////////////////////////////


# 0. Load libraries ------------------------------------------------------------

library(TCGAbiolinks) # Bioconductor - Load TCGA RNASeq data sets
library(DESeq2)       # Bioconductor - DGE analysis
library(fgsea)        # Bioconductor - Gene set enrichment
library(msigdbr)      # Bioconductor - Database of hallmark gene sets
library(tidyverse)    # CRAN         - Data manipulation
library(ggplot2)      # CRAN         - Plotting
library(umap)         # CRAN         - UMAP Projections

setwd("./BCH339N_Predicting_Tumor_Phenotypes_from_Gene_Expression")




# 1. Download & save data from TCGA --------------------------------------------

# Breast cancer (BRCA) data
brca.query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor")
    )

GDCdownload(
    query = brca.query,
    files.per.chunk = 100
    )

brca.exp <- GDCprepare(
    query = brca.query,
    save = TRUE,
    save.filename = "./data/brca.rda"
    )

# Skin Cutaneous Melanoma (SKCM) data
skcm.query <- GDCquery(
    project = "TCGA-SKCM",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor")
    )

GDCdownload(
    query = skcm.query,
    files.per.chunk = 100
    )

skcm.exp <- GDCprepare(
    query = skcm.query,
    save = TRUE,
    save.filename = "./data/skcm.rda"
    )

# Lower grade glioma cancer (LGG) data
lgg.query <- GDCquery(
    project = "TCGA-LGG",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor")
    )

GDCdownload(
    query = lgg.query,
    files.per.chunk = 100
    )

lgg.exp <- GDCprepare(
    query = lgg.query,
    save = TRUE,
    save.filename = "./data/lgg.rda"
    )

# Mesothelioma cancer (MESO) data
meso.query <- GDCquery(
    project = "TCGA-MESO",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor")
    )

GDCdownload(
    query = meso.query,
    files.per.chunk = 100
    )

meso.exp <- GDCprepare(
    query = meso.query,
    save = TRUE,
    save.filename = "./data/meso.rda"
    )



# 2. Construct a combined count matrix -----------------------------------------

# Save raw counts as df with samples as rows and genes as cols
skcm.counts = as.data.frame(t(assay(skcm.exp)))
lgg.counts = as.data.frame(t(assay(lgg.exp)))
meso.counts = as.data.frame(t(assay(meso.exp)))

# BRCA had 1000+ samples and was crashing DESeq Analyses
# Went back and chose a subset of 500 samples
# Note: However, UMAP Viz includes all BRCA samples
set.seed(1)
brca.counts = as.data.frame(t(assay(brca.exp)))
subset.brca.counts = runif(500, min = 1, max = nrow(brca.counts))
brca.counts = brca.counts[subset.brca.counts, ]

# Add column with tumor classification
brca.counts$tumor = "brca"
skcm.counts$tumor = "skcm"
lgg.counts$tumor = "lgg"
meso.counts$tumor = "meso"

# Merge rows into one data.frame
merged.counts = rbind(brca.counts, skcm.counts, lgg.counts, meso.counts)



# 3. Define Sample Conditions for DESeq Objects --------------------------------

# Transform counts so rows are genes and samples are columns and remove tumor type columns
merged.counts.genes.only = merged.counts[,1:60660]
cancer.counts = t(merged.counts.genes.only)

### For UMAP:
# Dataframe of sample tumor conditions
tumor.samples = NULL
tumor.samples$condition = merged.counts$tumor
tumor.samples$sample = rownames(merged.counts)
tumor.samples = as.data.frame(tumor.samples)

### For DGE Analysis:
# BRCA vs. not BRCA
brca.deseq = tumor.samples
brca.deseq[brca.deseq$condition != "brca", "condition"] = "not_brca"
brca.deseq$condition = factor(brca.deseq$condition, levels = c("brca", "not_brca"))
levels(brca.deseq$condition)

# SKCM vs. not SKCM
skcm.deseq = tumor.samples
skcm.deseq[skcm.deseq$condition != "skcm", "condition"] = "not_skcm"
skcm.deseq$condition = factor(skcm.deseq$condition, levels = c("skcm", "not_skcm"))
levels(skcm.deseq$condition)

# LGG vs. not LGG
lgg.deseq = tumor.samples
lgg.deseq[lgg.deseq$condition != "lgg", "condition"] = "not_lgg"
lgg.deseq$condition = factor(lgg.deseq$condition, levels = c("lgg", "not_lgg"))
levels(lgg.deseq$condition)

# MESO vs. not MESO
meso.deseq = tumor.samples
meso.deseq[meso.deseq$condition != "meso", "condition"] = "not_meso"
meso.deseq$condition = factor(meso.deseq$condition, levels = c("meso", "not_meso"))
levels(meso.deseq$condition)

# Convert back to factor
tumor.samples$condition = factor(merged.counts$tumor)



# 4. Make DESeq Objects  -------------------------------------------------------

# To normalize counts for UMAP Visualization
deseq.input =  DESeqDataSetFromMatrix(
    countData = cancer.counts,
    colData = tumor.samples,
    design = ~ condition
    )

# For BRCA DGE
brca.deseq.input = DESeqDataSetFromMatrix(
    countData = cancer.counts,
    colData = brca.deseq,
    design = ~ condition
    )

# For SKCM DGE
skcm.deseq.input = DESeqDataSetFromMatrix(
    countData = cancer.counts,
    colData = skcm.deseq,
    design = ~ condition
    )

# For LGG DGE
lgg.deseq.input  = DESeqDataSetFromMatrix(
    countData = cancer.counts,
    colData = lgg.deseq,
    design = ~ condition
    )

# For MESO DGE
meso.deseq.input = DESeqDataSetFromMatrix(
    countData = cancer.counts,
    colData = meso.deseq,
    design = ~ condition
    )



# 5. UMAP ----------------------------------------------------------------------

# Variance Stabilizing Transformation
vst.norm.counts = varianceStabilizingTransformation(deseq.input, blind = F)

# Transpose so rows are samples and genes are columns
vst.norm.counts.t = t(assay(vst.norm.counts))

# Perform UMAP on the normalized counts
umap.output = umap(vst.norm.counts.t)
umap.coordinates = as.data.frame(umap.output$layout)
umap.coordinates$sample = rownames(umap.coordinates)
umap.coordinates = inner_join(umap.coordinates, tumor.samples)

# Plot 2D Projection
umap.plot = umap.coordinates %>%
    ggplot() +
    geom_point(aes(x = V1, y = V2, col = condition)) +
    theme_minimal() +
    labs(
        title = "2-D UMAP Projection of Gene Expression",
        col = "",
        x = "UMAP_1",
        y = "UMAP_2"
        ) +
    scale_color_viridis_d(
        labels = c(
            "Breast Cancer",
            "Low Grade Glioma",
            "Mesothelioma",
            "Skin Cutaneous Melanoma"
            )
        ) +
    theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold")
    )

# Export plot
ggsave(
    umap.plot,
    "./figures/umap_plot.png"
    )


# 6. Differential Gene Expression ----------------------------------------------

# Run DESeq Analysis
brca.deseq.output = DESeq(brca.deseq.input)
skcm.deseq.output = DESeq(skcm.deseq.input)
lgg.deseq.output  = DESeq(lgg.deseq.input)
meso.deseq.output = DESeq(meso.deseq.input)

brca.results = results(brca.deseq.output, tidy = T)
skcm.results = results(skcm.deseq.output, tidy = T)
lgg.results  = results(lgg.deseq.output,  tidy = T)
meso.results = results(meso.deseq.output, tidy = T)



# 7. Gene Set Enrichment (GSE) -------------------------------------------------

# Remove gene name version suffix ".version" (ie, only ENSG#####)
brca.results$row = gsub("\\..+", "", brca.results$row)
skcm.results$row = gsub("\\..+", "", skcm.results$row)
lgg.results$row  = gsub("\\..+", "", lgg.results$row)
meso.results$row = gsub("\\..+", "", meso.results$row)

# Rank genes for fgsea input
brca.stat = brca.results %>%
    dplyr::select(row, stat) %>%
    na.omit() %>%
    distinct %>%
    group_by(row) %>%
    summarize(stat = mean(stat)) %>%
    mutate(rank = rank(stat, ties.method = "random")) %>%
    arrange(desc(rank)) %>%
    dplyr::select(row, stat)

skcm.stat = skcm.results %>%
    dplyr::select(row, stat) %>%
    na.omit() %>%
    distinct %>%
    group_by(row) %>%
    summarize(stat = mean(stat)) %>%
    mutate(rank = rank(stat, ties.method = "random")) %>%
    arrange(desc(rank)) %>%
    dplyr::select(row, stat)

lgg.stat = lgg.results %>%
    dplyr::select(row, stat) %>%
    na.omit() %>%
    distinct %>%
    group_by(row) %>%
    summarize(stat = mean(stat)) %>%
    mutate(rank = rank(stat, ties.method = "random")) %>%
    arrange(desc(rank)) %>%
    dplyr::select(row, stat)

meso.stat = meso.results %>%
    dplyr::select(row, stat) %>%
    na.omit() %>%
    distinct %>%
    group_by(row) %>%
    summarize(stat = mean(stat)) %>%
    mutate(rank = rank(stat, ties.method = "random")) %>%
    arrange(desc(rank)) %>%
    dplyr::select(row, stat)

# Convert ENSEMBL column to row names, then convert to list
brca.ranks = deframe(brca.stat)
skcm.ranks = deframe(skcm.stat)
lgg.ranks  = deframe(lgg.stat)
meso.ranks = deframe(meso.stat)

# Load hallmark gene sets (well-defined biological pathways)
hallmark.pathways = msigdbr(species = "human", category = "H")
hallmark.pathways = split(hallmark.pathways$human_ensembl_gene, hallmark.pathways$gs_name)



# 7.1 BRCA GSE
brca.fgsea = fgsea(pathways = hallmark.pathways, stats = brca.ranks)

# Clean pathway names
brca.fgsea$pathway = gsub("HALLMARK_", "", brca.fgsea$pathway)

# Plot BRCA GSE pathways
brca.gsea.plot = brca.fgsea %>%
    ggplot() +
    geom_col(
        aes(
            reorder(pathway, NES),
            NES,
            fill = padj < 0.05
            )
        ) +
    coord_flip() +
    labs(
        title = "Gene Set Enrichment Based on Differential Expression of BRCA Tumors",
        x = "Pathways",
        y = "Normalized Enrichment Score"
        ) +
    theme_minimal() +
    theme(title = element_text(hjust = 0.5, face = "bold"))

# Export plot
ggsave(
    brca.gsea.plot,
    "./figures/brca_gsea_plot.png"
    )


# 7.2 SKCM GSE
skcm.fgsea = fgsea(pathways = hallmark.pathways, stats = skcm.ranks)

# Clean pathway names
skcm.fgsea$pathway = gsub("HALLMARK_", "", skcm.fgsea$pathway)

# Plot SKCM GSE pathways
skcm.gsea.plot = skcm.fgsea %>%
    ggplot() +
    geom_col(
        aes(
            reorder(pathway, NES),
            NES,
            fill = padj < 0.05
            )
        ) +
    coord_flip() +
    labs(
        title = "Gene Set Enrichment Based on Differential Expression of SKCM Tumors",
        x = "Pathways",
        y = "Normalized Enrichment Score"
        ) +
    theme_minimal() +
    theme(title = element_text(hjust = 0.5, face = "bold"))

# Export plot
ggsave(
    skcm.gsea.plot,
    "./figures/skcm_gsea_plot.png"
)


# 7.3 LGG GSE
lgg.fgsea = fgsea(pathways = hallmark.pathways, stats = lgg.ranks)

# Clean pathway names
lgg.fgsea$pathway = gsub("HALLMARK_", "", lgg.fgsea$pathway)

# Plot LGG GSE pathways
lgg.gsea.plot = lgg.fgsea %>%
    ggplot() +
    geom_col(
        aes(
            reorder(pathway, NES),
            NES,
            fill = padj < 0.05)
        ) +
    coord_flip() +
    labs(
        title = "Gene Set Enrichment Based on Differential Expression of LGG Tumors",
        x = "Pathways",
        y = "Normalized Enrichment Score"
        ) +
    theme_minimal() +
    theme(title = element_text(hjust = 0.5, face = "bold"))

# Export plot
ggsave(
    lgg.gsea.plot,
    "./figures/lgg_gsea_plot.png"
)


# 7.4 MESO GSE
meso.fgsea = fgsea(pathways = hallmark.pathways, stats = meso.ranks)

# Clean pathway names
meso.fgsea$pathway = gsub("HALLMARK_", "", meso.fgsea$pathway)

# Plot MESO GSE pathways
meso.gsea.plot = meso.fgsea %>%
    ggplot() +
    geom_col(
        aes(
            reorder(pathway, NES),
            NES,
            fill = padj < 0.05
            )
        ) +
    coord_flip() +
    labs(
        title = "Gene Set Enrichment Based on Differential Expression of MESO Tumors",
        x = "Pathways",
        y = "Normalized Enrichment Score"
        ) +
    theme_minimal() +
    theme(title = element_text(hjust = 0.5, face = "bold"))

# Export plot
ggsave(
    meso.gsea.plot,
    "./figures/meso_gsea_plot.png"
)