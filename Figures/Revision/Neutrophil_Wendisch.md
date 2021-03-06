Neutrophil\_Wendisch
================
Sam Freeman

Neutrophil Wendisch et al. data analysis

Load the necessary libraries:

``` r
library(Seurat)
library(ggplot2)
library(batchelor)
library(leiden)
library(fgsea)
library(openxlsx)
library(plyr)
library(stringr)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)
```

``` r
prefix <- "/Volumes/GoogleDrive/.shortcut-targets-by-id/1SQXfCUGIenBLXc4w_p0n7Mufi3JwjFCN/COVID19_Neutrophils/Revision/Final_Steps/Github/"
```

# Download Wendisch et al. BAL scRNA data

# `{r} #con <- "https://nubes.helmholtz-berlin.de/s/XrM8igTzFTFSoio/download" #`

Load seurat object, plot UMAP with all cell types

``` r
object <- readRDS(paste0(prefix,"BAL.Rds"))

g_all_umap <- Seurat::DimPlot(object, group.by = "Celltype") +
  ggplot2::scale_color_brewer(palette = "Set1") +
  ggplot2::coord_fixed()

g_all_umap
```

![](Neutrophil_Wendisch_files/figure-markdown_github/unnamed-chunk-3-1.png)

Subset Seurat object to Neutrophils, renormalize and identify variable genes

``` r
#Subset
cells <- colnames(object)[which(
  object@meta.data$Celltype == "Neutrophils"
)]

tcounts <- GetAssayData(object@assays[["RNA"]], slot = "counts")[,cells]
genes.use <- rowSums(tcounts) > 0
genes <- names(genes.use[genes.use])

# Filter
object <- subset(object, cells = cells, features = genes)

# Normalize mRNA counts
object <- Seurat::NormalizeData(
  object               = object,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)

# Find HVGs
object <- Seurat::FindVariableFeatures(
  object           = object,
  selection.method = "vst", 
  nfeatures        = 3000
)
```

Run UMAP for Neutrophils

``` r
# Calculate low-dimensional embedding & clustering

# MNN-corrected PCA
set.seed(1993)
object@reductions[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = batchelor::fastMNN(
    object@assays$RNA@data[object@assays$RNA@var.features, ],
    batch = object@meta.data$patient,
    d     = 15
  )@int_colData$reducedDims$corrected,
  key        = "PC_",
  assay      = "RNA"
)
```

    ## Warning: No columnames present in cell embeddings, setting to 'PC_1:15'

``` r
# UMAP
object <- Seurat::RunUMAP(object = object, dims = 1:15, seed.use = 1993)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

``` r
# Clustering
set.seed(1993)
object <- Seurat::FindNeighbors(object = object, dims = 1:15)
object@meta.data$Cluster <- factor(leiden::leiden(
  object@graphs$RNA_snn, resolution_parameter = 0.6, seed = 1993
))

g_neu_umap <- Seurat::DimPlot(object, group.by = "Cluster") + ggplot2::scale_color_brewer(palette = "Set1") + ggplot2::coord_fixed()
#ggsave("wendisch.neu.umap.pdf",g_neu_umap)
g_neu_umap
```

![](Neutrophil_Wendisch_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
genepc <- read.delim(paste0(prefix,"Ensembl_to_Symbol.txt"))
bal_featureplot <- FeaturePlot(object, features = c(genepc[genepc$Gene.name=="S100A8",]$Gene.stable.ID,genepc[genepc$Gene.name=="S100A9",]$Gene.stable.ID,genepc[genepc$Gene.name=="S100A12",]$Gene.stable.ID,genepc[genepc$Gene.name=="FCGR2A",]$Gene.stable.ID,genepc[genepc$Gene.name=="CEACAM8",]$Gene.stable.ID,genepc[genepc$Gene.name=="CXCR2",]$Gene.stable.ID), raster=TRUE)
tmp_geneids <- c("S100A8","S100A9","S100A12","FCGR2A","CEACAM8","CXCR2")
for (i in 1: length(tmp_geneids)) bal_featureplot[[i]]$labels$title = tmp_geneids[i]
#ggsave("BAL_featureplot.pdf",bal_featureplot)
bal_featureplot
```

![](Neutrophil_Wendisch_files/figure-markdown_github/unnamed-chunk-5-2.png)

Score all BAL neutrophils by NMF signatures

``` r
gmt.file.ensg <- gmtPathways(paste0(prefix,"neutrophil_state_gene_sets_ensembl.gmt"))

tab <- object@assays$RNA@counts[object@assays$RNA@counts@Dimnames[[1]] %in% gmt.file.ensg$NMF1,]
tab_z <- apply(tab,1,scale)
tab_z <- t(tab_z)
nmf1_meanz <- as.numeric(colMeans(tab_z))
object <- AddMetaData(object, metadata = nmf1_meanz, col.name = "nmf1")
tab <- object@assays$RNA@counts[object@assays$RNA@counts@Dimnames[[1]] %in% gmt.file.ensg$NMF2,]
tab_z <- apply(tab,1,scale)
tab_z <- t(tab_z)
nmf2_meanz <- as.numeric(colMeans(tab_z))
object <- AddMetaData(object, metadata = nmf2_meanz, col.name = "nmf2")
tab <- object@assays$RNA@counts[object@assays$RNA@counts@Dimnames[[1]] %in% gmt.file.ensg$NMF3,]
tab_z <- apply(tab,1,scale)
tab_z <- t(tab_z)
nmf3_meanz <- as.numeric(colMeans(tab_z))
object <- AddMetaData(object, metadata = nmf3_meanz, col.name = "nmf3")
tab <- object@assays$RNA@counts[object@assays$RNA@counts@Dimnames[[1]] %in% gmt.file.ensg$NMF4,]
tab_z <- apply(tab,1,scale)
tab_z <- t(tab_z)
nmf4_meanz <- as.numeric(colMeans(tab_z))
object <- AddMetaData(object, metadata = nmf4_meanz, col.name = "nmf4")
tab <- object@assays$RNA@counts[object@assays$RNA@counts@Dimnames[[1]] %in% gmt.file.ensg$NMF5,]
tab_z <- apply(tab,1,scale)
tab_z <- t(tab_z)
nmf5_meanz <- as.numeric(colMeans(tab_z))
object <- AddMetaData(object, metadata = nmf5_meanz, col.name = "nmf5")
tab <- object@assays$RNA@counts[object@assays$RNA@counts@Dimnames[[1]] %in% gmt.file.ensg$NMF6,]
tab_z <- apply(tab,1,scale)
tab_z <- t(tab_z)
nmf6_meanz <- as.numeric(colMeans(tab_z))
object <- AddMetaData(object, metadata = nmf6_meanz, col.name = "nmf6")

#plot signatures scores on BAL Neutrophil UMAP
p0 <- DimPlot(object, reduction = "umap", group.by = "Cluster", label = TRUE, raster = TRUE)
p0 <- p0 + theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) + ggtitle("")
p0
```

![](Neutrophil_Wendisch_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
p1 <- FeaturePlot(object, reduction = "umap", features = c("nmf1"), max.cutoff = 1.5, cols = c("blue","yellow","red"), raster = TRUE) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds"))) 
p1 <- p1 + theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "none") + ggtitle("NMF1")
p2 <- FeaturePlot(object, reduction = "umap", features = c("nmf2"), max.cutoff = 1.5, cols = c("blue","yellow","red"), raster = TRUE) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds"))) 
p2 <- p2 + theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "none") + ggtitle("NMF2")
p3 <- FeaturePlot(object, reduction = "umap", features = c("nmf3"), max.cutoff = 1.5, cols = c("blue","yellow","red"), raster = TRUE) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
p3 <- p3 + theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "none") + ggtitle("NMF3")
p4 <- FeaturePlot(object, reduction = "umap", features = c("nmf4"), max.cutoff = 1.5, cols = c("blue","yellow","red"), raster = TRUE) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
p4 <- p4 + theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "none") + ggtitle("NMF4")
p5 <- FeaturePlot(object, reduction = "umap", features = c("nmf5"), max.cutoff = 1.5, cols = c("blue","yellow","red"), raster = TRUE) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))
p5 <- p5 + theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "none") + ggtitle("NMF5")
p6 <- FeaturePlot(object, reduction = "umap", features = c("nmf6"), max.cutoff = 1.5, cols = c("blue","yellow","red"), raster = TRUE) & scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds"))) 
p6 <- p6 + theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "none") + ggtitle("NMF6")

cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol=3)
```

![](Neutrophil_Wendisch_files/figure-markdown_github/unnamed-chunk-6-2.png)

``` r
#cowplot::save_plot("BAL_neutrophil_NMF_UMAP.pdf",cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol=3))

#plot NMF signature split by patient status (alive/dead)
#patients C19-82 and C19-83 had WHO status 8 (death) and patients C19-62,C19-85,C19-98,C19-120 and C19-136 had WHO status 7 (severe disease)
object <- AddMetaData(object,ifelse(object$patient %in% c("C19-82","C19-83"),"dead","alive"),col.name="status")

g_neu_umap_alivedead <- Seurat::DimPlot(object, group.by = "status") + ggplot2::scale_color_manual(values= c("blue","red")) + ggplot2::coord_fixed()
#ggsave("wendisch.neu.umap.alivedead.pdf",g_neu_umap_alivedead)
g_neu_umap_alivedead
```

![](Neutrophil_Wendisch_files/figure-markdown_github/unnamed-chunk-6-3.png)

``` r
df_bal <- data.frame(NMF1=object@meta.data$nmf1,
                     NMF2=object@meta.data$nmf2,
                     NMF3=object@meta.data$nmf3,
                     NMF4=object@meta.data$nmf4,
                     NMF5=object@meta.data$nmf5,
                     NMF6=object@meta.data$nmf6,
                     patient=object@meta.data$patient,
                     status=object@meta.data$status,stringsAsFactors=F)

p_nmf1_ad <- ggplot(df_bal, aes(x = factor(status), y = NMF1, fill = factor(status))) + geom_boxplot(outlier.shape = NA) + theme_bw() + xlab("Status") + ylab("NMF1") + scale_fill_manual(values = c("blue","red")) + coord_fixed(ratio = 0.5) + coord_cartesian(ylim = c(-0.5,1)) + stat_compare_means(label = "p.signif",label.y=0.5,label.x=1.4) + theme(legend.position="none",plot.title=element_text(hjust = 0.5)) + ggtitle("NMF1")
p_nmf2_ad <- ggplot(df_bal, aes(x = factor(status), y = NMF2, fill = factor(status))) + geom_boxplot(outlier.shape = NA) + theme_bw() + xlab("Status") + ylab("NMF2") + scale_fill_manual(values = c("blue","red")) + coord_fixed(ratio = 0.5) + coord_cartesian(ylim = c(-0.5,1)) + stat_compare_means(label = "p.signif",label.y=0.5,label.x=1.4) + theme(legend.position="none",plot.title=element_text(hjust = 0.5)) + ggtitle("NMF2")
p_nmf3_ad <- ggplot(df_bal, aes(x = factor(status), y = NMF3, fill = factor(status))) + geom_boxplot(outlier.shape = NA) + theme_bw() + xlab("Status") + ylab("NMF3") + scale_fill_manual(values = c("blue","red")) + coord_fixed(ratio = 0.5) + coord_cartesian(ylim = c(-0.5,1)) + stat_compare_means(label = "p.signif",label.y=0.5,label.x=1.4) + theme(legend.position="none",plot.title=element_text(hjust = 0.5)) + ggtitle("NMF3")
p_nmf4_ad <- ggplot(df_bal, aes(x = factor(status), y = NMF4, fill = factor(status))) + geom_boxplot(outlier.shape = NA) + theme_bw() + xlab("Status") + ylab("NMF4") + scale_fill_manual(values = c("blue","red")) + coord_fixed(ratio = 0.5) + coord_cartesian(ylim = c(-0.5,1)) + stat_compare_means(label = "p.signif",label.y=0.5,label.x=1.4) + theme(legend.position="none",plot.title=element_text(hjust = 0.5)) + ggtitle("NMF4")
p_nmf5_ad <- ggplot(df_bal, aes(x = factor(status), y = NMF5, fill = factor(status))) + geom_boxplot(outlier.shape = NA) + theme_bw() + xlab("Status") + ylab("NMF5") + scale_fill_manual(values = c("blue","red")) + coord_fixed(ratio = 0.5) + coord_cartesian(ylim = c(-0.5,1)) + stat_compare_means(label = "p.signif",label.y=0.5,label.x=1.4) + theme(legend.position="none",plot.title=element_text(hjust = 0.5)) + ggtitle("NMF5")
p_nmf6_ad <- ggplot(df_bal, aes(x = factor(status), y = NMF6, fill = factor(status))) + geom_boxplot(outlier.shape = NA) + theme_bw() + xlab("Status") + ylab("NMF6") + scale_fill_manual(values = c("blue","red")) + coord_fixed(ratio = 0.5) + coord_cartesian(ylim = c(-0.5,1)) + stat_compare_means(label = "p.signif",label.y=0.5,label.x=1.4) + theme(legend.position="none",plot.title=element_text(hjust = 0.5)) + ggtitle("NMF6")

g_grid <- grid.arrange(p_nmf1_ad,p_nmf2_ad,p_nmf3_ad,p_nmf4_ad,p_nmf5_ad,p_nmf6_ad,nrow=1)
```

![](Neutrophil_Wendisch_files/figure-markdown_github/unnamed-chunk-6-4.png)

``` r
#ggsave("BAL_NMF_signature_alive_vs_dead.pdf",g_grid,height=7,width=12)
g_grid
```

    ## TableGrob (1 x 6) "arrange": 6 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (1-1,2-2) arrange gtable[layout]
    ## 3 3 (1-1,3-3) arrange gtable[layout]
    ## 4 4 (1-1,4-4) arrange gtable[layout]
    ## 5 5 (1-1,5-5) arrange gtable[layout]
    ## 6 6 (1-1,6-6) arrange gtable[layout]
