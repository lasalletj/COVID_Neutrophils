Figure4\_S14-S16
================
Tom LaSalle

This document contains all the code necessary to generate the plots for
Figure 4 and related supplementary figures (S14-S16). Plots are
subsequently edited in Adobe Illustrator to produce the final figures.

Load the necessary libraries:

``` r
library(knitr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(DESeq2)
library(openxlsx)
library(cowplot)
library(fgsea)
library(ggpubr)
```

Load in the data:

``` r
prefix <- "~/Downloads/Github/"
metadata_long <- read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 4)
Count <- read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 7, rowNames = TRUE)
TPM <- read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 8, rowNames = TRUE)
qc_data <- read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 9)
genomic_signatures <- read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 10)
genepc <- read.delim(paste0(prefix,"Ensembl_to_Symbol.txt"))
logTPM <- log2(TPM + 1)

metadata_long <- metadata_long[which(metadata_long$Public.ID %in% qc_data$Public.ID),]
metadata_long <- merge(metadata_long, qc_data, all.y = TRUE)

metadata_filtered <- metadata_long[metadata_long$percent.mt < 20 & metadata_long$Genes.Detected > 10000 & metadata_long$Median.Exon.CV < 1 & metadata_long$Exon.CV.MAD < 0.75 & metadata_long$Exonic.Rate*100 > 25 & metadata_long$Median.3..bias < 0.9,]

logTPM_filtered <- logTPM[,colnames(logTPM) %in% metadata_filtered$Public.Sample.ID]
TPM_filtered <- TPM[,colnames(TPM) %in% metadata_filtered$Public.Sample.ID]
Count_filtered <- Count[,colnames(Count) %in% metadata_filtered$Public.Sample.ID]

tf <- rowSums(TPM_filtered > 0.1) > ncol(TPM_filtered)*.2
TPM_filtered <- TPM_filtered[tf,]
Count_filtered <- Count_filtered[tf,]
logTPM_filtered <- logTPM_filtered[tf,]
tf <- rowSums(Count_filtered >= 6) > ncol(Count_filtered)*.2
TPM_filtered <- TPM_filtered[tf,]
Count_filtered <- Count_filtered[tf,]
logTPM_filtered <- logTPM_filtered[tf,]

metadata_filtered <- merge(metadata_filtered, genomic_signatures)
metadata_filtered$Public.Sample.ID <- metadata_filtered$Public.Sample.ID
metadata_filtered$COVID <- mapvalues(metadata_filtered$COVID, from = c(0,1), to = c("Negative","Positive"))

# Color Palette
vermillion <- rgb(213,94,0,max=255)
bluishgreen <- rgb(0,158,115,max=255)
yellow <- rgb(240,228,66,max=255)
blue <- rgb(0,114,178,max=255)
orange <- rgb(230,159,0,max=255)
skyblue <- rgb(86,180,233,max=255)
lightgray <- rgb(211,211,211,max=255)
```

In Figure 4 we explore signatures of neutrophil effector functions. We
begin by defining a transcriptional signature of virally-induced NETosis
composed of the genes *PADI4*, *MPO*, *ELANE*, *TNF*, *CXCL8*, *GSDMD*,
*TLR3*. First we check the pairwise correlations between these genes.

``` r
netosisgenes <- c("PADI4_logTPM","MPO_logTPM","ELANE_logTPM","TNF_logTPM","CXCL8_logTPM","GSDMD_logTPM","TLR3_logTPM")

metadata_temp <- metadata_filtered[metadata_filtered$COVID == "Positive",]

my.cols <- brewer.pal(3,"RdBu")
upper.panel<-function(x, y){
  points(x,y, pch=19, col=mapvalues(metadata_temp$severity.max, from = c("non-severe","severe"), to = c(my.cols[2],my.cols[1])))
  r <- round(cor(as.numeric(x), as.numeric(y)), digits=2)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.5, 0.9, txt)
}
```

**Figure S14A:**

``` r
pairs(sapply(metadata_temp[,colnames(metadata_temp) %in% netosisgenes],as.numeric), pch = 19, upper.panel = NULL, cex = 0.5, col = mapvalues(metadata_temp$severity.max, from = c("non-severe","severe"), to = c(my.cols[3],my.cols[1])))
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Then we score each sample according to these genes.

``` r
source(paste0(prefix,"Pathway_scoring.R"))

netosisgenes <- c("PADI4","MPO","ELANE","TNF","CXCL8","GSDMD","TLR3")
NETosis.score <- Pathway_scoring(netosisgenes)
metadata_filtered$NETosis.score <- (NETosis.score)
```

As a sanity check, we score each sample based on a previously-defined
NETosis signature (Mukhopadhyay et al.) and check for a correlation with
our score.

``` r
gmt.file <- gmtPathways(paste0(prefix,"all_gene_sets.gmt"))
metadata_filtered$NETosis_Mukhopadhyay <- Pathway_scoring("NETOSIS_MUKHOPADHYAY")

metadata_temp <- metadata_filtered[metadata_filtered$COVID == "Positive",]

summary(lm(NETosis.score ~ NETosis_Mukhopadhyay, metadata_temp))
```

    ## 
    ## Call:
    ## lm(formula = NETosis.score ~ NETosis_Mukhopadhyay, data = metadata_temp)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.00050 -0.21420  0.00735  0.20974  1.29725 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          -0.13707    0.01409  -9.729   <2e-16 ***
    ## NETosis_Mukhopadhyay  1.55396    0.03062  50.752   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3216 on 621 degrees of freedom
    ## Multiple R-squared:  0.8057, Adjusted R-squared:  0.8054 
    ## F-statistic:  2576 on 1 and 621 DF,  p-value: < 2.2e-16

``` r
summary(lm(NETosis.score ~ NETosis_Mukhopadhyay, metadata_temp))$coefficients
```

    ##                        Estimate Std. Error   t value      Pr(>|t|)
    ## (Intercept)          -0.1370718 0.01408838 -9.729425  6.429878e-21
    ## NETosis_Mukhopadhyay  1.5539636 0.03061882 50.751910 3.931734e-223

**Figure S14B:**

``` r
plot(metadata_temp$NETosis.score, metadata_temp$NETosis_Mukhopadhyay, col = mapvalues(metadata_temp$severity.max, from = c("non-severe","severe"), to = c(my.cols[3],my.cols[1])), pch = 19, cex = 0.5, xlab = "NETosis Metagene Score", ylab = "Mukhopadhyay NETosis Metagene Score")
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

We then check how the NETosis metagene score varies across disease
severity, COVID status, disease acuity, and NMF cluster.

``` r
my.cols <- brewer.pal(3,"RdBu")

p1 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(Day), y = as.numeric(NETosis.score), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge()) + theme_bw() + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + theme(panel.grid = element_blank()) + xlab("Day") + ylab("NETosis Metagene Score") + stat_compare_means()
p1$labels$colour <- "Severity Max"
p2 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(cluster_neuhi), y = as.numeric(NETosis.score), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_bw() + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("Day") + ylab("NETosis Metagene Score") + stat_compare_means()
```

**Figure 4A:**

``` r
p1
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
p2
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
my.cols <- brewer.pal(3,"Set2")
p1 <- ggplot(metadata_filtered[metadata_filtered$Day == c("D0"),], aes(x = factor(COVID), y = as.numeric(NETosis.score), fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_bw() + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("COVID") + ylab("NETosis Metagene Score") + stat_compare_means() + coord_fixed(ratio = 1)
p2 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7") & metadata_filtered$severity.max == "severe",], aes(x = factor(Day), y = as.numeric(NETosis.score), fill = factor(Acuity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge()) + theme_bw() + scale_fill_manual(values = c("red","navy")) + theme(panel.grid = element_blank()) + xlab("Day") + ylab("NETosis Metagene Score") + stat_compare_means()
p2$labels$colour <- "Acuity Max"
p3 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0"),], aes(x = factor(cluster_neuhi), y = as.numeric(NETosis.score), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_bw() + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("Day") + ylab("NETosis Metagene Score") + stat_compare_means()
p4 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D3"),], aes(x = factor(cluster_neuhi), y = as.numeric(NETosis.score), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_bw() + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("Day") + ylab("NETosis Metagene Score") + stat_compare_means()
p5 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D7"),], aes(x = factor(cluster_neuhi), y = as.numeric(NETosis.score), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_bw() + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("Day") + ylab("NETosis Metagene Score") + stat_compare_means()
```

**Figure S14C:**

``` r
p1
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

**Figure S14D:**

``` r
p2
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

**Figure S14E:**

``` r
p3
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
p4
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

``` r
p5
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->

Many factors influencing NET release are post-transcriptional, so we
next searched for selected protein markers of NETosis in the matched
plasma samples. Here we read in the fully-processed selected Olink
proteins. The full code for pre-processing of the Olink data is included
in the code for the proteomics section, Figure 6. We compare the protein
levels across disease severity by day,

``` r
my.cols <- brewer.pal(3,"RdBu")
p1 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(Day), y = as.numeric(MPO_NPX), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.1) + theme_bw() + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("") + ylab("MPO NPX") + stat_compare_means() + coord_fixed(ratio = 1)
p2 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(Day), y = as.numeric(CXCL8_NPX), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.1) + theme_bw() + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("") + ylab("CXCL8 NPX") + stat_compare_means() + coord_fixed(ratio = .6)
p3 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(Day), y = as.numeric(TNF_NPX), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.1) + theme_bw() + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("") + ylab("TNF NPX") + stat_compare_means() + coord_fixed(ratio = .8)
p4 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(Day), y = as.numeric(PADI4_NPX), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.1) + theme_bw() + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("") + ylab("PADI4 NPX") + stat_compare_means() + coord_fixed(ratio = .8)
p5 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(Day), y = as.numeric(HGF_NPX), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.1) + theme_bw() + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("") + ylab("HGF NPX") + stat_compare_means() + coord_fixed(ratio = .65)
p6 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(Day), y = as.numeric(CD177_NPX), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.1) + theme_bw() + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("") + ylab("CD177 NPX") + stat_compare_means() + coord_fixed(ratio = .7)
```

**Figure 4B:**

``` r
cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol=3)
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
p1 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(cluster_neuhi), y = as.numeric(MPO_NPX), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.05) + theme_bw() + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("") + ylab("MPO NPX") + stat_compare_means()
p2 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(cluster_neuhi), y = as.numeric(CXCL8_NPX), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.05) + theme_bw() + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("") + ylab("CXCL8 NPX") + stat_compare_means()
p3 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(cluster_neuhi), y = as.numeric(TNF_NPX), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.05) + theme_bw() + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("") + ylab("TNF NPX") + stat_compare_means()
p4 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(cluster_neuhi), y = as.numeric(PADI4_NPX), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.05) + theme_bw() + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("") + ylab("PADI4 NPX") + stat_compare_means()
p5 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(cluster_neuhi), y = as.numeric(HGF_NPX), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.05) + theme_bw() + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("") + ylab("HGF NPX") + stat_compare_means()
p6 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),], aes(x = factor(cluster_neuhi), y = as.numeric(CD177_NPX), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.05) + theme_bw() + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + theme(panel.grid = element_blank(), legend.position = "none") + xlab("") + ylab("CD177 NPX") + stat_compare_means()
```

**Figure S14F:**

``` r
cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol=3)
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

We also looked at the concentration of cell-free DNA (cfDNA) as a marker
of NETosis. Recent studies have shown that the majority of cfDNA in the
plasma of severe COVID-19 patients originates from hematopoietic cells,
specifically neutrophils. cfDNA measurements were available from all
plasma samples regardless of if there was enough volume to isolate
neutrophils.

``` r
metadata_temp <- metadata_long[metadata_long$COVID == "1" & metadata_long$Day %in% c("D0","D3","D7"),]
my.cols <- brewer.pal(3, "RdBu")
p1 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(cell_free_dna), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("cfDNA Concentration (ng/uL)") + xlab("Day") + scale_y_log10(expand = c(0.01,0.01)) + scale_fill_manual(values = my.cols[c(3,1)]) + coord_fixed(ratio = 1.25) + stat_compare_means()
p1$labels$fill <- "Severity Max"

my.cols <- brewer.pal(3,"Set2")
p2 <- ggplot(metadata_long[metadata_long$Day == c("D0"),], aes(x = factor(COVID), y = as.numeric(cell_free_dna), fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_bw() + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + theme(panel.grid = element_blank(), legend.position = "none") + scale_y_log10(expand = c(0.01,0.01)) + xlab("COVID") + ylab("cfDNA Concentration (ng/uL)") + stat_compare_means() + coord_fixed(ratio = 2.5)

my.cols <- brewer.pal(9,"Blues")
metadata_temp <- metadata_long[metadata_long$Day == "D0" & metadata_long$COVID == "1",]
metadata_temp <- metadata_temp[complete.cases(metadata_temp$ANC.matched),]
p3 <- ggplot(metadata_temp, aes(x = factor(ANC.matched), y = as.numeric(cell_free_dna), fill = factor(ANC.matched))) + geom_boxplot(outlier.shape = NA) + theme_bw() + geom_jitter(height = 0, alpha = 0.3) + theme(legend.position = "none", panel.grid = element_blank()) + xlab("ANC Quintile") + ylab("cfDNA Concentration (ng/uL)") + stat_compare_means() + scale_y_log10(expand = c(0.01,0.01)) + coord_fixed(ylim = c(min((as.numeric(metadata_temp$cell_free_dna)),na.rm = TRUE),max((as.numeric(metadata_temp$cell_free_dna)),na.rm = TRUE)), ratio = 3.2) + scale_fill_manual(values = c(my.cols[c(1,3,5,7,9)]))
cor.test(y = as.numeric(metadata_temp$cell_free_dna), x = as.numeric(metadata_temp$ANC.matched), method = "kendall")
```

    ## 
    ##  Kendall's rank correlation tau
    ## 
    ## data:  as.numeric(metadata_temp$ANC.matched) and as.numeric(metadata_temp$cell_free_dna)
    ## z = 5.2324, p-value = 1.673e-07
    ## alternative hypothesis: true tau is not equal to 0
    ## sample estimates:
    ##       tau 
    ## 0.2350802

**Figure 4C:**

``` r
p1
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

**Figure 4D:**

``` r
p2
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

**Figure 4E:**

``` r
p3
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
metadata_temp <- metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7") & metadata_filtered$severity.max == "severe",]
my.cols <- brewer.pal(3, "RdBu")
p1 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(cell_free_dna), fill = factor(Acuity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("cfDNA Concentration (ng/uL)") + xlab("Day") + scale_y_log10(expand = c(0.01,0.01)) + scale_fill_manual(values = c("red","navy")) + coord_fixed(ratio = 1.25) + stat_compare_means()
p1$labels$fill <- "Acuity Max"

metadata_temp <- metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7") & metadata_filtered$cluster_neuhi %in% c(1,4),]
p2 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(cell_free_dna), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("cfDNA Concentration (ng/uL)") + xlab("Day") + scale_y_log10(expand = c(0.01,0.01)) + scale_fill_manual(values = c(orange,yellow)) + coord_fixed(ratio = 1.25) + stat_compare_means()
p2$labels$fill <- "Acuity Max"

my.cols <- brewer.pal(9,"Blues")
metadata_temp <- metadata_filtered[metadata_filtered$Day == "D3" & metadata_filtered$COVID == "Positive",]
metadata_temp <- metadata_temp[complete.cases(metadata_temp$ANC.matched),]
p3 <- ggplot(metadata_temp, aes(x = factor(ANC.matched), y = as.numeric(cell_free_dna), fill = factor(ANC.matched))) + geom_boxplot(outlier.shape = NA) + theme_bw() + geom_jitter(height = 0, alpha = 0.3) + theme(legend.position = "none", panel.grid = element_blank()) + xlab("ANC Quintile") + ylab("cfDNA Concentration (ng/uL)") + stat_compare_means() + scale_y_log10(expand = c(0.01,0.01)) + coord_fixed(ylim = c(min((as.numeric(metadata_temp$cell_free_dna)),na.rm = TRUE),max((as.numeric(metadata_temp$cell_free_dna)),na.rm = TRUE)), ratio = 3.2) + scale_fill_manual(values = c(my.cols[c(1,3,5,7,9)]))
cor.test(y = as.numeric(metadata_temp$cell_free_dna), x = as.numeric(metadata_temp$ANC.matched), method = "kendall")
```

    ## 
    ##  Kendall's rank correlation tau
    ## 
    ## data:  as.numeric(metadata_temp$ANC.matched) and as.numeric(metadata_temp$cell_free_dna)
    ## z = 5.6118, p-value = 2.002e-08
    ## alternative hypothesis: true tau is not equal to 0
    ## sample estimates:
    ##       tau 
    ## 0.3224663

``` r
my.cols <- brewer.pal(9,"Blues")
metadata_temp <- metadata_filtered[metadata_filtered$Day == "D7" & metadata_filtered$COVID == "Positive",]
metadata_temp <- metadata_temp[complete.cases(metadata_temp$ANC.matched),]
p4 <- ggplot(metadata_temp, aes(x = factor(ANC.matched), y = as.numeric(cell_free_dna), fill = factor(ANC.matched))) + geom_boxplot(outlier.shape = NA) + theme_bw() + geom_jitter(height = 0, alpha = 0.3) + theme(legend.position = "none", panel.grid = element_blank()) + xlab("ANC Quintile") + ylab("cfDNA Concentration (ng/uL)") + stat_compare_means() + scale_y_log10(expand = c(0.01,0.01)) + coord_fixed(ylim = c(min((as.numeric(metadata_temp$cell_free_dna)),na.rm = TRUE),max((as.numeric(metadata_temp$cell_free_dna)),na.rm = TRUE)), ratio = 3.2) + scale_fill_manual(values = c(my.cols[c(1,3,5,7,9)]))
cor.test(y = as.numeric(metadata_temp$cell_free_dna), x = as.numeric(metadata_temp$ANC.matched), method = "kendall")
```

    ## 
    ##  Kendall's rank correlation tau
    ## 
    ## data:  as.numeric(metadata_temp$ANC.matched) and as.numeric(metadata_temp$cell_free_dna)
    ## z = 4.5787, p-value = 4.679e-06
    ## alternative hypothesis: true tau is not equal to 0
    ## sample estimates:
    ##       tau 
    ## 0.3456522

**Figure S15A:**

``` r
p1
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

**Figure S15B:**

``` r
p2
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

**Figure S15C:**

``` r
p3
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
p4
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-25-2.png)<!-- -->

We also check for cfDNA differences in patients with detectable troponin
in their bloodstream within 72 hours of hospitalization.

``` r
metadata_temp <- metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),]
p1 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(cell_free_dna), fill = factor(Trop_72h))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("cfDNA Concentration (ng/uL)") + xlab("Day") + scale_y_log10(expand = c(0.01,0.01)) + scale_fill_manual(values = c("tan",vermillion)) + coord_fixed(ratio = 1.25) + stat_compare_means()
p1$labels$fill <- "Trop_72h"
```

**Figure S15D:**

``` r
p1
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

The next neutrophil effector function we explore is degranulation. There
are three types of neutrophil granules: azurophil, specific, and
tertiary. We will define a generic neutrophil degranulation score as
well as scores for each type of granule.

``` r
metadata_filtered$REACTOME_NEUTROPHIL_DEGRANULATION <- Pathway_scoring("REACTOME_NEUTROPHIL_DEGRANULATION")
metadata_filtered$GO_AZUROPHIL_GRANULE <- Pathway_scoring("GO_AZUROPHIL_GRANULE")
metadata_filtered$GO_SPECIFIC_GRANULE <- Pathway_scoring("GO_SPECIFIC_GRANULE")
metadata_filtered$GO_TERTIARY_GRANULE <- Pathway_scoring("GO_TERTIARY_GRANULE")
```

We check how the degranulation score breaks down by severity, NMF
cluster, and acuity.

``` r
metadata_temp <- metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),]
my.cols <- brewer.pal(3, "RdBu")
p1 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(REACTOME_NEUTROPHIL_DEGRANULATION), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("Neutrophil Degranulation Metagene Score") + xlab("Day") + scale_fill_manual(values = my.cols[c(3,1)]) + stat_compare_means()
p1$labels$fill <- "Severity Max"

p2 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$cluster_neuhi %in% c(1,4),], aes(x = factor(cluster_neuhi), y = as.numeric(REACTOME_NEUTROPHIL_DEGRANULATION), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.3) + theme_bw() + ylab("Neutrophil Degranulation Metagene Score") + xlab("") + scale_fill_manual(values = c(orange,yellow)) + stat_compare_means() + coord_fixed(ratio = 3)
p2$labels$fill <- "NMF Cluster"
```

**Figure 4F:**

``` r
p1
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
p2
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-30-2.png)<!-- -->

``` r
metadata_temp <- metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7") & metadata_filtered$severity.max == "severe",]
my.cols <- brewer.pal(3, "RdBu")
p1 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(REACTOME_NEUTROPHIL_DEGRANULATION), fill = factor(Acuity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("Neutrophil Degranulation Metagene Score") + xlab("Day") + scale_fill_manual(values = c("red","navy")) + stat_compare_means()
p1$labels$fill <- "Acuity Max"

p2 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive",], aes(x = factor(cluster_neuhi), y = as.numeric(REACTOME_NEUTROPHIL_DEGRANULATION), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.3) + theme_bw() + ylab("Neutrophil Degranulation Metagene Score") + xlab("") + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + stat_compare_means() + coord_fixed(ratio = 3)
p2$labels$fill <- "NMF Cluster"
```

**Figure S16A:**

``` r
p1
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

**Figure S16B:**

``` r
p2
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

Then we break down each of the three types of granule.

``` r
metadata_temp <- metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),]
my.cols <- brewer.pal(3, "RdBu")
p1 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(GO_AZUROPHIL_GRANULE), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("Azurophil Granule Metagene Score") + xlab("Day") + scale_fill_manual(values = my.cols[c(3,1)]) + stat_compare_means()
p1$labels$fill <- "Severity Max"

p2 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(GO_SPECIFIC_GRANULE), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("Specific Granule Metagene Score") + xlab("Day") + scale_fill_manual(values = my.cols[c(3,1)]) + stat_compare_means()
p2$labels$fill <- "Severity Max"

p3 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(GO_TERTIARY_GRANULE), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("Tertiary Granule Metagene Score") + xlab("Day") + scale_fill_manual(values = my.cols[c(3,1)]) + stat_compare_means()
p3$labels$fill <- "Severity Max"
```

**Figure 4G:**

``` r
p1
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

**Figure S16C:**

``` r
p2
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

**Figure S16D:**

``` r
p3
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

Finally, we investigate genes known to contribute to T cell suppression.
These include *ARG1*, *CD274*, *NECTIN2*, *IDO1*, and *SLC7A11*.

``` r
metadata_temp <- metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),]
my.cols <- brewer.pal(3, "RdBu")
p1 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(ARG1_logTPM), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("log2(TPM+1) Expression") + xlab("Day") + scale_fill_manual(values = my.cols[c(3,1)]) + stat_compare_means() + ggtitle("ARG1")
p1$labels$fill <- "Severity Max"

p2 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(CD274_logTPM), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("log2(TPM+1) Expression") + xlab("Day") + scale_fill_manual(values = my.cols[c(3,1)]) + stat_compare_means() + ggtitle("CD274")
p2$labels$fill <- "Severity Max"
```

**Figure 4H:**

``` r
p1
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

``` r
p2
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-39-2.png)<!-- -->

``` r
p1 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive",], aes(x = factor(cluster_neuhi), y = as.numeric(ARG1_logTPM), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.3) + theme_bw() + ylab("log2(TPM+1) Expression") + xlab("") + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + stat_compare_means() + ggtitle("ARG1")
p1$labels$fill <- "NMF Cluster"

p2 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive",], aes(x = factor(cluster_neuhi), y = as.numeric(CD274_logTPM), fill = factor(cluster_neuhi))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.3) + theme_bw() + ylab("log2(TPM+1) Expression") + xlab("") + scale_fill_manual(values = c(orange, skyblue, bluishgreen, yellow, blue, vermillion, "#c64dd1")) + stat_compare_means() + ggtitle("CD274")
p2$labels$fill <- "NMF Cluster"
```

**Figure S16E:**

``` r
p1
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

**Figure S16F:**

``` r
p2
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

``` r
metadata_temp <- metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),]
my.cols <- brewer.pal(3, "RdBu")
p1 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(NECTIN2_logTPM), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("log2(TPM+1) Expression") + xlab("Day") + scale_fill_manual(values = my.cols[c(3,1)]) + stat_compare_means() + ggtitle("NECTIN2")
p1$labels$fill <- "Severity Max"

p2 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(IDO1_logTPM), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("log2(TPM+1) Expression") + xlab("Day") + scale_fill_manual(values = my.cols[c(3,1)]) + stat_compare_means() + ggtitle("IDO1")
p2$labels$fill <- "Severity Max"

p3 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(SLC7A11_logTPM), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_bw() + ylab("log2(TPM+1) Expression") + xlab("Day") + scale_fill_manual(values = my.cols[c(3,1)]) + stat_compare_means() + ggtitle("SLC7A11")
p3$labels$fill <- "Severity Max"
```

**Figure S16G:**

``` r
p1
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

**Figure S16H:**

``` r
p2
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

**Figure S16I:**

``` r
p3
```

![](Figure4_S14-S16_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggpubr_0.4.0                fgsea_1.18.0               
    ##  [3] cowplot_1.1.1               openxlsx_4.2.4             
    ##  [5] DESeq2_1.32.0               SummarizedExperiment_1.22.0
    ##  [7] Biobase_2.52.0              MatrixGenerics_1.4.3       
    ##  [9] matrixStats_0.60.1          GenomicRanges_1.44.0       
    ## [11] GenomeInfoDb_1.28.2         IRanges_2.26.0             
    ## [13] S4Vectors_0.30.0            BiocGenerics_0.38.0        
    ## [15] dplyr_1.0.7                 plyr_1.8.6                 
    ## [17] RColorBrewer_1.1-2          ggrepel_0.9.1              
    ## [19] ggplot2_3.3.5               knitr_1.33                 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-7           bit64_4.0.5            httr_1.4.2            
    ##  [4] tools_4.1.1            backports_1.2.1        utf8_1.2.2            
    ##  [7] R6_2.5.1               DBI_1.1.1              colorspace_2.0-2      
    ## [10] withr_2.4.2            tidyselect_1.1.1       gridExtra_2.3         
    ## [13] bit_4.0.4              curl_4.3.2             compiler_4.1.1        
    ## [16] DelayedArray_0.18.0    labeling_0.4.2         scales_1.1.1          
    ## [19] genefilter_1.74.0      stringr_1.4.0          digest_0.6.27         
    ## [22] foreign_0.8-81         rmarkdown_2.10         rio_0.5.27            
    ## [25] XVector_0.32.0         pkgconfig_2.0.3        htmltools_0.5.2       
    ## [28] highr_0.9              fastmap_1.1.0          readxl_1.3.1          
    ## [31] rlang_0.4.11           RSQLite_2.2.8          farver_2.1.0          
    ## [34] generics_0.1.0         BiocParallel_1.26.2    zip_2.2.0             
    ## [37] car_3.0-11             RCurl_1.98-1.4         magrittr_2.0.1        
    ## [40] GenomeInfoDbData_1.2.6 Matrix_1.3-4           Rcpp_1.0.7            
    ## [43] munsell_0.5.0          fansi_0.5.0            abind_1.4-5           
    ## [46] lifecycle_1.0.0        stringi_1.7.4          yaml_2.2.1            
    ## [49] carData_3.0-4          zlibbioc_1.38.0        grid_4.1.1            
    ## [52] blob_1.2.2             forcats_0.5.1          crayon_1.4.1          
    ## [55] lattice_0.20-44        Biostrings_2.60.2      haven_2.4.3           
    ## [58] splines_4.1.1          annotate_1.70.0        hms_1.1.0             
    ## [61] KEGGREST_1.32.0        locfit_1.5-9.4         pillar_1.6.2          
    ## [64] ggsignif_0.6.2         codetools_0.2-18       geneplotter_1.70.0    
    ## [67] fastmatch_1.1-3        XML_3.99-0.7           glue_1.4.2            
    ## [70] evaluate_0.14          data.table_1.14.0      png_0.1-7             
    ## [73] vctrs_0.3.8            cellranger_1.1.0       gtable_0.3.0          
    ## [76] purrr_0.3.4            tidyr_1.1.3            assertthat_0.2.1      
    ## [79] cachem_1.0.6           xfun_0.25              xtable_1.8-4          
    ## [82] broom_0.7.9            rstatix_0.7.0          survival_3.2-13       
    ## [85] tibble_3.1.4           AnnotationDbi_1.54.1   memoise_2.0.0         
    ## [88] ellipsis_0.3.2
