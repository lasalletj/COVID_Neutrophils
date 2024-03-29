Figure1_S1
================
Tom LaSalle

This document contains all the code necessary to generate the plots for
Figure 1 and related figure S1. Plots are subsequently edited in Adobe
Illustrator to produce the final figures.

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
library(corrplot)
library(ggmosaic)
library(cowplot)
library(Seurat)
library(heatmap3)
library(ggpubr)
library(umap)
library(fgsea)
```

Import the metadata and keep only data for which neutrophil enrichment
was performed:

``` r
prefix <- "~/Downloads/COVID19_Neutrophil_Code/" #Adapt as necessary
metadata_bypatient <- read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 2)
qc_data <- read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 7)
metadata_bypatient <- metadata_bypatient[which(metadata_bypatient$Public.ID %in% qc_data$Public.ID),]
genomic_signatures <- read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 10)
metadata_long <- read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 4)
```

Due to an institutional IRB-approved waiver of informed consent, all
clinical data are reported in quintiles, and not all clinical
characteristics are publicly available. Correlations were calculated
between the following clinical variables: Age, Sex, BMI, Heart
condition, Lung condition, Kidney condition, Diabetes, Hypertension,
Immunocompromised, Symptom duration, Respiratory symptoms, Fever
symptoms, GI Symptoms, Acuity (D0, D3, D7, Max within D7-D28, D28, Max
Overall), ANC (D0, D3, D7), ALC (D0, D3, D7), AMC (D0, D3, D7),
Creatinine (D0, D3, D7), CRP (D0, D3, D7), D-dimer (D0, D3, D7), LDH
(D0, D3, D7), Troponin detected within 72 hours, Intubation status, Race
(White, Black, Asian, Other), Hispanic ethnicity, CXR infiltrates, and
Smoking history. Multiple hypothesis correction was performed to
determine significant correlations, and then those parameters with
significant correlations were plotted in a correlation heatmap. Though
it is not possible to show all of the data, the code to generate the
figure is below and is written as if the additional columns were present
in the metadata file:

``` r
# # Reverse factor levels of Acuity such that positive correlations with Acuity indicate correlations with worse disease severity
# metadata_bypatient$Acuity.0 <- mapvalues(metadata_bypatient$Acuity.0, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
# metadata_bypatient$Acuity.3 <- mapvalues(metadata_bypatient$Acuity.3, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
# metadata_bypatient$Acuity.7 <- mapvalues(metadata_bypatient$Acuity.7, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
# metadata_bypatient$Acuity.28 <- mapvalues(metadata_bypatient$Acuity.28, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
# metadata_bypatient$Acuity.7.to.28 <- mapvalues(metadata_bypatient$Acuity.7.to.28, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
# metadata_bypatient$Acuity.max <- mapvalues(metadata_bypatient$Acuity.max, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
# 
# rownames(metadata_bypatient) <- metadata_bypatient$Public.ID
# 
# metadata_bypatient <- metadata_bypatient[,-which(colnames(metadata_bypatient) %in% c("Study_ID","Public.ID","D0_draw","D3_draw","D7_draw","DE_draw","D0_Neu_Sample","D3_Neu_Sample","D7_Neu_Sample","DE_Neu_Sample"))]
# metadata_bypatient <- metadata_bypatient[metadata_bypatient$COVID == "1",]
# metadata_bypatient <- metadata_bypatient[,-which(colnames(metadata_bypatient) %in% c("COVID"))]
# 
# cormatrix <- matrix(0L, nrow = ncol(metadata_bypatient), ncol = ncol(metadata_bypatient))
# pmatrix <- matrix(0L, nrow = ncol(metadata_bypatient), ncol = ncol(metadata_bypatient))
# rownames(cormatrix) <- rownames(pmatrix) <- colnames(cormatrix) <- colnames(pmatrix) <- colnames(metadata_bypatient)
# 
# for (i in 1:nrow(cormatrix)){
#   for (j in 1:nrow(cormatrix)){
#     stats <- cor.test(x = as.numeric(metadata_bypatient[,i]), y = as.numeric(metadata_bypatient[,j]), use = "pairwise.complete.obs", method = "kendall")
#     cormatrix[i,j] <- stats$estimate
#     pmatrix[i,j] <- stats$p.value
#   }
# }

cormatrix <- as.matrix(read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 5, rowNames = TRUE))
pmatrix <- as.matrix(read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 6, rowNames = TRUE))

ptable <- matrix(0L, nrow = nrow(pmatrix)^2, ncol = 3)
colnames(ptable) <- c("Row","Column","pvalue")
ptable[,1] <- rep(colnames(pmatrix), each = nrow(pmatrix))
ptable[,2] <- rep(colnames(pmatrix), nrow(pmatrix))
for (i in 1:nrow(ptable)){
  ptable[i,3] <- pmatrix[rownames(pmatrix) %in% ptable[i,1],colnames(pmatrix) %in% ptable[i,2]]
}
ptable <- as.data.frame(ptable[complete.cases(ptable[,3]),])
ptable$fdr <- p.adjust(ptable[,3], method = "fdr")
ptable <- ptable[ptable[,4] < 0.05,]
ptable <- ptable[ptable[,1] %in% c("abs_neut_0_cat","abs_neut_3_cat","abs_neut_7_cat") | ptable[,2] %in% c("abs_neut_0_cat","abs_neut_3_cat","abs_neut_7_cat"),]

cormatrix <- cormatrix[rownames(cormatrix) %in% unique(c(ptable[,1],ptable[,2])),colnames(cormatrix) %in% unique(c(ptable[,1],ptable[,2]))]
pmatrix <- pmatrix[rownames(pmatrix) %in% unique(c(ptable[,1],ptable[,2])),colnames(pmatrix) %in% unique(c(ptable[,1],ptable[,2]))]

rownames(cormatrix) <- colnames(cormatrix) <- rownames(pmatrix) <- colnames(pmatrix) <- mapvalues(rownames(cormatrix), from = c("BMI.cat","abs_mono_0_cat","abs_mono_3_cat","abs_mono_7_cat","Age.cat","creat_0_cat","creat_3_cat","creat_7_cat","ldh_0_cat","ldh_3_cat","ldh_7_cat","crp_0_cat","crp_3_cat","crp_7_cat","Acuity.28","Acuity.7","Acuity.3","Acuity.0","Acuity.max","abs_neut_0_cat","abs_neut_3_cat","abs_neut_7_cat","ddimer_0_cat","ddimer_3_cat","ddimer_7_cat"), to = c("BMI","AMC, D0", "AMC, D3","AMC, D7","Age","Creatinine, D0","Creatinine, D3","Creatinine, D7","LDH, D0","LDH, D3","LDH, D7","CRP, D0","CRP, D3","CRP, D7","Acuity, D28","Acuity, D7","Acuity, D3","Acuity, D0","Acuity, Max","ANC, D0","ANC, D3","ANC, D7","D-dimer, D0","D-dimer, D3","D-dimer, D7"))
```

**Figure 1B:**

``` r
corrplot(cormatrix, method = "square", type = "lower", order = "hclust", hclust.method = "ward.D", tl.col="black", col=colorRampPalette(c("blue","white","red"))(200), p.mat = pmatrix, sig.level = 0.05, insig = "blank")
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Next we highlight the correlations between ANC on each individual day
and Acuity Max, using Kendall’s tau correlation:

``` r
metadata_long <- metadata_long[which(metadata_long$Public.ID %in% qc_data$Public.ID),]
metadata_temp <- metadata_long[complete.cases(metadata_long$ANC.matched) & metadata_long$COVID == 1,]

metadata_temp$Acuity.max <- factor(metadata_temp$Acuity.max, levels = c("5","4","3","2","1"))
corvals0 <- cor.test(x = as.numeric(metadata_temp[metadata_temp$Day == "D0",]$Acuity.max), y = as.numeric(metadata_temp[metadata_temp$Day == "D0",]$ANC.matched), use = "pairwise.complete.obs", method = "kendall")
corvals3 <- cor.test(x = as.numeric(metadata_temp[metadata_temp$Day == "D3",]$Acuity.max), y = as.numeric(metadata_temp[metadata_temp$Day == "D3",]$ANC.matched), use = "pairwise.complete.obs", method = "kendall")
corvals7 <- cor.test(x = as.numeric(metadata_temp[metadata_temp$Day == "D7",]$Acuity.max), y = as.numeric(metadata_temp[metadata_temp$Day == "D7",]$ANC.matched), use = "pairwise.complete.obs", method = "kendall")

metadata_temp$ANC.matched <- factor(metadata_temp$ANC.matched, levels = c("5","4","3","2","1"))

my.cols <- brewer.pal(9, "RdBu")
p1 <- ggplot(data = metadata_temp[metadata_temp$Day == "D0",]) + geom_mosaic(aes(x = product(ANC.matched, Acuity.max), fill = ANC.matched), offset = 0.005) + coord_fixed(ratio = 1) + theme_bw() + scale_fill_manual(values = (my.cols[(c(1,2,3,7,9))])) + ggtitle("Day 0") + annotate(geom="text", x=0.5, y=.8, label=paste0("τ =",round(cormatrix[rownames(cormatrix) == "ANC, D0",colnames(cormatrix) == "Acuity, Max"],3),"\npadj =",round(pmatrix[rownames(pmatrix) == "ANC, D0",colnames(pmatrix) == "Acuity, Max"],12)), color="black") + theme(legend.key.size=unit(4,'mm'), legend.text=element_text(size=5), legend.title=element_text(size=5), axis.title = element_text(size = 7))
p2 <- ggplot(data = metadata_temp[metadata_temp$Day == "D3",]) + geom_mosaic(aes(x = product(ANC.matched, Acuity.max), fill = ANC.matched), offset = 0.005) + coord_fixed(ratio = 1) + theme_bw() + scale_fill_manual(values = (my.cols[(c(1,2,3,7,9))])) + ggtitle("Day 3") + annotate(geom="text", x=0.5, y=.8, label=paste0("τ =",round(cormatrix[rownames(cormatrix) == "ANC, D3",colnames(cormatrix) == "Acuity, Max"],3),"\npadj =",round(pmatrix[rownames(pmatrix) == "ANC, D3",colnames(pmatrix) == "Acuity, Max"],11)), color="black") + theme(legend.key.size=unit(4,'mm'), legend.text=element_text(size=5), legend.title=element_text(size=5), axis.title = element_text(size = 7))
p3 <- ggplot(data = metadata_temp[metadata_temp$Day == "D7",]) + geom_mosaic(aes(x = product(ANC.matched, Acuity.max), fill = ANC.matched), offset = 0.005) + coord_fixed(ratio = 1) + theme_bw() + scale_fill_manual(values = (my.cols[(c(1,2,3,7,9))])) + ggtitle("Day 7") + annotate(geom="text", x=0.5, y=.8, label=paste0("τ =",round(cormatrix[rownames(cormatrix) == "ANC, D7",colnames(cormatrix) == "Acuity, Max"],3),"\npadj =",round(pmatrix[rownames(pmatrix) == "ANC, D7",colnames(pmatrix) == "Acuity, Max"],11)), color="black") + theme(legend.key.size=unit(4,'mm'), legend.text=element_text(size=5), legend.title=element_text(size=5), axis.title = element_text(size = 7))
```

**Figure 1C:**

``` r
plot_grid(p1,p2,p3,ncol=3)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Now we begin to look at our transcriptomic data. We import the count and
TPM data from RSEM, as well as the Ensembl ID to gene symbol
conversions, and log-transform the TPM values:

``` r
Count <- read.table(gzfile(paste0(prefix,"Neutrophil_RNAseq_Count_Matrix.txt.gz")),sep="\t")
colnames(Count) <- Count[1,]
Count <- Count[-1,]
Count <- Count[,-2]
rownames(Count) <- Count[,1]
nams <- Count[,1]
Count <- Count[,-1]
Count <- as.data.frame(apply(Count,2,as.numeric))
rownames(Count) <- nams
TPM <- read.table(gzfile(paste0(prefix,"Neutrophil_RNAseq_TPM_Matrix.txt.gz")),sep="\t")
colnames(TPM) <- TPM[1,]
TPM <- TPM[-1,]
TPM <- TPM[,-2]
rownames(TPM) <- TPM[,1]
nams <- TPM[,1]
TPM <- TPM[,-1]
TPM <- as.data.frame(apply(TPM,2,as.numeric))
rownames(TPM) <- nams
```

``` r
genepc <- read.delim(paste0(prefix,"Ensembl_to_Symbol.txt"))
logTPM <- log2(TPM + 1)
metadata_long <- merge(metadata_long, qc_data)
rownames(metadata_long) <- metadata_long$Public.Sample.ID
```

To get an overall sense of the data we start with a PCA of the
protein-coding genes:

``` r
proteincoding <- genepc$Gene.stable.ID[genepc$Gene.type == "protein_coding"]
logPCG <- logTPM[which(rownames(logTPM) %in% proteincoding),]
logPCG <- logPCG[which(rowSums(logPCG) > 0),]
pcaResults <- prcomp(t(logPCG), center = T, scale. = T)
PCs <- data.frame(pcaResults$x[,1:2])
metadata_long$PC1 <- PCs$PC1
metadata_long$PC2 <- PCs$PC2
```

The main contributions to PC1 are Exonic rate and Median 3’ bias:

``` r
myPalette <- colorRampPalette((brewer.pal(9, "RdYlBu")))
sc <- scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(0,70))
p1 <- ggplot(data = metadata_long, aes(x = -1*PC1, y = PC2, color = Exonic.Rate*100)) + geom_point(size = 0.9) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank())  + ggtitle("PCA") + sc + coord_fixed(ratio = 2.2)
p1$labels$colour <- "Exonic Rate"

myPalette <- colorRampPalette((brewer.pal(9, "RdYlGn")))
sc <- scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(min(metadata_long$Median.3..bias),1))
p2 <- ggplot(data = metadata_long, aes(x = -1*PC1, y = PC2, color = Median.3..bias)) + geom_point(size = 0.9) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank())  + ggtitle("PCA") + sc + coord_fixed(ratio = 2.2)
p2$labels$colour <- "Median 3' Bias"
```

**Figure Not Included:**

``` r
plot_grid(p1,p2)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

We apply quality control filtration to the data and re-examine the PCA:

``` r
#Calculate mitochondrial percentages
mitorows <- grepl("^MT-",genepc$Gene.name)
mitoids <- genepc$Gene.stable.ID[mitorows]
mitocounts <- Count[rownames(Count) %in% mitoids,]
mitosums <- colSums(mitocounts)
totalcounts <- colSums(Count)
percent.mt <- mitosums/totalcounts*100
metadata_long$percent.mt <- percent.mt

#Perform quality control filtration based on RNASeQC parameters
metadata_filtered <- metadata_long[metadata_long$percent.mt < 20 & metadata_long$Genes.Detected > 10000 & metadata_long$Median.Exon.CV < 1 & metadata_long$Exon.CV.MAD < 0.75 & metadata_long$Exonic.Rate*100 > 25 & metadata_long$Median.3..bias < 0.9,]

logPCG_filtered <- logPCG[,colnames(logPCG) %in% rownames(metadata_filtered)]
logTPM_filtered <- logTPM[,colnames(logTPM) %in% rownames(metadata_filtered)]
TPM_filtered <- TPM[,colnames(TPM) %in% rownames(metadata_filtered)]
Count_filtered <- Count[,colnames(Count) %in% rownames(metadata_filtered)]

tf <- rowSums(TPM_filtered > 0.1) > ncol(TPM_filtered)*.2
TPM_filtered <- TPM_filtered[tf,]
Count_filtered <- Count_filtered[tf,]
logTPM_filtered <- logTPM_filtered[tf,]
tf <- rowSums(Count_filtered >= 6) > ncol(Count_filtered)*.2
TPM_filtered <- TPM_filtered[tf,]
Count_filtered <- Count_filtered[tf,]
logTPM_filtered <- logTPM_filtered[tf,]

#Repeat PCA
logPCG_filtered <- logPCG_filtered[which(rowSums(logPCG_filtered) > 0),]
pcaResults <- prcomp(t(logPCG_filtered), center = T, scale. = T)
PCs <- data.frame(pcaResults$x[,1:2])
metadata_filtered$PC1 <- PCs$PC1
metadata_filtered$PC2 <- PCs$PC2

myPalette <- colorRampPalette((brewer.pal(9, "RdYlBu")))
sc <- scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(0,70))
p1 <- ggplot(data = metadata_filtered, aes(x = -1*PC1, y = -1*PC2, color = Exonic.Rate*100)) + geom_point(size = 0.9) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank())  + ggtitle("PCA") + sc + coord_fixed(ratio = 1.3)
p1$labels$colour <- "Exonic Rate"

myPalette <- colorRampPalette((brewer.pal(9, "RdYlGn")))
sc <- scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(min(metadata_filtered$Median.3..bias),1))
p2 <- ggplot(data = metadata_filtered, aes(x = -1*PC1, y = -1*PC2, color = Median.3..bias)) + geom_point(size = 0.9) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=12),axis.title=element_text(size=14), plot.title = element_text(lineheight=.8, face="bold", size = 16), axis.ticks.y=element_blank())  + ggtitle("PCA") + sc + coord_fixed(ratio = 1.3)
p2$labels$colour <- "Median 3' Bias"
```

**Figure Not Included:**

``` r
plot_grid(p1,p2)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

We check that quality control filtration has not introduced any biases
in terms of patient demographics.

``` r
metadata_long$QC <- rownames(metadata_long) %in% rownames(metadata_filtered)
metadata_long$QC <- mapvalues(metadata_long$QC, from = c("FALSE","TRUE"), to = c("FAIL","PASS"))

histmaker <- function(parameter){
  metadata_temp <- metadata_long[,colnames(metadata_long) %in% c(parameter,"QC")]
  metadata_temp[,1] <- factor(metadata_temp[,1])
  histtable <- as.data.frame(matrix(0L, nrow = nlevels(metadata_temp[,1])*2, ncol = 3))
  colnames(histtable) <- c("Param","QC","Freq")
  histtable[,1] <- rep(levels(metadata_temp[,1]), each = 2)
  histtable[,2] <- rep(c("FAIL","PASS"),nlevels(metadata_temp[,1]))
  for (i in 1:nrow(histtable)){
     histtable$Freq[i] <- sum(metadata_temp[,1] == histtable$Param[i] & metadata_temp$QC == histtable$QC[i])
  }
  return(histtable)
}

parameter <- "COVID"
fishertest <- fisher.test(table(metadata_long[,colnames(metadata_long) == "QC"],metadata_long[,colnames(metadata_long) == parameter]))
histtable <- histmaker(parameter)
p1 <- ggplot(histtable, aes(x = Param, y = Freq, fill = QC)) + geom_bar(position = "dodge",stat = "identity") + theme_bw() + ylab("Count") + xlab(parameter) + coord_fixed(ratio = .0033) + scale_fill_manual(values = c("darkred","forestgreen")) + annotate("text", x=2, y=max(histtable$Freq), label= paste0("Fisher: ",round(fishertest$p.value, digits = 3)), colour = "black")

parameter <- "Acuity.max"
fishertest <- fisher.test(table(metadata_long[,colnames(metadata_long) == "QC"],metadata_long[,colnames(metadata_long) == parameter]))
histtable <- histmaker(parameter)
p2 <- ggplot(histtable, aes(x = Param, y = Freq, fill = QC)) + geom_bar(position = "dodge",stat = "identity") + theme_bw() + ylab("Count") + xlab(parameter) + coord_fixed(ratio = .018) + scale_fill_manual(values = c("darkred","forestgreen")) + annotate("text", x=2, y=max(histtable$Freq), label= paste0("Fisher: ",round(fishertest$p.value, digits = 3)), colour = "black")

parameter <- "severity.max"
fishertest <- fisher.test(table(metadata_long[,colnames(metadata_long) == "QC"],metadata_long[,colnames(metadata_long) == parameter]))
histtable <- histmaker(parameter)
p3 <- ggplot(histtable, aes(x = Param, y = Freq, fill = QC)) + geom_bar(position = "dodge",stat = "identity") + theme_bw() + ylab("Count") + xlab(parameter) + coord_fixed(ratio = .008) + scale_fill_manual(values = c("darkred","forestgreen")) + annotate("text", x=2, y=max(histtable$Freq), label= paste0("Fisher: ",round(fishertest$p.value, digits = 3)), colour = "black")

parameter <- "Age"
fishertest <- fisher.test(table(metadata_long[,colnames(metadata_long) == "QC"],metadata_long[,colnames(metadata_long) == parameter]))
histtable <- histmaker(parameter)
p4 <- ggplot(histtable, aes(x = Param, y = Freq, fill = QC)) + geom_bar(position = "dodge",stat = "identity") + theme_bw() + ylab("Count") + xlab(parameter) + coord_fixed(ratio = .024) + scale_fill_manual(values = c("darkred","forestgreen")) + annotate("text", x=2, y=max(histtable$Freq), label= paste0("Fisher: ",round(fishertest$p.value, digits = 3)), colour = "black")

#parameter <- "sex"
#fishertest <- fisher.test(table(metadata_long[,colnames(metadata_long) == "QC"],metadata_long[,colnames(metadata_long) == parameter]))
#histtable <- histmaker(parameter)
#p5 <- ggplot(histtable, aes(x = Param, y = Freq, fill = QC)) + geom_bar(position = "dodge",stat = "identity") + theme_bw() + ylab("Count") + xlab(parameter) + coord_fixed(ratio = .008) + scale_fill_manual(values = c("darkred","forestgreen")) + annotate("text", x=2, y=max(histtable$Freq), label= paste0("Fisher: ",round(fishertest$p.value, digits = 3)), colour = "black")
p5 <- ""

#parameter <- "race"
#fishertest <- fisher.test(table(metadata_long[,colnames(metadata_long) == "QC"],metadata_long[,colnames(metadata_long) == parameter]))
#histtable <- histmaker(parameter)
#p6 <- ggplot(histtable, aes(x = Param, y = Freq, fill = QC)) + geom_bar(position = "dodge",stat = "identity") + theme_bw() + ylab("Count") + xlab(parameter) + coord_fixed(ratio = .008) + scale_fill_manual(values = c("darkred","forestgreen")) + annotate("text", x=2, y=max(histtable$Freq), label= paste0("Fisher: ",round(fishertest$p.value, digits = 3)), colour = "black")
p6 <- ""

#parameter <- "ethnicity"
#fishertest <- fisher.test(table(metadata_long[,colnames(metadata_long) == "QC"],metadata_long[,colnames(metadata_long) == parameter]))
#histtable <- histmaker(parameter)
#p7 <- ggplot(histtable, aes(x = Param, y = Freq, fill = QC)) + geom_bar(position = "dodge",stat = "identity") + theme_bw() + ylab("Count") + xlab(parameter) + coord_fixed(ratio = .008) + scale_fill_manual(values = c("darkred","forestgreen")) + annotate("text", x=2, y=max(histtable$Freq), label= paste0("Fisher: ",round(fishertest$p.value, digits = 3)), colour = "black")
p7 <- ""

parameter <- "BMI"
fishertest <- fisher.test(table(metadata_long[,colnames(metadata_long) == "QC"],metadata_long[,colnames(metadata_long) == parameter]))
histtable <- histmaker(parameter)
p8 <- ggplot(histtable, aes(x = Param, y = Freq, fill = QC)) + geom_bar(position = "dodge",stat = "identity") + theme_bw() + ylab("Count") + xlab(parameter) + coord_fixed(ratio = .024) + scale_fill_manual(values = c("darkred","forestgreen")) + annotate("text", x=2, y=max(histtable$Freq), label= paste0("Fisher: ",round(fishertest$p.value, digits = 3)), colour = "black")
```

**Figures Not Included:**

``` r
p1
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
p2
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
p3
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->

``` r
p4
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

``` r
p8
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-15-5.png)<!-- -->

We do observe that more samples are discarded from later time points:

``` r
parameter <- "Day"
fishertest <- fisher.test(table(metadata_long[,colnames(metadata_long) == "QC"],metadata_long[,colnames(metadata_long) == parameter]))
histtable <- histmaker(parameter)
p1 <- ggplot(histtable, aes(x = Param, y = Freq, fill = QC)) + geom_bar(position = "dodge",stat = "identity") + theme_bw() + ylab("Count") + xlab(parameter) + coord_fixed(ratio = .014) + scale_fill_manual(values = c("darkred","forestgreen")) + annotate("text", x=2, y=max(histtable$Freq), label= paste0("Fisher: p = ",round(fishertest$p.value, digits = 8)), colour = "black")
```

**Figure Not Included:**

``` r
p1
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

This may be partially attributed to longer sample processing times.

**Figure Not Included:**

``` r
wilcoxtest <- wilcox.test(Processing.Time ~ QC, metadata_long)
ggplot(metadata_long, aes(x = factor(QC), y = as.numeric(Processing.Time), fill = factor(QC))) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_bw() + ylab("Processing Time (h)") + scale_fill_manual(values = c("darkred","forestgreen")) + coord_fixed(ratio = .08) + annotate("text", x=1, y=23, label= paste0("Wilcoxon: p = ",round(wilcoxtest$p.value, digits = 5)), colour = "black")
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

However, there are no significant differences in processing time between
samples that pass or fail within the same day.

**Figure Not Included:**

``` r
ggplot(metadata_long, aes(x = factor(Day), y = Processing.Time, fill = QC)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.height = 0, jitter.width = 0.2), size = 0.5, alpha = .5) + theme_bw() + scale_fill_manual(values = c("darkred","forestgreen")) + xlab("Day") + ylab("Processing Time (h)") + coord_fixed(ratio = 0.21)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

Next we want to get an estimate of the sample purity since we performed
bulk RNA-seq of enriched neutrophils and 100% purity is not guaranteed
(see flow cytometry data from Figure S1A-B). In order to do so, we need
a reference single-cell RNA-seq dataset for deconvolution using
CIBERSORTx. We will use the Bonn Cohort 2 from Schulte-Schrepping et
al. We start by importing the data and collapsing the cluster labels
into major hematopoeitic cell lineages.

``` r
seuratwb <- readRDS(paste0(prefix,"seurat_COVID19_freshWB-PBMC_cohort2_rhapsody_jonas_FG_2020-08-18.rds"))

barcodecluster <- as.character(seuratwb$RNA_snn_res.0.8)
temp <- seuratwb$RNA_snn_res.0.8
names(barcodecluster) <- names(temp)
barcodecluster <- cbind(rownames(barcodecluster),barcodecluster)
barcodecluster <- as.data.frame(barcodecluster)
barcodecluster$barcode <- as.character(rownames(barcodecluster))
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("0","1","4","8")] <- "Mature_Neutrophil"
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("12","15")] <- "Immature_Neutrophil"
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("3","10","14","18")] <- "Monocyte"
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("2","5","6","13","16")] <- "T_NK"
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("7","22")] <- "B"
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("19")] <- "Plasmablast"
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("9","11","17","20","21","23","24")] <- "Other"
seuratwb <- AddMetaData(seuratwb, metadata = subset(barcodecluster, select = c("barcode")), col.name = "Cell.Type")
```

**Figure S1C:**

``` r
DimPlot(seuratwb, reduction = "umap", group.by = "Cell.Type", label = FALSE) + scale_color_manual(values = c("grey60","tomato","forestgreen","skyblue","grey90","black","slateblue3")) + theme(axis.line.x = element_blank(), axis.line.y = element_blank(), axis.title.x = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

For CIBERSORTx, we will generate a cell type signature matrix based on
“pseudobulk” samples from a single cell type. We will add up the counts
from all the cells of a given lineage for each patient and treat them as
bulk samples.

``` r
pseudobulk <- matrix(0L, nrow = nrow(seuratwb@assays$RNA@counts), ncol = (length(levels(factor(seuratwb$donor))))*6)
rownames(pseudobulk) <- (seuratwb@assays$RNA@data@Dimnames[[1]])
idx <- (levels(factor(seuratwb$donor)))
pseudobulk <- as.data.frame(pseudobulk)
for (i in 0:(ncol(pseudobulk)/6 - 1)){
  k <- i*6+1
  colnames(pseudobulk)[k] <- paste(as.character(idx[i+1]),"_Mature_Neu",sep = "")
  colnames(pseudobulk)[k+1] <- paste(as.character(idx[i+1]),"_Immature_Neu",sep = "")
  colnames(pseudobulk)[k+2] <- paste(as.character(idx[i+1]),"_Mono",sep = "")
  colnames(pseudobulk)[k+3] <- paste(as.character(idx[i+1]),"_TNK",sep = "")
  colnames(pseudobulk)[k+4] <- paste(as.character(idx[i+1]),"_B",sep = "")
  colnames(pseudobulk)[k+5] <- paste(as.character(idx[i+1]),"_Plasmablast",sep = "")
}

#Patient BN-31 does not have all the cell types so they are excluded
for (i in 0:17){
  k <- i*6+1
  temp <- subset(seuratwb, subset = donor == as.character(idx[i+1]) & Cell.Type == "Mature_Neutrophil")
  pseudobulk[,k] <- rowSums(temp@assays$RNA@counts)
  temp <- subset(seuratwb, subset = donor == as.character(idx[i+1]) & Cell.Type == "Immature_Neutrophil")
  pseudobulk[,k+1] <- rowSums(temp@assays$RNA@counts)
  temp <- subset(seuratwb, subset = donor == as.character(idx[i+1]) & Cell.Type == "Monocyte")
  pseudobulk[,k+2] <- rowSums(temp@assays$RNA@counts)
  temp <- subset(seuratwb, subset = donor == as.character(idx[i+1]) & Cell.Type == "T_NK")
  pseudobulk[,k+3] <- rowSums(temp@assays$RNA@counts)
  temp <- subset(seuratwb, subset = donor == as.character(idx[i+1]) & Cell.Type == "B")
  pseudobulk[,k+4] <- rowSums(temp@assays$RNA@counts)
  temp <- subset(seuratwb, subset = donor == as.character(idx[i+1]) & Cell.Type == "Plasmablast")
  pseudobulk[,k+5] <- rowSums(temp@assays$RNA@counts)
}
for (i in 19:(ncol(pseudobulk)/6 - 1)){
  k <- i*6+1
  temp <- subset(seuratwb, subset = donor == as.character(idx[i+1]) & Cell.Type == "Mature_Neutrophil")
  pseudobulk[,k] <- rowSums(temp@assays$RNA@counts)
  temp <- subset(seuratwb, subset = donor == as.character(idx[i+1]) & Cell.Type == "Immature_Neutrophil")
  pseudobulk[,k+1] <- rowSums(temp@assays$RNA@counts)
  temp <- subset(seuratwb, subset = donor == as.character(idx[i+1]) & Cell.Type == "Monocyte")
  pseudobulk[,k+2] <- rowSums(temp@assays$RNA@counts)
  temp <- subset(seuratwb, subset = donor == as.character(idx[i+1]) & Cell.Type == "T_NK")
  pseudobulk[,k+3] <- rowSums(temp@assays$RNA@counts)
  temp <- subset(seuratwb, subset = donor == as.character(idx[i+1]) & Cell.Type == "B")
  pseudobulk[,k+4] <- rowSums(temp@assays$RNA@counts)
  temp <- subset(seuratwb, subset = donor == as.character(idx[i+1]) & Cell.Type == "Plasmablast")
  pseudobulk[,k+5] <- rowSums(temp@assays$RNA@counts)
}
pseudobulk <- pseudobulk[rowSums(pseudobulk) > 0,]
pseudobulk <- pseudobulk[-grep("\\.",rownames(pseudobulk)),]
pseudobulk <- pseudobulk[-grep("^MT",rownames(pseudobulk)),]
```

Then we convert the Ensembl gene IDs to symbols for the neutrophil data
and filter the pseudobulk matrix to include only the overlapping genes
between the two datasets:

``` r
genepctemp <- genepc[genepc$Gene.stable.ID %in% rownames(Count_filtered),]
Count_temp <- Count
for (i in 1:nrow(Count_temp)){
  if (length(which(genepctemp$Gene.stable.ID == rownames(Count_temp)[i])) == 1){
    if (!(genepctemp$Gene.name[which(genepctemp$Gene.stable.ID == rownames(Count_temp)[i])] %in% rownames(Count_temp))){
      rownames(Count_temp)[i] <- genepctemp$Gene.name[which(genepctemp$Gene.stable.ID == rownames(Count_temp)[i])]
    }
  }
}
pseudobulk_temp <- pseudobulk[rownames(pseudobulk) %in% rownames(Count_temp),]
```

The matrix that was used to generate the signature matrix is included in
the Github.

Now we are able to run CIBERSORTx Generate Signature Matrix, setting
limits of 50 to 100 genes per cell type and filtering for only
hematopoietic genes. The CIBERSORTx cell type signature matrix:

**Figure S1D:**

``` r
matrix <- read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 13, rowNames = TRUE)
par(mar=c(1,1,1,1))
heatmap3(matrix, Colv = NA, cexRow = 0.2)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

Using this cell type signature matrix, we generate estimated cell type
fractions in CIBERSORTx using default parameters. An overview of the
cell type proportions:

``` r
rownames(genomic_signatures) <- genomic_signatures$Public.Sample.ID
cibersort_temp <- genomic_signatures[,colnames(genomic_signatures) %in% c("Mature_Neutrophil","Immature_Neutrophil","Monocyte","T_NK","B","Plasmablast","Neutrophil_total")]
cibersort_temp <- cibersort_temp[rownames(cibersort_temp) %in% rownames(metadata_filtered),]
cibersort_sorted <- cibersort_temp[order(cibersort_temp$Neutrophil_total),]
cibersort_sorted <- cibersort_sorted[,-7]
data_percentage <- t(cibersort_sorted*100)
coul <- c("forestgreen","tomato","skyblue","slateblue3","gray","black")
```

**Figure S1E:**

``` r
barplot(data_percentage, col=coul , border=NA, xlab=NA, axisnames = FALSE)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

Search for differences in CIBERSORTx estimated cell fractions across
COVID status for all samples:

``` r
metadata_filtered <- merge(metadata_filtered, genomic_signatures)
rownames(metadata_filtered) <- metadata_filtered$Public.Sample.ID
metadata_filtered$COVID <- mapvalues(metadata_filtered$COVID, from = c(0,1), to = c("Negative","Positive"))

my.cols <- brewer.pal(3, "Set2")

wilcoxtest <- wilcox.test(as.numeric(Neutrophil_total)*100 ~ factor(COVID), metadata_filtered)
p1 <- ggplot(metadata_filtered, aes(x = factor(COVID), y = as.numeric(Neutrophil_total)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 3)), colour = "black") + ggtitle("Total")

wilcoxtest <- wilcox.test(as.numeric(Mature_Neutrophil)*100 ~ factor(COVID), metadata_filtered)
p2 <- ggplot(metadata_filtered, aes(x = factor(COVID), y = as.numeric(Mature_Neutrophil)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 8)), colour = "black") + ggtitle("Mature")

wilcoxtest <- wilcox.test(as.numeric(Immature_Neutrophil)*100 ~ factor(COVID), metadata_filtered)
p3 <- ggplot(metadata_filtered, aes(x = factor(COVID), y = as.numeric(Immature_Neutrophil)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 8)), colour = "black") + ggtitle("Immature")

wilcoxtest <- wilcox.test(as.numeric(Monocyte)*100 ~ factor(COVID), metadata_filtered)
p4 <- ggplot(metadata_filtered, aes(x = factor(COVID), y = as.numeric(Monocyte)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 8)), colour = "black") + ggtitle("Monocyte")

wilcoxtest <- wilcox.test(as.numeric(T_NK)*100 ~ factor(COVID), metadata_filtered)
p5 <- ggplot(metadata_filtered, aes(x = factor(COVID), y = as.numeric(T_NK)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 4)), colour = "black") + ggtitle("T/NK")

wilcoxtest <- wilcox.test(as.numeric(B)*100 ~ factor(COVID), metadata_filtered)
p6 <- ggplot(metadata_filtered, aes(x = factor(COVID), y = as.numeric(B)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 3)), colour = "black") + ggtitle("B")

wilcoxtest <- wilcox.test(as.numeric(Plasmablast)*100 ~ factor(COVID), metadata_filtered)
p7 <- ggplot(metadata_filtered, aes(x = factor(COVID), y = as.numeric(Plasmablast)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 10)), colour = "black") + ggtitle("Plasmablast")
```

**Figure Not Included:**

``` r
plot_grid(p1,p2,p3,p4,p5,p6,p7,ncol = 7)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

Search for differences in CIBERSORTx estimated cell fractions across
disease severity for all samples regardless of COVID status:

``` r
metadata_filtered <- merge(metadata_filtered, genomic_signatures)
metadata_temp <- metadata_filtered[metadata_filtered$severity.max %in% c("severe","non-severe"),]

my.cols <- brewer.pal(3, "RdBu")

wilcoxtest <- wilcox.test(as.numeric(Neutrophil_total)*100 ~ factor(severity.max), metadata_temp)
p1 <- ggplot(metadata_temp, aes(x = factor(severity.max), y = as.numeric(Neutrophil_total)*100, fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 18)), colour = "black") + ggtitle("Total")

wilcoxtest <- wilcox.test(as.numeric(Mature_Neutrophil)*100 ~ factor(severity.max), metadata_temp)
p2 <- ggplot(metadata_temp, aes(x = factor(severity.max), y = as.numeric(Mature_Neutrophil)*100, fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 3)), colour = "black") + ggtitle("Mature")

wilcoxtest <- wilcox.test(as.numeric(Immature_Neutrophil)*100 ~ factor(severity.max), metadata_temp)
p3 <- ggplot(metadata_temp, aes(x = factor(severity.max), y = as.numeric(Immature_Neutrophil)*100, fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 8)), colour = "black") + ggtitle("Immature")

wilcoxtest <- wilcox.test(as.numeric(Monocyte)*100 ~ factor(severity.max), metadata_temp)
p4 <- ggplot(metadata_temp, aes(x = factor(severity.max), y = as.numeric(Monocyte)*100, fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 5)), colour = "black") + ggtitle("Monocyte")

wilcoxtest <- wilcox.test(as.numeric(T_NK)*100 ~ factor(severity.max), metadata_temp)
p5 <- ggplot(metadata_temp, aes(x = factor(severity.max), y = as.numeric(T_NK)*100, fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 24)), colour = "black") + ggtitle("T/NK")

wilcoxtest <- wilcox.test(as.numeric(B)*100 ~ factor(severity.max), metadata_temp)
p6 <- ggplot(metadata_temp, aes(x = factor(severity.max), y = as.numeric(B)*100, fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 3)), colour = "black") + ggtitle("B")

wilcoxtest <- wilcox.test(as.numeric(Plasmablast)*100 ~ factor(severity.max), metadata_temp)
p7 <- ggplot(metadata_temp, aes(x = factor(severity.max), y = as.numeric(Plasmablast)*100, fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 4)), colour = "black") + ggtitle("Plasmablast")
```

**Figure Not Included:**

``` r
plot_grid(p1,p2,p3,p4,p5,p6,p7,ncol = 7)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

We repeat this analysis on Day 0 samples only as there were no
COVID-19-negative samples collected after Day 0:

``` r
metadata_temp <- metadata_filtered[metadata_filtered$Day == "D0",]

my.cols <- brewer.pal(3, "Set2")

wilcoxtest <- wilcox.test(as.numeric(Neutrophil_total)*100 ~ factor(COVID), metadata_temp)
p1 <- ggplot(metadata_temp, aes(x = factor(COVID), y = as.numeric(Neutrophil_total)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 3)), colour = "black") + ggtitle("Total")

wilcoxtest <- wilcox.test(as.numeric(Mature_Neutrophil)*100 ~ factor(COVID), metadata_temp)
p2 <- ggplot(metadata_temp, aes(x = factor(COVID), y = as.numeric(Mature_Neutrophil)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 8)), colour = "black") + ggtitle("Mature")

wilcoxtest <- wilcox.test(as.numeric(Immature_Neutrophil)*100 ~ factor(COVID), metadata_temp)
p3 <- ggplot(metadata_temp, aes(x = factor(COVID), y = as.numeric(Immature_Neutrophil)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 8)), colour = "black") + ggtitle("Immature")

wilcoxtest <- wilcox.test(as.numeric(Monocyte)*100 ~ factor(COVID), metadata_temp)
p4 <- ggplot(metadata_temp, aes(x = factor(COVID), y = as.numeric(Monocyte)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 8)), colour = "black") + ggtitle("Monocyte")

wilcoxtest <- wilcox.test(as.numeric(T_NK)*100 ~ factor(COVID), metadata_temp)
p5 <- ggplot(metadata_temp, aes(x = factor(COVID), y = as.numeric(T_NK)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 4)), colour = "black") + ggtitle("T/NK")

wilcoxtest <- wilcox.test(as.numeric(B)*100 ~ factor(COVID), metadata_temp)
p6 <- ggplot(metadata_temp, aes(x = factor(COVID), y = as.numeric(B)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 3)), colour = "black") + ggtitle("B")

wilcoxtest <- wilcox.test(as.numeric(Plasmablast)*100 ~ factor(COVID), metadata_temp)
p7 <- ggplot(metadata_temp, aes(x = factor(COVID), y = as.numeric(Plasmablast)*100, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(wilcoxtest$p.value, digits = 15)), colour = "black") + ggtitle("Plasmablast")
```

**Figure Not Included:**

``` r
plot_grid(p1,p2,p3,p4,p5,p6,p7,ncol = 7)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

Next we want to understand determinants of sample purity. We can observe
correlations between ANC and CIBERSORTx Total neutrophil fraction, as
well as ALC and CIBERSORTx T/NK fraction.

``` r
metadata_temp <- metadata_filtered[metadata_filtered$COVID == "Positive",]
metadata_temp <- metadata_temp[complete.cases(metadata_temp$AMC.matched) & complete.cases(metadata_temp$ALC.matched),]

my.cols <- brewer.pal(3, "Set2")
p1 <- ggplot(metadata_temp, aes(x = factor(ANC.matched), y = (Neutrophil_total)*100)) + geom_boxplot(outlier.shape = NA, fill = "blue3") + geom_jitter(height = 0, alpha = 0.3) + theme_bw() + stat_compare_means() + ylab("CIBERSORTx Total Neutrophil (%)") + xlab("ANC Category") + scale_fill_manual(values = my.cols[c(1,2)]) + scale_y_continuous(limits = c(0,max(metadata_temp$Neutrophil_total)*100)) + coord_fixed(ratio = 1/14) + theme(legend.position = "none")

p2 <- ggplot(metadata_temp, aes(x = factor(ALC.matched), y = (T_NK)*100)) + geom_boxplot(outlier.shape = NA, fill = "mediumpurple4") + geom_jitter(height = 0, alpha = 0.3) + theme_bw() + stat_compare_means() + ylab("CIBERSORTx T/NK (%)") + xlab("ALC Category") + scale_fill_manual(values = my.cols[c(1,2)]) + scale_y_continuous(limits = c(0,max(metadata_temp$T_NK)*100)) + coord_fixed(ratio = 1/11) + theme(legend.position = "none")

p3 <- ggplot(metadata_temp, aes(x = factor(AMC.matched), y = (Monocyte)*100)) + geom_boxplot(outlier.shape = NA, fill = "skyblue") + geom_jitter(height = 0, alpha = 0.3) + theme_bw() + stat_compare_means() + ylab("CIBERSORTx Monocyte (%)") + xlab("AMC Category") + scale_fill_manual(values = my.cols[c(1,2)]) + scale_y_continuous(limits = c(0,max(metadata_temp$T_NK)*100)) + coord_fixed(ratio = 1/11) + theme(legend.position = "none")
```

**Figure Not Included:**

``` r
plot_grid(p1,p2,p3,ncol=3)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

We see that lower purity samples tend to have lower ANC and higher ALC.
The converse is also true.

``` r
metadata_lo <- metadata_filtered[order(metadata_filtered$Neutrophil_total),]
metadata_lo <- metadata_lo[1:round(nrow(metadata_lo)*.25),]
comparisondf <- matrix(0L, nrow = 25, ncol = 3)
colnames(comparisondf) <- c("Var1","Var2","Freq")
comparisondf[,1] <- c(rep(5,5),rep(4,5),rep(3,5),rep(2,5),rep(1,5))
comparisondf[,2] <- c(rep(1:5,5))
for (i in 1:nrow(comparisondf)){
  comparisondf[i,3] <- sum(metadata_lo$ANC.matched == comparisondf[i,1] & metadata_lo$ALC.matched == comparisondf[i,2], na.rm = TRUE)
}
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.grid = element_blank())   
}
my.cols <- brewer.pal(7, "Greens")
p1 <- ggplot(as.data.frame(comparisondf), aes(x = factor(Var1), y = Freq, fill = factor(Var2))) + geom_bar(position="stack", stat="identity") + scale_fill_manual(values = (my.cols[3:7])) + xlab("ANC Category") + ylab("Count") + theme_bw() + coord_fixed(ratio = .1)
p1$labels$fill <- "ALC Quintile"

metadata_hi <- metadata_filtered[order(metadata_filtered$Neutrophil_total),]
metadata_hi <- metadata_hi[round(nrow(metadata_hi)*.75):nrow(metadata_hi),]
comparisondf <- matrix(0L, nrow = 25, ncol = 3)
colnames(comparisondf) <- c("Var1","Var2","Freq")
comparisondf[,1] <- c(rep(5,5),rep(4,5),rep(3,5),rep(2,5),rep(1,5))
comparisondf[,2] <- c(rep(1:5,5))
for (i in 1:nrow(comparisondf)){
  comparisondf[i,3] <- sum(metadata_hi$ANC.matched == comparisondf[i,1] & metadata_hi$ALC.matched == comparisondf[i,2], na.rm = TRUE)
}
my.cols <- brewer.pal(7, "Greens")
p2 <- ggplot(as.data.frame(comparisondf), aes(x = factor(Var1), y = Freq, fill = factor(Var2))) + geom_bar(position="stack", stat="identity") + scale_fill_manual(values = (my.cols[3:7])) + xlab("ANC Category") + ylab("Count") + theme_bw() + coord_fixed(ratio = .1)
p2$labels$fill <- "ALC Quintile"

metadata_filtered$Purity <- -1*as.numeric(rownames(metadata_filtered) %in% rownames(metadata_lo)) + as.numeric(rownames(metadata_filtered) %in% rownames(metadata_hi))
metadata_filtered$Purity <- mapvalues(metadata_filtered$Purity, from = c(-1,0,1), to = c("lo","mid","hi"))
metadata_temp <- metadata_filtered[metadata_filtered$Purity %in% c("lo","hi"),]
fisher.test(table(metadata_temp$ANC.matched, metadata_temp$Purity))
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  table(metadata_temp$ANC.matched, metadata_temp$Purity)
    ## p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
fisher.test(table(metadata_temp$ALC.matched, metadata_temp$Purity))
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  table(metadata_temp$ALC.matched, metadata_temp$Purity)
    ## p-value = 0.01286
    ## alternative hypothesis: two.sided

**Figure Not Included:**

``` r
p1
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

``` r
p2
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-36-2.png)<!-- -->

Next, we view the CIBERSORTx estimated neutrophil fractions over time in
COVID-19-positive patients:

``` r
my.cols <- brewer.pal(3, "Set2")
metadata_temp <- metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0","D3","D7"),]

kwtest <- kruskal.test(g = metadata_temp$Day, x = metadata_temp$Neutrophil_total)
p1 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(Neutrophil_total)*100, fill = factor(Day))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2],my.cols[3])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(kwtest$p.value, digits = 15)), colour = "black") + ggtitle("Total")

kwtest <- kruskal.test(g = metadata_temp$Day, x = metadata_temp$Mature_Neutrophil)
p2 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(Mature_Neutrophil)*100, fill = factor(Day))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2],my.cols[3])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(kwtest$p.value, digits = 37)), colour = "black") + ggtitle("Mature")

kwtest <- kruskal.test(g = metadata_temp$Day, x = metadata_temp$Immature_Neutrophil)
p3 <- ggplot(metadata_temp, aes(x = factor(Day), y = as.numeric(Immature_Neutrophil)*100, fill = factor(Day))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.08) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .04) + scale_fill_manual(values = c(my.cols[1],my.cols[2],my.cols[3])) + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + annotate("text", x=1, y=95, label= paste0("p=",round(kwtest$p.value, digits = 34)), colour = "black") + ggtitle("Immature")
```

**Figure 1D:**

``` r
plot_grid(p1,p2,p3,ncol=3)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

We also check if CIBERSORTx estimated neutrophil-related cell type
fractions are associated with any clinical parameters.

``` r
my.cols <- brewer.pal(5,"Reds")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$Creatinine.matched),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p1 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(Creatinine.matched), y = Neutrophil_total, fill = factor(Creatinine.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p2 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(Creatinine.matched), y = Neutrophil_total, fill = factor(Creatinine.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p3 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(Creatinine.matched), y = Neutrophil_total, fill = factor(Creatinine.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)

my.cols <- brewer.pal(5,"Greens")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$CRP.matched),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p4 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(CRP.matched), y = Neutrophil_total, fill = factor(CRP.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p5 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(CRP.matched), y = Neutrophil_total, fill = factor(CRP.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p6 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(CRP.matched), y = Neutrophil_total, fill = factor(CRP.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)

my.cols <- brewer.pal(5,"Blues")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$Ddimer.matched),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p7 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(Ddimer.matched), y = Neutrophil_total, fill = factor(Ddimer.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p8 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(Ddimer.matched), y = Neutrophil_total, fill = factor(Ddimer.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p9 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(Ddimer.matched), y = Neutrophil_total, fill = factor(Ddimer.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)

my.cols <- brewer.pal(5,"Purples")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$LDH.matched),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p10 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(LDH.matched), y = Neutrophil_total, fill = factor(LDH.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p11 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(LDH.matched), y = Neutrophil_total, fill = factor(LDH.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p12 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(LDH.matched), y = Neutrophil_total, fill = factor(LDH.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)

# my.cols <- brewer.pal(10,"Paired")
# metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$intubated),]
# metadata_temp <- metadata_temp[metadata_temp$COVID == "1",]
# p13 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(intubated), y = Neutrophil_total, fill = factor(intubated))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = (my.cols[5:6]))
# p14 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(intubated), y = Neutrophil_total, fill = factor(intubated))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = (my.cols[5:6]))
# p15 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(intubated), y = Neutrophil_total, fill = factor(intubated))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = (my.cols[5:6]))

my.cols <- brewer.pal(10,"Paired")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$CXR.infiltrates),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p16 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(CXR.infiltrates), y = Neutrophil_total, fill = factor(CXR.infiltrates))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[7:8])
p17 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(CXR.infiltrates), y = Neutrophil_total, fill = factor(CXR.infiltrates))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[7:8])
p18 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(CXR.infiltrates), y = Neutrophil_total, fill = factor(CXR.infiltrates))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[7:8])

my.cols <- brewer.pal(12,"Paired")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$Trop_72h),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p19 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(Trop_72h), y = Neutrophil_total, fill = factor(Trop_72h))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[11:12])
p20 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(Trop_72h), y = Neutrophil_total, fill = factor(Trop_72h))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[11:12])
p21 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(Trop_72h), y = Neutrophil_total, fill = factor(Trop_72h))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[11:12])
```

**Figure Not Included (Total Neutrophil):**

``` r
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=3)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

``` r
plot_grid(#p13,p14,p15,
          p16,p17,p18,p19,p20,p21,ncol=3)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-40-2.png)<!-- -->

``` r
my.cols <- brewer.pal(5,"Reds")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$Creatinine.matched),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p1 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(Creatinine.matched), y = Mature_Neutrophil, fill = factor(Creatinine.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p2 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(Creatinine.matched), y = Mature_Neutrophil, fill = factor(Creatinine.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p3 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(Creatinine.matched), y = Mature_Neutrophil, fill = factor(Creatinine.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)

my.cols <- brewer.pal(5,"Greens")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$CRP.matched),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p4 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(CRP.matched), y = Mature_Neutrophil, fill = factor(CRP.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p5 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(CRP.matched), y = Mature_Neutrophil, fill = factor(CRP.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p6 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(CRP.matched), y = Mature_Neutrophil, fill = factor(CRP.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)

my.cols <- brewer.pal(5,"Blues")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$Ddimer.matched),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p7 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(Ddimer.matched), y = Mature_Neutrophil, fill = factor(Ddimer.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p8 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(Ddimer.matched), y = Mature_Neutrophil, fill = factor(Ddimer.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p9 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(Ddimer.matched), y = Mature_Neutrophil, fill = factor(Ddimer.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)

my.cols <- brewer.pal(5,"Purples")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$LDH.matched),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p10 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(LDH.matched), y = Mature_Neutrophil, fill = factor(LDH.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p11 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(LDH.matched), y = Mature_Neutrophil, fill = factor(LDH.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p12 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(LDH.matched), y = Mature_Neutrophil, fill = factor(LDH.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)

# my.cols <- brewer.pal(10,"Paired")
# metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$intubated),]
# metadata_temp <- metadata_temp[metadata_temp$COVID == "1",]
# p13 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(intubated), y = Mature_Neutrophil, fill = factor(intubated))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = (my.cols[5:6]))
# p14 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(intubated), y = Mature_Neutrophil, fill = factor(intubated))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = (my.cols[5:6]))
# p15 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(intubated), y = Mature_Neutrophil, fill = factor(intubated))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = (my.cols[5:6]))

my.cols <- brewer.pal(10,"Paired")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$CXR.infiltrates),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p16 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(CXR.infiltrates), y = Mature_Neutrophil, fill = factor(CXR.infiltrates))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[7:8])
p17 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(CXR.infiltrates), y = Mature_Neutrophil, fill = factor(CXR.infiltrates))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[7:8])
p18 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(CXR.infiltrates), y = Mature_Neutrophil, fill = factor(CXR.infiltrates))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[7:8])

my.cols <- brewer.pal(12,"Paired")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$Trop_72h),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p19 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(Trop_72h), y = Mature_Neutrophil, fill = factor(Trop_72h))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[11:12])
p20 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(Trop_72h), y = Mature_Neutrophil, fill = factor(Trop_72h))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[11:12])
p21 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(Trop_72h), y = Mature_Neutrophil, fill = factor(Trop_72h))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[11:12])
```

**Figure Not Included (Mature Neutrophil):**

``` r
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=3)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

``` r
plot_grid(#p13,p14,p15,
          p16,p17,p18,p19,p20,p21,ncol=3)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-42-2.png)<!-- -->

``` r
my.cols <- brewer.pal(5,"Reds")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$Creatinine.matched),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p1 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(Creatinine.matched), y = Immature_Neutrophil, fill = factor(Creatinine.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p2 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(Creatinine.matched), y = Immature_Neutrophil, fill = factor(Creatinine.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p3 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(Creatinine.matched), y = Immature_Neutrophil, fill = factor(Creatinine.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)

my.cols <- brewer.pal(5,"Greens")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$CRP.matched),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p4 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(CRP.matched), y = Immature_Neutrophil, fill = factor(CRP.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p5 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(CRP.matched), y = Immature_Neutrophil, fill = factor(CRP.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p6 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(CRP.matched), y = Immature_Neutrophil, fill = factor(CRP.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)

my.cols <- brewer.pal(5,"Blues")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$Ddimer.matched),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p7 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(Ddimer.matched), y = Immature_Neutrophil, fill = factor(Ddimer.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p8 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(Ddimer.matched), y = Immature_Neutrophil, fill = factor(Ddimer.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p9 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(Ddimer.matched), y = Immature_Neutrophil, fill = factor(Ddimer.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)

my.cols <- brewer.pal(5,"Purples")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$LDH.matched),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p10 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(LDH.matched), y = Immature_Neutrophil, fill = factor(LDH.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p11 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(LDH.matched), y = Immature_Neutrophil, fill = factor(LDH.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)
p12 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(LDH.matched), y = Immature_Neutrophil, fill = factor(LDH.matched))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 3) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols)

# my.cols <- brewer.pal(10,"Paired")
# metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$intubated),]
# metadata_temp <- metadata_temp[metadata_temp$COVID == "1",]
# p13 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(intubated), y = Immature_Neutrophil, fill = factor(intubated))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = (my.cols[5:6]))
# p14 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(intubated), y = Immature_Neutrophil, fill = factor(intubated))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = (my.cols[5:6]))
# p15 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(intubated), y = Immature_Neutrophil, fill = factor(intubated))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = (my.cols[5:6]))

my.cols <- brewer.pal(10,"Paired")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$CXR.infiltrates),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p16 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(CXR.infiltrates), y = Immature_Neutrophil, fill = factor(CXR.infiltrates))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[7:8])
p17 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(CXR.infiltrates), y = Immature_Neutrophil, fill = factor(CXR.infiltrates))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[7:8])
p18 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(CXR.infiltrates), y = Immature_Neutrophil, fill = factor(CXR.infiltrates))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[7:8])

my.cols <- brewer.pal(12,"Paired")
metadata_temp <- metadata_filtered[complete.cases(metadata_filtered$Trop_72h),]
metadata_temp <- metadata_temp[metadata_temp$COVID == "Positive",]
p19 <- ggplot(metadata_temp[metadata_temp$Day == "D0",], aes(x = factor(Trop_72h), y = Immature_Neutrophil, fill = factor(Trop_72h))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[11:12])
p20 <- ggplot(metadata_temp[metadata_temp$Day == "D3",], aes(x = factor(Trop_72h), y = Immature_Neutrophil, fill = factor(Trop_72h))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[11:12])
p21 <- ggplot(metadata_temp[metadata_temp$Day == "D7",], aes(x = factor(Trop_72h), y = Immature_Neutrophil, fill = factor(Trop_72h))) + geom_violin() + theme_bw() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + coord_fixed(ratio = 1.2) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + stat_compare_means(label.x = 1, label.y = .1) + scale_fill_manual(values = my.cols[11:12])
```

**Figure Not Included (Immature Neutrophil):**

``` r
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=3)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

``` r
plot_grid(#p13,p14,p15,
          p16,p17,p18,p19,p20,p21,ncol=3)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-44-2.png)<!-- -->

We next generate UMAPs for sample visualization in Figure 1E.

``` r
set.seed(10101)
#tpm.umap <- umap(t(logTPM_filtered))
#metadata_filtered$umap1 <- tpm.umap$layout[,1]
#metadata_filtered$umap2 <- tpm.umap$layout[,2]

metadata_filtered$severity.covid <- matrix(0L, nrow = nrow(metadata_filtered), ncol = 1)
for (i in 1:nrow(metadata_filtered)){
  if (metadata_filtered$COVID[i] == "Positive"){
    metadata_filtered$severity.covid[i] <- metadata_filtered$severity.max[i]
  }
  if (metadata_filtered$COVID[i] == "Negative"){
    if (metadata_filtered$severity.max[i] == "H"){
      metadata_filtered$severity.covid[i] <- "H"
    }
    else {
      metadata_filtered$severity.covid[i] <- "non_COVID"
    }
  }
}

my.cols <- brewer.pal(3, "Set2")
p1 <- ggplot(metadata_filtered, aes(x = umap1, y = umap2)) + theme_bw() + geom_point(aes(colour = factor(COVID)),size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) + scale_color_manual(values = c(my.cols[1], my.cols[2]))
p1$labels$colour <- "COVID"

my.cols.2 <- brewer.pal(3, "RdBu")
p2 <- ggplot(metadata_filtered, aes(x = umap1, y = umap2)) + theme_bw() + geom_point(aes(colour = factor(severity.covid)),size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) + scale_color_manual(values = c("grey", my.cols[1], my.cols.2[3],  my.cols[2]))
p2$labels$colour <- "Severity"

a <- "Mature_Neutrophil"
myPalette <- colorRampPalette((brewer.pal(9, "RdYlBu")))
sc <- scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(min(metadata_filtered[,colnames(metadata_filtered) == a])*100,max(metadata_filtered[,colnames(metadata_filtered) == a])*100)) 
p3 <- ggplot(metadata_filtered, aes(x = umap1, y = umap2, colour = (as.numeric(Mature_Neutrophil)*100))) + theme_bw() + geom_point(size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) + sc
p3$labels$colour <- "Mature Neu"

a <- "Immature_Neutrophil"
myPalette <- colorRampPalette((brewer.pal(9, "RdYlBu")))
sc <- scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(min(metadata_filtered[,colnames(metadata_filtered) == a])*100,max(metadata_filtered[,colnames(metadata_filtered) == a])*100)) 
p4 <- ggplot(metadata_filtered, aes(x = umap1, y = umap2, colour = (as.numeric(Immature_Neutrophil)*100))) + theme_bw() + geom_point(size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) + sc
p4$labels$colour <- "Immature Neu"
```

**Figure 1E:**

``` r
plot_grid(p1,p2,p3,p4,ncol=2)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

We perform differential expression to test for genes and pathways
associated with COVID-19 infection. We use the helper script
Neutrophil_DESeq2.R to import the function. As we will see, if we do not
make a correction for plasmablast contamination, we will have many
highly expressed immunoglobulin genes in the COVID+ side.

``` r
source(paste0(prefix,"Neutrophil_DESeq2.R"))
```

``` r
rownames(metadata_filtered) <- metadata_filtered$Public.Sample.ID
DESeq2_list <- Neutrophil_DESeq2(counts = Count_filtered, mdata = metadata_filtered, day = "D0")
dds <- DESeqDataSetFromMatrix(countData = DESeq2_list$Count_select, colData = DESeq2_list$coldata, design = ~ COVID)
dds <- DESeq(dds)

res <- as.data.frame(results(dds, name="COVID_Positive_vs_Negative"))
filenam <- "Day0_COVIDpos_vs_COVIDneg_uncorrected"
temp <- genepc[which(genepc$Gene.stable.ID %in% rownames(res)),]
res$symbol <- matrix(0L, nrow = nrow(res))
for (i in 1:nrow(res)){
  if (rownames(res)[i] %in% temp$Gene.stable.ID){
    res$symbol[i] <- temp$Gene.name[which(rownames(res)[i] == temp$Gene.stable.ID)]
  } else {
    res$symbol[i] <- rownames(res)[i]
  }
}
res$rank <- sign(res$log2FoldChange)*(-1)*log10(res$pvalue)
res <- res[complete.cases(res),]
res_sig <- res[res$padj < 0.05,]
#write.table(res,paste0(prefix,"DESeq2/",filenam,".txt"),sep = "\t")

resordered <- res[order(res$rank),]

log2fc <- as.numeric(resordered$log2FoldChange)
log10p <- as.numeric(-1*log10(resordered$pvalue))
padj <- resordered$padj
rank <- resordered$rank
symbol <- resordered$symbol
combo <- cbind(log2fc,log10p,padj,symbol,rank)
colnames(combo) <- c("log2fc","log10p","padj","symbol","rank")
rownames(combo) <- rownames(resordered)
combo <- as.data.frame(combo)
combo$rank <- as.numeric(combo$rank)
combo$log10p <- as.numeric(combo$log10p)
combo$log2fc <- as.numeric(combo$log2fc)
combo$significance <- as.numeric(combo$padj < 0.05)
combo$significance <- as.factor(combo$significance)

IGgenes <- c("IGKC","IGKV4-1","IGKV5-2","IGKV6-21","IGKV3D-20","IGKV3D-11","IGLV4-69","IGLV8-61","IGLV4-60","IGLV6-57","IGLV10-54","IGLV1-51","IGLV1-47","IGLV7-46","IGLV5-45","IGLV1-44","IGLV7-43","IGLV1-40","IGLV1-36","IGLV3-27","IGLV3-25","IGLV2-23","IGLV3-21","IGLV3-19","IGLV2-18","IGLV3-16","IGLV2-14","IGLV2-11","IGLV3-10","IGLV3-9","IGLV3-1","IGLC1","IGLC2","IGLC3","IGLC7","IGHA2","IGHG4","IGHG2","IGHA1","IGHG1","IGHG3","IGHD","IGHM","IGHV6-1","IGHV1-2","IGHV1-3","IGHV2-5","IGHV3-7","IGHV3-11","IGHV3-13","IGHV3-15","IGHV1-18","IGHV3-20","IGHV3-21","IGHV3-23","IGHV1-24","IGHV2-26","IGHV3-33","IGHV4-34","IGHV4-39","IGHV1-46","IGHV3-48","IGHV3-49","IGHV5-51","IGHV3-53","IGHV1-58","IGHV4-61","IGHV3-66","IGHV1-69","IGHV2-70D","IGHV3-73","IGLV9-49","IGHV3-64","IGKV3D-15","IGHV4-59","IGHV3-74","IGHV3-72","IGHV4-31","IGHV3-43","IGKV2D-30","IGKV1-6","IGKV3-20","IGKV1D-33","IGKV1-17","IGKV1-8","IGKV1-16","IGKV1D-16","IGKV2-24","IGKV3-11","IGKV1-9","IGKV1-33","IGKV1-39","IGKV2D-28","IGKV1D-17","IGKV2-30","IGKV2D-29","IGKV1-12","IGKV1-5","IGKV2-28","IGKV3-15","IGKV1-27","IGKV2D-40","IGKV1D-39","IGLVI-70","IGKV2-29","IGLL5","IGHV3-30","IGHV2-70","IGKV1D-13","IGHV4-4","IGLV2-8","IGHV1-69D","IGHV7-4-1","IGHV3-64D","IGHV5-10-1")
combo$color <- 0
combo$color[combo$symbol %in% IGgenes] <- 1

combo$labels <- 0
combo$labels[combo$symbol %in% c("MZB1","JCHAIN")] <- 1
combo$labels <- as.factor(combo$labels)

options(ggrepel.max.overlaps = Inf)
my.cols <- brewer.pal(3,"Set2")
plot1 <- ggplot(combo, aes(x = log2fc, y = log10p)) + geom_point(data = subset(combo, color == 1), colour = my.cols[2]) + geom_point(data = subset(combo, color == 0), colour = "grey") + geom_text_repel(data = subset(combo, labels == 1), aes(label = as.character(symbol)))  + theme_bw() + ylab("-Log10(p-value)") + xlab("Log2(FC)") + annotate("text", x=1.3, y=0, label= "COVID+", colour = my.cols[2]) + annotate("text", x=-1.2, y=0, label= "COVID-", colour = my.cols[1]) + coord_fixed(ratio = .13) + theme(panel.grid = element_blank())
```

**Figure S1F:**

``` r
plot1
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->

Thus we will use CIBERSORTx estimated cell type fractions as well as an
additional immunoglobulin score to regress non-neutrophil contamination
out of the data.

``` r
source(paste0(prefix,"Pathway_scoring.R"))

IG.score <- Pathway_scoring(IGgenes)
metadata_filtered$IG.score <- (IG.score)
```

We check the correlations with this score with *MZB1*, a plasmablast
marker, and *MS4A1*, a B cell marker.

``` r
goi <- "MZB1"
metadata_filtered$GOI <- t(logTPM_filtered[genepc$Gene.stable.ID[which(genepc$Gene.name == goi)],])
metadata_filtered$GOI <- as.numeric(metadata_filtered$GOI)
summry <- lm(GOI ~ IG.score, metadata_filtered)
p1 <- ggplot(metadata_filtered, aes(x = as.numeric(IG.score), y = GOI)) + geom_point() + theme_bw() + geom_smooth(method = "lm") + xlab("Immunoglobulin Score") + ylab("MZB1 Log2(TPM+1)") + annotate("text", x=-3, y=8.75, label= paste0("R2 =",round(summary(summry)$r.squared,3)), colour = "black") + annotate("text", x=-3, y=8, label= paste0("p =",round(summary(summry)$coefficients[,4],308)), colour = "black")
```

``` r
goi <- "MS4A1"
metadata_filtered$GOI <- t(logTPM_filtered[genepc$Gene.stable.ID[which(genepc$Gene.name == goi)],])
metadata_filtered$GOI <- as.numeric(metadata_filtered$GOI)
summry <- lm(GOI ~ IG.score, metadata_filtered)
p2 <- ggplot(metadata_filtered, aes(x = as.numeric(IG.score), y = GOI)) + geom_point() + theme_bw() + geom_smooth(method = "lm") + xlab("Immunoglobulin Score") + ylab("MS4A1 Log2(TPM+1)") + annotate("text", x=-3, y=5, label= paste0("R2 =",round(summary(summry)$r.squared,3)), colour = "black") + annotate("text", x=-3, y=4.5, label= paste0("p =",round(summary(summry)$coefficients[,4],18)), colour = "black")
```

**Figure S1G:**

``` r
p1
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

``` r
p2
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-53-2.png)<!-- -->

View the immunoglobulin score on a UMAP:

``` r
a = "IG.score"
myPalette <- colorRampPalette((brewer.pal(9, "RdYlBu")))
sc <- scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(min(metadata_filtered[,colnames(metadata_filtered) == a]),max(metadata_filtered[,colnames(metadata_filtered) == a]))) 
p1 <- ggplot(metadata_filtered, aes(x = umap1, y = umap2, colour = (as.numeric(IG.score)))) + theme_bw() + geom_point(size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) + sc
p1$labels$colour <- "IG Score"
```

**Figure S1H:**

``` r
p1
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-55-1.png)<!-- -->

We check the IG score across COVID status, Day, and SeverityMax.

``` r
my.cols <- brewer.pal(3,"Set2")
p1 <- ggplot(metadata_filtered, aes(x = factor(COVID), y = IG.score, fill = factor(COVID))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.15) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .5) + scale_fill_manual(values = c(my.cols[1],my.cols[2])) + stat_compare_means()
```

**Figure S1I:**

``` r
p1
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

``` r
my.cols <- brewer.pal(5,"Set2")
p2 <- ggplot(metadata_filtered, aes(x = factor(Day), y = IG.score, fill = factor(Day))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.15) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .5) + scale_fill_manual(values = c(my.cols[1],my.cols[2],my.cols[3],my.cols[4],my.cols[5])) + stat_compare_means()
```

**Figure Not Included:**

``` r
p2
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

``` r
my.cols <- brewer.pal(3,"RdBu")
p3 <- ggplot(metadata_filtered[metadata_filtered$severity.max %in% c("severe","non-severe"),], aes(x = factor(severity.max), y = IG.score, fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.15) + theme_bw() + theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) + xlab("") + ylab("") + coord_fixed(ratio = .5) + scale_fill_manual(values = c(my.cols[3],my.cols[1])) + stat_compare_means()
```

**Figure Not Included:**

``` r
p3
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->

Now we can perform the COVID+ vs. COVID- differential expression
analysis.

``` r
DESeq2_list <- Neutrophil_DESeq2(counts = Count_filtered, mdata = metadata_filtered, day = "D0")
dds <- DESeqDataSetFromMatrix(countData = DESeq2_list$Count_select, colData = DESeq2_list$coldata, design = ~ Neutrophil_total + T_NK_factor + Monocyte_factor + IG_factor + Plasmablast_factor + COVID)
dds <- DESeq(dds)

res <- as.data.frame(results(dds, name="COVID_Positive_vs_Negative"))
filenam <- "Day0_COVIDpos_vs_COVIDneg_correct-NeuCont+TNK+Monocyte+Plasmablast+IG"
temp <- genepc[which(genepc$Gene.stable.ID %in% rownames(res)),]
res$symbol <- matrix(0L, nrow = nrow(res))
for (i in 1:nrow(res)){
  if (rownames(res)[i] %in% temp$Gene.stable.ID){
    res$symbol[i] <- temp$Gene.name[which(rownames(res)[i] == temp$Gene.stable.ID)] 
  } else {
    res$symbol[i] <- rownames(res)[i]
  }
}
res$rank <- sign(res$log2FoldChange)*(-1)*log10(res$pvalue)
res <- res[complete.cases(res),]
res_sig <- res[res$padj < 0.05,]
#write.table(res,paste0(prefix,"DESeq2/",filenam,".txt"),sep = "\t")

resordered <- res[order(res$rank),]

log2fc <- as.numeric(resordered$log2FoldChange)
log10p <- as.numeric(-1*log10(resordered$pvalue))
pvalue <- as.numeric(resordered$pvalue)
padj <- resordered$padj
rank <- resordered$rank
symbol <- resordered$symbol
combo <- cbind(log2fc,log10p,padj,pvalue,symbol,rank)
colnames(combo) <- c("log2fc","log10p","padj","pvalue","symbol","rank")
rownames(combo) <- rownames(resordered)
combo <- as.data.frame(combo)
combo$rank <- as.numeric(combo$rank)
combo$log10p <- as.numeric(combo$log10p)
combo$log2fc <- as.numeric(combo$log2fc)
combo$pvalue <- as.numeric(combo$pvalue)
combo$significance <- as.numeric(combo$padj < 0.05)
combo$significance <- as.factor(combo$significance)

combo$color <- 0
combo$color[combo$pvalue < 1e-4 & combo$log2fc > 0.5] <- 1
combo$color[combo$pvalue < 1e-4 & combo$log2fc < -0.5] <- -1

combo$labels <- 0
combo$labels[rownames(combo) %in% rownames(resordered)[1:20]] <- 1
combo$labels[rownames(combo) %in% rev(rownames(resordered))[1:20]] <- 1
combo$labels <- as.factor(combo$labels)

options(ggrepel.max.overlaps = Inf)
my.cols <- brewer.pal(3,"Set2")
plot1 <- ggplot(combo, aes(x = log2fc, y = log10p)) + geom_point(data = subset(combo, color == 0), colour = "grey") + geom_point(data = subset(combo, color == 1), colour = my.cols[2]) + geom_point(data = subset(combo, color == -1), colour = my.cols[1]) + geom_text_repel(data = subset(combo, labels == 1), aes(label = as.character(symbol)))  + theme_bw() + ylab("-Log10(p-value)") + xlab("Log2(FC)") + annotate("text", x=1.3, y=0, label= "COVID+", colour = my.cols[2]) + annotate("text", x=-1.2, y=0, label= "COVID-", colour = my.cols[1]) + coord_fixed(ratio = .21) + theme(panel.grid = element_blank())
```

**Figure 1F:**

``` r
plot1
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-63-1.png)<!-- -->

Next we perform Gene Set Enrichment Analysis.

``` r
gmt.file <- gmtPathways(paste0(prefix,"all_gene_sets.gmt"))
res <- read.xlsx(paste0(prefix,"Tables/TableS1.xlsx"), sheet = 16)

ranking <- res[,"rank"]
names(ranking) <- res$symbol
set.seed(15001)
fgseaRes <- fgsea(pathways = gmt.file, 
                  stats = ranking,
                  minSize=25,
                  maxSize=1000,
                  eps = 0)
#write.table(fgseaRes[,1:7], file = paste0(prefix,"GSEA_",filenam,".txt"), sep = "\t")
```

The running enrichment score data points are generated using the
following code.

``` r
getEnrichmentDataframe <- function (pathway, stats, gseaParam = 1, ticksSize = 0.2) {
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    combo <- as.data.frame(cbind(tops,bottoms))
    combo$average <- matrix(0L, nrow = nrow(combo), ncol = 1)
    for (p in 1:nrow(combo)){
      combo$average[p] <- (combo$tops[p]+combo$bottoms[p])/2
    }
    combo <- cbind(combo,pathway)
    return(combo)
}

signalingpathways <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE","GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY","HALLMARK_TNFA_SIGNALING_VIA_NFKB","GO_INTERLEUKIN_1_BETA_PRODUCTION","GO_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","HALLMARK_IL6_JAK_STAT3_SIGNALING","GO_INTERLEUKIN_8_PRODUCTION","GO_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION","GO_INTERLEUKIN_6_PRODUCTION","GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY","GO_RESPONSE_TO_EPIDERMAL_GROWTH_FACTOR")

cellularpathways <- c("GO_DEFENSE_RESPONSE_TO_VIRUS","GO_INNATE_IMMUNE_RESPONSE","GO_REGULATION_OF_VIRAL_GENOME_REPLICATION","GO_OXIDATIVE_PHOSPHORYLATION","GO_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS","GO_REGULATION_OF_CARBOHYDRATE_BIOSYNTHETIC_PROCESS","GO_RIBOSOME_ASSEMBLY","GO_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE")

dataframe <- getEnrichmentDataframe(gmt.file[[signalingpathways[1]]], ranking)
signalingdf <- cbind(dataframe$pathway,dataframe$average,rep(signalingpathways[1],nrow(dataframe)))
colnames(signalingdf) <- c("rank","enrichment","pathway")
for (i in 2:length(signalingpathways)){
  dataframe <- getEnrichmentDataframe(gmt.file[[signalingpathways[i]]],
               ranking)
  temp <- cbind(dataframe$pathway,dataframe$average,rep(signalingpathways[i],nrow(dataframe)))
  colnames(temp) <- c("rank","enrichment","pathway")
  signalingdf <- rbind(signalingdf,temp)
}

dataframe <- getEnrichmentDataframe(gmt.file[[cellularpathways[1]]], ranking)
cellulardf <- cbind(dataframe$pathway,dataframe$average,rep(cellularpathways[1],nrow(dataframe)))
colnames(cellulardf) <- c("rank","enrichment","pathway")
for (i in 2:length(cellularpathways)){
  dataframe <- getEnrichmentDataframe(gmt.file[[cellularpathways[i]]],
               ranking)
  temp <- cbind(dataframe$pathway,dataframe$average,rep(cellularpathways[i],nrow(dataframe)))
  colnames(temp) <- c("rank","enrichment","pathway")
  cellulardf <- rbind(cellulardf,temp)
}
```

``` r
my.cols <- brewer.pal(12,"Paired")
p1 <- ggplot(as.data.frame(signalingdf), aes(x = as.numeric(rank), y = as.numeric(enrichment), colour = pathway)) + geom_point() + theme_bw() + scale_colour_manual(values = rev(my.cols)) + coord_fixed(ratio = 10000) + theme(panel.grid.minor = element_blank(), legend.text=element_text(size=6)) + xlab("Rank in Gene List") + ylab("Running Enrichment Score") + ggtitle("GSEA: Signaling Pathways")
my.cols <- brewer.pal(8,"Dark2")
p2 <- ggplot(as.data.frame(cellulardf), aes(x = as.numeric(rank), y = as.numeric(enrichment), colour = pathway)) + geom_point() + theme_bw() + scale_colour_manual(values = rev(my.cols)) + coord_fixed(ratio = 10000) + theme(panel.grid.minor = element_blank(), legend.text=element_text(size=6)) + xlab("Rank in Gene List") + ylab("Running Enrichment Score") + ggtitle("GSEA: Cellular Processes")
```

**Figure 1G:**

``` r
p1
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

**Figure 1H:**

``` r
p2
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-68-1.png)<!-- -->

Since MHC genes were differentially expressed, we look for MHC metagene
associations with CIBERSORTx estimated cell type fractions. The results
are included in Table S1.

``` r
metadata_filtered$MHC_Class_I <- as.numeric(metadata_filtered$MHC_Class_I)
summry <- lm(MHC_Class_I ~ Neutrophil_total, metadata_filtered)
metadata_filtered$MHC_Class_I <- as.numeric(metadata_filtered$MHC_Class_I)
summry <- lm(MHC_Class_I ~ Mature_Neutrophil, metadata_filtered)
metadata_filtered$MHC_Class_I <- as.numeric(metadata_filtered$MHC_Class_I)
summry <- lm(MHC_Class_I ~ Immature_Neutrophil, metadata_filtered)
metadata_filtered$MHC_Class_I <- as.numeric(metadata_filtered$MHC_Class_I)
summry <- lm(MHC_Class_I ~ T_NK, metadata_filtered)
metadata_filtered$MHC_Class_I <- as.numeric(metadata_filtered$MHC_Class_I)
summry <- lm(MHC_Class_I ~ Monocyte, metadata_filtered)
metadata_filtered$MHC_Class_I <- as.numeric(metadata_filtered$MHC_Class_I)
summry <- lm(MHC_Class_I ~ B, metadata_filtered)
metadata_filtered$MHC_Class_I <- as.numeric(metadata_filtered$MHC_Class_I)
summry <- lm(MHC_Class_I ~ Plasmablast, metadata_filtered)

summry <- lm(MHC_Class_II ~ Neutrophil_total, metadata_filtered)
summry <- lm(MHC_Class_II ~ Mature_Neutrophil, metadata_filtered)
summry <- lm(MHC_Class_II ~ Immature_Neutrophil, metadata_filtered)
summry <- lm(MHC_Class_II ~ T_NK, metadata_filtered)
summry <- lm(MHC_Class_II ~ Monocyte, metadata_filtered)
summry <- lm(MHC_Class_II ~ B, metadata_filtered)
summry <- lm(MHC_Class_II ~ Plasmablast, metadata_filtered)
```

We end by comparing the estimated CIBERSORTx cell type proportions with
severity.

``` r
my.cols <- brewer.pal(3, "RdBu")
p1 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0"),], aes(x = factor(severity.max), y = as.numeric(Neutrophil_total*100), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_continuous(limits = c(0,100))
p2 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D3"),], aes(x = factor(severity.max), y = as.numeric(Neutrophil_total*100), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_continuous(limits = c(0,100))
p3 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D7"),], aes(x = factor(severity.max), y = as.numeric(Neutrophil_total*100), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_continuous(limits = c(0,100))
p4 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0"),], aes(x = factor(severity.max), y = as.numeric(Mature_Neutrophil*100), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_continuous(limits = c(0,100))
p5 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D3"),], aes(x = factor(severity.max), y = as.numeric(Mature_Neutrophil*100), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_continuous(limits = c(0,100))
p6 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D7"),], aes(x = factor(severity.max), y = as.numeric(Mature_Neutrophil*100), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_continuous(limits = c(0,100))
p7 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0"),], aes(x = factor(severity.max), y = as.numeric(Immature_Neutrophil*100), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_continuous(limits = c(0,100))
p8 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D3"),], aes(x = factor(severity.max), y = as.numeric(Immature_Neutrophil*100), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_continuous(limits = c(0,100))
p9 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D7"),], aes(x = factor(severity.max), y = as.numeric(Immature_Neutrophil*100), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_continuous(limits = c(0,100))
```

**Figure 1I:**

``` r
cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=3)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->

We also check the other CIBERSORTx estimated cell type proportions to
see if there is any confounding.

``` r
my.cols <- brewer.pal(3, "RdBu")
p1 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0"),], aes(x = factor(severity.max), y = as.numeric(T_NK*100+1e-5), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2, shape = 16, size = 0.5) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_log10(limits = c(1e-6,100), breaks = c(1e-2,1e-1,1,10,100))
p2 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D3"),], aes(x = factor(severity.max), y = as.numeric(T_NK*100+1e-5), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2, shape = 16) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_log10(limits = c(1e-6,100), breaks = c(1e-2,1e-1,1,10,100))
p3 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D7"),], aes(x = factor(severity.max), y = as.numeric(T_NK*100+1e-5), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2, shape = 16) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_log10(limits = c(1e-6,100), breaks = c(1e-2,1e-1,1,10,100))
p4 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0"),], aes(x = factor(severity.max), y = as.numeric(Monocyte*100+1e-5), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2, shape = 16) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_log10(limits = c(1e-6,100), breaks = c(1e-2,1e-1,1,10,100))
p5 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D3"),], aes(x = factor(severity.max), y = as.numeric(Monocyte*100+1e-5), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2, shape = 16) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_log10(limits = c(1e-6,100), breaks = c(1e-2,1e-1,1,10,100))
p6 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D7"),], aes(x = factor(severity.max), y = as.numeric(Monocyte*100+1e-5), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2, shape = 16) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_log10(limits = c(1e-6,100), breaks = c(1e-2,1e-1,1,10,100))
p7 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0"),], aes(x = factor(severity.max), y = as.numeric(B*100+1e-5), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2, shape = 16) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_log10(limits = c(1e-6,100), breaks = c(1e-2,1e-1,1,10,100))
p8 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D3"),], aes(x = factor(severity.max), y = as.numeric(B*100+1e-5), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2, shape = 16) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_log10(limits = c(1e-6,100), breaks = c(1e-2,1e-1,1,10,100))
p9 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D7"),], aes(x = factor(severity.max), y = as.numeric(B*100+1e-5), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2, shape = 16) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_log10(limits = c(1e-6,100), breaks = c(1e-2,1e-1,1,10,100))
p10 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D0"),], aes(x = factor(severity.max), y = as.numeric(Plasmablast*100+1e-5), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2, shape = 16) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_log10(limits = c(1e-6,100), breaks = c(1e-2,1e-1,1,10,100))
p11 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D3"),], aes(x = factor(severity.max), y = as.numeric(Plasmablast*100+1e-5), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2, shape = 16) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_log10(limits = c(1e-6,100), breaks = c(1e-2,1e-1,1,10,100))
p12 <- ggplot(metadata_filtered[metadata_filtered$COVID == "Positive" & metadata_filtered$Day %in% c("D7"),], aes(x = factor(severity.max), y = as.numeric(Plasmablast*100+1e-5), fill = factor(severity.max))) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.2, shape = 16) + theme_bw() + scale_fill_manual(values = my.cols[c(3,1)]) + theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_blank()) + xlab("") + ylab("") + stat_compare_means() + scale_y_log10(limits = c(1e-6,100), breaks = c(1e-2,1e-1,1,10,100))
```

**Figure Not Included:**

``` r
cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol=3)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-73-1.png)<!-- -->

``` r
cowplot::plot_grid(p7,p8,p9,p10,p11,p12,ncol=3)
```

![](Figure1_S1_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 4.2.0 (2022-04-22)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur/Monterey 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] fgsea_1.22.0                umap_0.2.8.0               
    ##  [3] ggpubr_0.4.0                heatmap3_1.1.9             
    ##  [5] sp_1.5-0                    SeuratObject_4.1.0         
    ##  [7] Seurat_4.1.1                cowplot_1.1.1              
    ##  [9] ggmosaic_0.3.3              corrplot_0.92              
    ## [11] openxlsx_4.2.5              DESeq2_1.36.0              
    ## [13] SummarizedExperiment_1.26.1 Biobase_2.56.0             
    ## [15] MatrixGenerics_1.8.1        matrixStats_0.62.0         
    ## [17] GenomicRanges_1.48.0        GenomeInfoDb_1.32.3        
    ## [19] IRanges_2.30.1              S4Vectors_0.34.0           
    ## [21] BiocGenerics_0.42.0         dplyr_1.0.9                
    ## [23] plyr_1.8.7                  RColorBrewer_1.1-3         
    ## [25] ggrepel_0.9.1               ggplot2_3.3.6              
    ## [27] knitr_1.39                 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.2.2             reticulate_1.25        tidyselect_1.1.2      
    ##   [4] RSQLite_2.2.16         AnnotationDbi_1.58.0   htmlwidgets_1.5.4     
    ##   [7] grid_4.2.0             BiocParallel_1.30.3    Rtsne_0.16            
    ##  [10] munsell_0.5.0          codetools_0.2-18       ica_1.0-3             
    ##  [13] future_1.27.0          miniUI_0.1.1.1         withr_2.5.0           
    ##  [16] spatstat.random_2.2-0  colorspace_2.0-3       progressr_0.10.1      
    ##  [19] highr_0.9              rstudioapi_0.13        ROCR_1.0-11           
    ##  [22] ggsignif_0.6.3         tensor_1.5             listenv_0.8.0         
    ##  [25] labeling_0.4.2         GenomeInfoDbData_1.2.8 polyclip_1.10-0       
    ##  [28] bit64_4.0.5            farver_2.1.1           parallelly_1.32.1     
    ##  [31] vctrs_0.4.1            generics_0.1.3         xfun_0.32             
    ##  [34] fastcluster_1.2.3      productplots_0.1.1     R6_2.5.1              
    ##  [37] locfit_1.5-9.6         bitops_1.0-7           spatstat.utils_2.3-1  
    ##  [40] cachem_1.0.6           DelayedArray_0.22.0    assertthat_0.2.1      
    ##  [43] promises_1.2.0.1       scales_1.2.1           rgeos_0.5-9           
    ##  [46] gtable_0.3.0           globals_0.16.0         goftest_1.2-3         
    ##  [49] rlang_1.0.4            genefilter_1.78.0      splines_4.2.0         
    ##  [52] rstatix_0.7.0          lazyeval_0.2.2         spatstat.geom_2.4-0   
    ##  [55] broom_1.0.0            yaml_2.3.5             reshape2_1.4.4        
    ##  [58] abind_1.4-5            backports_1.4.1        httpuv_1.6.5          
    ##  [61] tools_4.2.0            ellipsis_0.3.2         spatstat.core_2.4-4   
    ##  [64] ggridges_0.5.3         Rcpp_1.0.9             zlibbioc_1.42.0       
    ##  [67] purrr_0.3.4            RCurl_1.98-1.8         rpart_4.1.16          
    ##  [70] openssl_2.0.2          deldir_1.0-6           pbapply_1.5-0         
    ##  [73] zoo_1.8-10             cluster_2.1.3          magrittr_2.0.3        
    ##  [76] data.table_1.14.2      RSpectra_0.16-1        scattermore_0.8       
    ##  [79] lmtest_0.9-40          RANN_2.6.1             fitdistrplus_1.1-8    
    ##  [82] patchwork_1.1.2        mime_0.12              evaluate_0.16         
    ##  [85] xtable_1.8-4           XML_3.99-0.10          gridExtra_2.3         
    ##  [88] compiler_4.2.0         tibble_3.1.8           KernSmooth_2.23-20    
    ##  [91] crayon_1.5.1           htmltools_0.5.3        mgcv_1.8-40           
    ##  [94] later_1.3.0            tidyr_1.2.0            geneplotter_1.74.0    
    ##  [97] DBI_1.1.3              MASS_7.3-58.1          Matrix_1.4-1          
    ## [100] car_3.1-0              cli_3.3.0              parallel_4.2.0        
    ## [103] igraph_1.3.4           pkgconfig_2.0.3        plotly_4.10.0         
    ## [106] spatstat.sparse_2.1-1  annotate_1.74.0        XVector_0.36.0        
    ## [109] stringr_1.4.1          digest_0.6.29          sctransform_0.3.4     
    ## [112] RcppAnnoy_0.0.19       spatstat.data_2.2-0    Biostrings_2.64.1     
    ## [115] rmarkdown_2.15         leiden_0.4.2           fastmatch_1.1-3       
    ## [118] uwot_0.1.11            shiny_1.7.2            lifecycle_1.0.1       
    ## [121] nlme_3.1-159           jsonlite_1.8.0         carData_3.0-5         
    ## [124] viridisLite_0.4.1      askpass_1.1            fansi_1.0.3           
    ## [127] pillar_1.8.1           lattice_0.20-45        KEGGREST_1.36.3       
    ## [130] fastmap_1.1.0          httr_1.4.4             survival_3.4-0        
    ## [133] glue_1.6.2             zip_2.2.0              png_0.1-7             
    ## [136] bit_4.0.4              stringi_1.7.8          blob_1.2.3            
    ## [139] memoise_2.0.1          irlba_2.3.5            future.apply_1.9.0
