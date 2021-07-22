#Authors: Benjamin Goldstein, Melyssa Minto, Dr. Martine Tremblay
#Summary: This script performs the differential analysis on the RNA-seq data from the HTSeq gene counts, by processing the data into a normalized counts matrix between samples, as well as creating dataframes of comparisons between conditions (including log2FC of normalized counts and adjusted p-values).
  #The processed data is then visualized in a number of forms, including volcano plots and heatmaps. The end result should be an automation of the differential analysis and initial visualization of the data.
#Note: Processed data to be read in for differential analysis can be found in GEO under accession number GSE180609.

#Sections: Loading packages, reading in data, defining functions, processing data, data visualization

# Loading packages --------------------------------------------------------

#Note: some packages may need to be installed through biocManager::install
library(DESeq2)
library(edgeR)
library(tidyverse)
library(dplyr)
library(ggrepel)
library(ggforce)
library(pheatmap)
library(dendextend)
library(fgsea)
library(RRHO)

# Reading in Data ---------------------------------------------------------

#Dir should be changed to the directory in which the processed data files reside on the user's computer
dir = "gene_reads/"
setwd(dir)

#Construction of an object with sample names, paths to sequencing counts files, and conditions
sampleFiles=grep('.txt',list.files(dir),value=TRUE)
sampleFiles
sampleNames=c("MTGFP2", "MTGFP3", "MTGFP4", "MTMUT1", "MTMUT2", "MTMUT3", "MTMUT4", "MTWT1", "MTWT2", "MTWT3", "MTWT4")
sampleCondition=c("GFP", "GFP", "GFP", "MUT", "MUT", "MUT", "MUT", "WT", "WT", "WT", "WT")
sampleTable=data.frame(sampleName=sampleNames, fileName=sampleFiles, condition=sampleCondition)
head(sampleTable)
factor(sampleTable$condition)

#The gtf file for mRatBN7.2, only used to construct tpm counts.
mRatBN7_2_gtf <- read_table2("ncbi-genomes-2021-04-08/GCF_015227675.2_mRatBN7.2_genomic.gtf", col_names = FALSE, comment = "#")

#Processing the gtf into an object with gene names and exon lengths
gtf = mRatBN7_2_gtf %>% 
  filter(X3 %in% "exon") %>% 
  mutate(exon_length = as.numeric(X5) - as.numeric(X4)) %>% 
  mutate(gene = sub('.*?"([^;"]+)"', "\\1", X10) ) %>% 
  mutate(gene = gsub(";", "", gene)) %>% 
  select(gene, exon_length, X2) %>% 
  group_by(gene) %>%
  summarise_at(vars(exon_length), list(length = sum)) 

# Defining Functions ------------------------------------------------------

#A function to extract normalized tpm from a raw counts file
tpm<-function(DESeqObj){
  # getting the raw counts
  raw_counts = as.data.frame(counts(DESeqObj) )%>% 
    mutate(gene = rownames(.)) %>% 
    left_join(., gtf, by="gene")
  
  # scale to gene(mRNA) length
  rpk = raw_counts[,1:11]/raw_counts$length
  # scale by kilobase per million
  scaleFactor= colSums(rpk, na.rm = TRUE)/1000000
  tpm_counts = rpk/ scaleFactor
  data.frame(genes = raw_counts$gene, tpm_counts)
}

#A function to format a table with differential expression information with tpm counts
formatNormalizedResults <-function(DESeqRes, DESeqObj, p.value=0.05, lfc=1){
  normalized_counts = tpm(DESeqObj)
  data.frame( genes = normalized_counts$genes, 
              DESeqRes, 
              significant = ifelse(DESeqRes$padj<p.value & abs(DESeqRes$log2FoldChange)>1, "yes", "no"),
              normalized_counts[,-1])
}

#A function to create a volcano plot from a DESeq object, colored by significance according to a specified pvalue and lfc cutoff
VolcanoPlot <- function(DESeqObj, p.value, lfc)
{
  as.data.frame(DESeqObj) %>% 
    mutate(label = case_when( padj > p.value & abs(log2FoldChange) < lfc ~ "NS",
                              padj > p.value & abs(log2FoldChange) > lfc  ~ paste0("|Log2FC|>",lfc),
                              padj < p.value & abs(log2FoldChange) < lfc  ~ paste0("FDR < ", p.value),
                              padj < p.value & abs(log2FoldChange) > lfc  ~ paste0("FDR < ", p.value, " & |Log2FC|>",lfc),
                              is.na(padj) ~ "NS")) %>% 
    mutate(point_label= ifelse(label %in% paste0("FDR < ", p.value, " & |Log2FC|>",lfc), genes, ""))%>% 
    ggplot(aes(y=-log(padj), x=log2FoldChange, color=label))+
    geom_point(alpha=0.5, size=2) +
    geom_text_repel(aes(label=point_label), color= "black", size=3)+
    geom_hline(yintercept=-log(p.value), linetype=2, color = "grey50")+
    geom_vline(xintercept=-c(lfc,-lfc), linetype=2, color = "grey50")+
    scale_color_manual(values = c("green3", "royalblue3", "red", "grey30"))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank(),
          legend.position = "top", 
          legend.direction="horizontal", 
          legend.title=element_blank(),
          legend.key=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA))+
    guides(colour = guide_legend(override.aes = list(size=2)))
}

#A function to generate the color scheme of the annotations in the clustered heatmap, according to the number of clusters
generate_color_scheme <- function(nclusters)
{
  colors <- scales::seq_gradient_pal("#E8D3C5", "#583923")(seq(0,1,length.out=nclusters))
  names(colors) <- 1:nclusters
  list(
    Condition = c(MUT = '#EBEBEB', WT = "#7A7A7A", GFP = '#0A0A0A'),
    Change = c(DOWN = "#D6747F", UP = "#83ECAC"),
    Cluster = colors
  )
}

#A function to split the heatmap into a specified number of clusters according to hierarachical clustering, and label genes according to what cluster they are in
label_clusters <- function(DESeqObj, nclusters)
{
  SigObj <- DESeqObj %>%
    filter(Significant != "not sig")
  HeatmapObj <- as.data.frame(logcounts[SigObj$genes, ])
  ClusterHeatmap <- pheatmap(HeatmapObj, silent = TRUE)
  HeatmapObj %>%
    mutate(cluster = cutree(ClusterHeatmap$tree_row, k = nclusters), reg = SigObj$Significant)
}

# DE Analysis & Data Processing -------------------------------------------

#DESeq2 from HTseq Gene Count

ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=dir, design= ~condition)
colData(ddsHTSeq)

#Pre-filtering Data
genes <- row.names(ddsHTSeq)
keep <- rowSums(counts(ddsHTSeq)) >= 20
dds <- ddsHTSeq[keep,]
genes <- genes[keep ]

#Additional filtering (filters out 0s)
ddsHTSeq_dataframe_filtered <- as.data.frame(counts(dds)) %>%
  filter(MTGFP2 != 0 & MTGFP3 != 0 & MTGFP4 != 0 & MTMUT1 != 0 & MTMUT2 != 0 & MTMUT3 != 0 & MTMUT4 != 0 & MTWT1 != 0 & MTWT2 != 0 & MTWT3 != 0 & MTWT4 != 0)
genes <- rownames(ddsHTSeq_dataframe_filtered)
dds <- ddsHTSeq[genes, ]

#Filtering out unnamed/unannotated genes
keep <- !(grepl("Gm[0-9][0-9]",genes) | grepl("Rik",genes) |  grepl("LOC",genes))
dds <- dds[keep,]
genes <- genes[keep ]
cat( "keeping ", sum(keep), " genes")

#Generate raw counts table
countsTable <- dds@assays@data@listData[["counts"]]

#Performing the analysis
dds <- DESeq(dds)
MUTvWT <- results(dds, contrast=c("condition","MUT","WT"))
MUTvGFP <- results(dds, contrast=c("condition","MUT","GFP"))
WTvGFP <- results(dds, contrast=c("condition","WT","GFP"))
head(MUTvWT)
head(MUTvGFP)
head(WTvGFP)

tpmTable <- tpm(dds)
dataNormalizedFormatted <- formatNormalizedResults(MUTvWT, dds, p.value = 0.05, lfc = 1)

#Counts after DESeq2
vsd <- vst(dds, blind =FALSE)
normalized <- assay(vsd)
head(as.data.frame(normalized))

#Summary of Results, Filtering for padj <= 0.05
summary(MUTvWT)
summary(MUTvGFP)
summary(WTvGFP)

#Converting to data frames so the object can be worked with in tidyverse
MUTvWT <- data.frame(MUTvWT@listData)
MUTvGFP <- data.frame(MUTvGFP@listData)
WTvGFP <- data.frame(WTvGFP@listData)

#Outputs number of differentially expressed genes
diff_expressed_MUTvWT = MUTvWT %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)
diff_expressed_MUTvGFP = MUTvGFP %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)
diff_expressed_WTvGFP = WTvGFP %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

#Tidy preparation of data for visualization
MUTvWT = MUTvWT %>% 
  mutate(genes = genes,
         log_mean_counts = log10(baseMean),
         Significant = case_when(
           padj < 0.05 & log2FoldChange > 1 ~ "UP",
           padj < 0.05 & log2FoldChange < -1 ~ "DOWN",
           log2FoldChange < 1 & log2FoldChange > -1 ~ "not sig",
           padj > 0.05 ~ "not sig",
           is.na(padj) ~ "not sig")) %>%
  filter(!is.na(padj))
MUTvGFP = MUTvGFP %>%
  mutate(genes = genes,
         log_mean_counts = log10(baseMean),
         Significant = case_when(
           padj < 0.05 & log2FoldChange > 1 ~ "UP",
           padj < 0.05 & log2FoldChange < -1 ~ "DOWN",
           log2FoldChange < 1 & log2FoldChange > -1 ~ "not sig",
           padj > 0.05 ~ "not sig",
           is.na(padj) ~ "not sig")) %>%
  filter(!is.na(padj))
WTvGFP = WTvGFP %>%
  mutate(genes = genes,
         log_mean_counts = log10(baseMean),
         Significant = case_when(
           padj < 0.05 & log2FoldChange > 1 ~ "UP",
           padj < 0.05 & log2FoldChange < -1 ~ "DOWN",
           log2FoldChange < 1 & log2FoldChange > -1 ~ "not sig",
           padj > 0.05 ~ "not sig",
           is.na(padj) ~ "not sig")) %>%
  filter(!is.na(padj))

# Data Visualization -------------------------------------------------------

#Numbers of genes differentially expressed between the three conditions
table(MUTvWT$Significant)
table(MUTvGFP$Significant)
table(WTvGFP$Significant)

#Volcano plots
MUTvWT_Volcano <- VolcanoPlot(MUTvWT, 0.05, 1) + ggtitle("MUT v WT")
MUTvGFP_Volcano <- VolcanoPlot(MUTvGFP, 0.05, 1) + ggtitle("MUT vs GFP")
WTvGFP_Volcano <- VolcanoPlot(WTvGFP, 0.05, 1) + ggtitle("WT vs GFP")

MUTvWT_Volcano
MUTvGFP_Volcano
WTvGFP_Volcano

#PCA of normalized counts by condition
mat <- normalized
mat.pca<-prcomp(t(mat))
scores <-as.data.frame(mat.pca$x)

pcaPlt = ggplot(scores,aes(x=PC1,y=PC2,color=sampleTable$condition,label=sampleTable$sampleName))+
  geom_point(size = 3) + labs(color='', x = paste0("PCA (", round(summary(mat.pca)$importance[2]*100,2), "%)"), y = paste0("PCA (", round(summary(mat.pca)$importance[5]*100,2), "%)")) 

pcaPlt

#Preparing Heatmaps with Hierarchical Clustering, Annotation of Sample Conditions

#Building logcounts function, selecting most variable genes
logcounts <- cpm(countsTable, log=TRUE)
my_sample_col <- data.frame(Condition = rep(c("GFP", "MUT", "WT"), c(3,4,4)))

#Heatmap of significant genes, and clustering (9 clusters)
color_scheme_9 <- generate_color_scheme(9)
MUTvWT_categories_9 <- label_clusters(MUTvWT, 9)
gene_cluster_ann_9 <- data.frame(Cluster = as.character(MUTvWT_categories_9$cluster), Change = MUTvWT_categories_9$reg)
row.names(gene_cluster_ann_9) <- row.names(MUTvWT_categories_9)
row.names(my_sample_col) <- sampleTable$sampleName
Heatmap_categories_9 <- pheatmap(MUTvWT_categories_9[,-c(12,13)], annotation_col = my_sample_col, annotation_row = gene_cluster_ann_9, annotation_colors = color_scheme_9, cutree_rows = 9, show_colnames = FALSE, show_rownames = FALSE, main = "MUTvWT sig")

#Preparing objects for RRHO (using signed -log10 FDRs for plot)
#I use dplyr:: before select() because I think there's a package conflict with this command
MvW_ready <- MUTvWT %>%
  mutate(sgnlogpadj = -log10(padj) * sign(log2FoldChange)) %>%
  dplyr::select(c(genes, sgnlogpadj))
MvG_ready <- MUTvGFP %>%
  mutate(sgnlogpadj = -log10(padj) * sign(log2FoldChange)) %>%
  dplyr::select(c(genes, sgnlogpadj))
WvG_ready <- WTvGFP %>%
  mutate(sgnlogpadj = -log10(padj) * sign(log2FoldChange)) %>%
  dplyr::select(c(genes, sgnlogpadj))

#Performing the RRHO
#NOTE: RRHO is set to save to home directory with this code.
#Output directory (the argument "outputdir") should be changed according to user preference.
RRHO(MvW_ready, MvG_ready, BY=TRUE, alternative='enrichment', log10.ind = TRUE, plot = TRUE, outputdir = "~", labels = c("Mutant vs Wildtype (rank)", "Mutant vs GFP (rank)"))
##########DONE##########