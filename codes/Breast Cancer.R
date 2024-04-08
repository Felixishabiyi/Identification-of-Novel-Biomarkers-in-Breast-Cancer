#Importing the count matrix
bres <- read.csv("C:/Users/Administrator/Documents/R/Dataset/GSE233242.csv")

#checking the dimensions of the count matrix
dim(bres)

#checking the column and row names 
colnames(bres)
rownames(bres)

View(bres)

#Loading libraries
library(tidyverse)
library(dplyr)

#subsetting the dataframe from the count matrix
bresss <- bres %>%
  select(-1,-2,-3)

View(bresss)

colSums(bresss)

#loading the geoquery library
library(GEOquery)

bres_meta <- getGEO("GSE233242",GSEMatrix = TRUE)

bres_meta

#Extracting the metadata
meta_bres <- pData(phenoData(bres_meta[[1]]))

View(meta_bres)

#Manipulating/cleaning and preprocessing of the metadata
meta_bres <- meta_bres %>%
  select(1,40,41)%>%
  rename(Title = title) %>%
  rename(Subtype = `subtype:ch1`)%>%
  rename(Tissue = `tissue:ch1`)

#mutate(Tissue = gsub("Human normal breast tissue","Normal",Tissue)) %>%
#mutate(Tissue = gsub("Human tumor breast tissue","Tumor",Tissue))  

View(meta_bres)

meta_bres <- meta_bres %>%
  mutate(Tissue = gsub("Human normal breast tissue","Normal",Tissue)) %>%
  mutate(Tissue = gsub("Human tumor breast tissue","Tumor",Tissue))

View(meta_bres)

library(ggplot2)

#Generating a barplot for the metadata
g <- ggplot(meta_bres, aes(Tissue)) + geom_bar(aes(fill=Subtype), 
                                               position = 'dodge')
g

#creating a matrix for the gene names
bres_row_names <- bres$gene_name

#lableing the count matrix with the genes as rownames
rownames(bresss) <- bres_row_names
View(bresss)

#claculating for the most variable genes
bresss_vari <- apply(bresss, 1, var)

#selecting the top50 variable genes
bresss_top50 <- names(bresss_vari[order(bresss_vari,decreasing = T)][1:50])

library(ggplot2)
library(pheatmap)

#selecting the tissue column in the metadata
meta_bres_annot <- meta_bres %>%
  select(3)

colnames(bresss) <- rownames(meta_bres)

#heatmap for the expression profile for the top50 variable genes in breast CA
pheatmap(bresss[bresss_top50,],scale = 'row',annotation_col = meta_bres[2:3],
         show_colnames = F,cluster_rows = T)

class(meta_bres)

dim(meta_bres)
library(dplyr)

#Extracting meta data for Luminal A
lumina_a <- meta_bres[meta_bres$Subtype %in% c('Luminal A', 'Normal'),]

View(lumina_a)

#rownames of luminal A dataset
rownames(lumina_a)

#Count matrix for Luminal A
lum_a_count <- bresss[rownames(lumina_a)]
View(lum_a_count)
dim(lum_a_count)

#Heatmap Luminal A
Lum_a_var <- apply(lum_a_count, 1, var)

#Top10 genes for the most variable genes in Luminal A 
luma_top10 <- names(Lum_a_var[order(Lum_a_var,decreasing = T)][1:10])

#Heatmap for the Top 10 most variable gene in Luminal A
pheatmap(lum_a_count[luma_top10,],scale = 'row',annotation_col = lumina_a[2:3],
         show_colnames = F,cluster_rows = T)


#Extracting metadata for Luminal B
Lumina_b <- meta_bres[meta_bres$Subtype %in% c('Luminal B', 'Normal'),]
View(Lumina_b)
rownames(Lumina_b)

#Count Matrix for Luminal B
lum_b_count <- bresss[rownames(Lumina_b)]
View(lum_b_count)

#Heatmap Luminal B
Lum_b_var <- apply(lum_b_count, 1, var)

lumb_top10 <- names(Lum_b_var[order(Lum_b_var,decreasing = T)][1:10])

#Top 10 most variable genes in Luminal B
pheatmap(lum_b_count[lumb_top10,],scale = 'row',annotation_col = Lumina_b[2:3],
         show_colnames = F,cluster_rows = T)


#Extracting the metadata for TNBC breast cancer subtype
TNBC <- meta_bres[meta_bres$Subtype %in% c('TNBC', 'Normal'),]
View(TNBC)
rownames(TNBC)

tnbc_count <- bresss[rownames(TNBC)]
View(tnbc_count)

#Heatmap TNBC
tnbc_var <- apply(tnbc_count, 1, var)

#Top10 most variable genes in TNBC
tnbc_top10 <- names(tnbc_var[order(tnbc_var,decreasing = T)][1:10])

#Heatmap for the Top 10 most variable gene in TNBC
pheatmap(tnbc_count[tnbc_top10,],scale = 'row',annotation_col = TNBC[2:3],
         show_colnames = F,cluster_rows = T)

HER2 <- meta_bres[meta_bres$Subtype %in% c('HER2', 'Normal'),]
View(HER2)
rownames(HER2)

her2_count <- bresss[rownames(HER2)]
View(her2_count)

#Heatmap HER2
her2_var <- apply(her2_count, 1, var)

#Heatmap for the Top 10 most variable gene in HER2
her2_top10 <- names(her2_var[order(her2_var,decreasing = T)][1:10])

#Heatmap for the Top 10 most variable gene in HER2
pheatmap(her2_count[her2_top10,],scale = 'row',annotation_col = HER2[2:3],
         show_colnames = F,cluster_rows = T)

View(bresss)

#Converting the count data to matrix
countData <- as.matrix(bresss)

#loading the Deseq2 Library
library(DESeq2)
library(stats)

#Indicating the parameter in the metadata to be the design for Deseq2
designFormula <- "~Tissue"

# meta_bres['Tissue'] <- factor(meta_bres$Tissue)

#forming the Deseq2 object
dds <- DESeqDataSetFromMatrix(countData = round(countData),
                              colData = meta_bres,
                              design = as.formula(designFormula))
print(dds)
row.names(dds)
colnames(dds)

#Remobing rows with sums > 1
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]

#Implementing the Deseq2 function
dds <- DESeq(dds)

#Assinging the DEG to the DE variable
DE = results(dds, contrast = c('Tissue', 'Tumor','Normal'))

#Ordering the results based on the Pvalue
DE <- DE[order(DE$pvalue),]

#Printing the results
print(DE)

#Loading the enhanced volcano package for the volcano plots
library(EnhancedVolcano)

#Implemeneting the package for the volcano plots
EnhancedVolcano(DE, x = 'log2FoldChange', y = 'pvalue', 
                lab = rownames(DE))

#Ploting the Diagnostic plot
DESeq2::plotMA(object = dds, ylim = c(-5, 5))

#Loading the Gprofiler2 package
library(gprofiler2)
library(knitr)

DE <- results(dds, contrast = c('Tissue', 'Tumor', 'Normal'))

#Removong the empty rows with empty values in the Padj column in the DE 
DEr <- DE[!is.na(DE$padj),]

#setting the padj value - cut-off mark
DEr <- DEr[DEr$padj < 0.05,]

#setting the log2foldchange value - cut-off mark
DEr <- DEr[abs(DEr$log2FoldChange) > 2,]

#getting the rownames for the genes within the cutoff point
genesOfInterest <- rownames(DEr)

#Implementing the gprofiler2 package for gene ontology
gostres <- gost(query = genesOfInterest, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE)

#plotting the gene ontology results
gostplot(gostres, capped = TRUE, interactive = TRUE)

View(genesOfInterest)

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

#Extracting gene names 
genes_expressed <- rownames(DE[DE$log2FoldChange > 2,])

#Implementing EnrichedGO for the Biological Processes
Go_result <- enrichGO(gene = genes_expressed, OrgDb = "org.Hs.eg.db",
                      keyType = "SYMBOL", ont = 'BP')

as.data.frame(Go_result)

#Barplot to reveal the BP
plot(barplot(Go_result, showCategory = 20))

#Implementing EnrichedGO for the Cellular Component
Go_result_BP <- enrichGO(gene = genes_expressed, OrgDb = "org.Hs.eg.db",
                      keyType = "SYMBOL", ont = 'CC')
as.data.frame(Go_result_BP)

#Barplot to reveal the CC
plot(barplot(Go_result_BP, showCategory = 20))

#Implementing EnrichedGO for the Molecular Fuction
Go_result_MF <- enrichGO(gene = genes_expressed, OrgDb = "org.Hs.eg.db",
                         keyType = "SYMBOL", ont = 'MF')

as.data.frame(Go_result_MF)

#Barplot to reveal the MF
plot(barplot(Go_result_MF, showCategory = 20))
