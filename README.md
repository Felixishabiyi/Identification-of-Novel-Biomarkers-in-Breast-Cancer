# Biomarkers-in-Breast-Cancer

Breast Cancer is referred to as the anomalous growth of cells metamorphizing into tumors; Breast Cancer is a common disease that affects women globally. The evolution of cells from their normal state to the cancerous state can be characterized through biological markers. Biomarkers help to identify biochemical parameters that can be used for diagnosis, prognosis, and therapeutic target identification. In this project, we performed the RNA-seq analysis of the breast cancer data set to identify the expressed genes and also the probable biomarkers that are associated with such genes. This information can further help in diagnostic and therapeutic measures.

# Steps

1.	Data [GSE233242](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233242) was obtained from NCBI
2.	Data preprocessing, cleaning, and manipulation were done using tidyverse and dplyr package.
3.	DeSeq2 Package was implemented to identify the up and down regulated genes
4.	The Diagnostic Plot was derived to have an overview of the normalized DEG data
5.	Functional Enrichment was conducted using G Profiler2 and GO packages.

The Results and plot are in the Images folder of this repository
