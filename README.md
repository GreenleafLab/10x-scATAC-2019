# Massively parallel single-cell chromatin landscapes of human immune cell development and intratumoral T cell exhaustion. Nature Biotechnology (Satpathy*, Granja* et al. 2019)

## **Link** : https://www.nature.com/articles/s41587-019-0206-z

## Please cite : Satpathy*, Granja* et al. , Massively parallel single-cell chromatin landscapes of human immune cell development and intratumoral T cell exhaustion. Nature Biotechnology (2019) <br/>

![](Figure1.png)

# Downsampled test data for PBMCs is available (~500 MB)

https://jeffgranja.s3.amazonaws.com/10x-scATAC-share/10x-scATAC-Downsampled-PBMC-hg19-data.zip

# Links To Supplementary Data

## Notes

**.rds** file is an R binarized object to read into R use readRDS(filename)

**SummarizedExperiment** is a class in R see : <br/>https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html

**deviations** (TF chromVAR) is a class in R see : <br/>https://bioconductor.org/packages/release/bioc/html/chromVAR.html

## scATAC-seq Hematopoiesis

**scATAC Summarized Experiment** :<br/>https://changseq.s3.amazonaws.com/Jeff/10x_ScATAC/scATAC_Heme_All_SummarizedExperiment.final.rds

**chromVAR Summarized Experiment** :
<br/>https://changseq.s3.amazonaws.com/Jeff/10x_ScATAC/chromVAR_Heme_All_SummarizedExperiment.final.rds

**Cicero Log2 Gene Acitvity Scores Summarized Experiment** :
<br/>https://changseq.s3.amazonaws.com/Jeff/10x_ScATAC/Log2_Gene_Activity_Heme_All_SummarizedExperiment.final.rds

## scATAC-seq CD34 Hematopoiesis

**scATAC Summarized Experiment** : <br/>https://changseq.s3.amazonaws.com/Jeff/10x_ScATAC/scATAC_CD34_BM_SummarizedExperiment.final.rds

## scATAC-seq BCC Tumor Microenvironment

**scATAC Summarized Experiment** :<br/>https://changseq.s3.amazonaws.com/Jeff/10x_ScATAC/scATAC_TME_All_SummarizedExperiment.final.rds

**chromVAR Summarized Experiment** :
<br/>https://changseq.s3.amazonaws.com/Jeff/10x_ScATAC/chromVAR_TME_All_SummarizedExperiment.final.rds

**scATAC Summarized Experiment** :
<br/>https://changseq.s3.amazonaws.com/Jeff/10x_ScATAC/Log2_Gene_Activity_TME_All_SummarizedExperiment.final.rds

## scATAC-seq BCC Tcells (Exhaustion)

**Cicero Log2 Gene Acitvity Scores Summarized Experiment** :<br/>https://changseq.s3.amazonaws.com/Jeff/10x_ScATAC/scATAC_TME_TCells_SummarizedExperiment.final.rds

# Getting 10x scATAC-seq Bam Files

## 1. Go to NIH GEO Page : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785&holding=F1000&otool=stanford.

![](Step0.png)

## 2. Click on a sample (I am showing SU001_Tumor_Immune_Post).

![](Step1.png)

## 3. Navigate down to the bottom and click on SRA link.

![](Step2.png)

## 4. Navigate down and fine the run SRR.

![](Step3.png)

## 5. Click on "Data Access" Tab and then navigate to "Original Format"

![](Step4.png)




