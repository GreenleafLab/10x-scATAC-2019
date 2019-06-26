#Running chromVAR on single-cell summarized experiment
#05/02/19
#Cite Satpathy*, Granja*, et al. 
#Massively parallel single-cell chromatin landscapes of human immune 
#cell development and intratumoral T cell exhaustion (2019)
#Created by Jeffrey Granja
library(chromVAR)
library(SummarizedExperiment)
library(chromVARmotifs)
library(motifmatchr)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
register(SerialParam())
set.seed(1)

#-----------------
# Read Inputs
#-----------------
genome <- BSgenome.Hsapiens.UCSC.hg19
se <- readRDS("results/scATAC-Summarized-Experiment.rds") # single-cell summarized experiment rowRanges as peaks
se <- addGCBias(se, genome = genome)
data("human_pwms_v1")
matches <- matchMotifs(human_pwms_v1, rowRanges(se), genome = "BSgenome.Hsapiens.UCSC.hg19")

#compute deviations
dev <- computeDeviations(object = se, annotations = matches)

#compute variability
metadata(dev)$Variability <- computeVariability(dev)

#add se
metadata(dev)$SummarizedExperiment <- se

#add matches
metadata(dev)$motifMatches <- matches

saveRDS(dev, "results/chromVAR-Summarized-Experiment.rds")
