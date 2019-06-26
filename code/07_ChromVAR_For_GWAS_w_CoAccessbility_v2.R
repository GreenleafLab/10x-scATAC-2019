#Running gwas chromVAR on single-cell summarized experiment using co-accessibility
#05/02/19
#Cite Satpathy*, Granja*, et al. 
#Massively parallel single-cell chromatin landscapes of human immune 
#cell development and intratumoral T cell exhaustion (2019)
#Created by Jeffrey Granja
library(chromVAR)
library(Matrix)
library(MatrixStats)
library(SummarizedExperiment)
library(magrittr)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
register(SerialParam())
set.seed(1)

se <- readRDS("results/scATAC-Summarized-Experiment.rds")
gr <- readRDS("data/PICS_GWAS_SNPS.gr.rds")
conn <- readRDS("results/Peaks-Co-Accessibility.rds")
conn <- conn[conn[,3] >= 0.35,]
peaknames <- paste(seqnames(se),start(se),end(se),sep="_")
conn[,4] <- match(paste0(conn[,1]), peaknames)
conn[,5] <- match(paste0(conn[,2]), peaknames)
connMat <- Matrix::sparseMatrix(i=conn[,4],j=conn[,5],x=rep(1,nrow(conn)),dims=c(nrow(se),nrow(se)))

#Add Bias
genome <- BSgenome.Hsapiens.UCSC.hg19
se <- addGCBias(se, genome = genome)

#Overlap GWAS SNPs
o <- lapply(split(gr, gr$disease), function(x){
	#extend snp +- 500 bp
	overlapsAny(se, resize(x,1001,"center"), ignore.strand = TRUE)
}) %>% Reduce("cbind",.)
ow <- which(o > 0, arr.ind=TRUE)
matches <- Matrix::sparseMatrix(i=ow[,1],j=ow[,2],x=o[cbind(ow[,1],ow[,2])], dims = c(nrow(se),length(unique(gr$disease))))
colnames(matches) <- names(split(gr, gr$disease))

#Use connections mat!
matches2 <- matches
idx <- which(Matrix::rowSums(matches2)>0) #which peaks have a snp
for(i in seq_along(idx)){
	if(i %% 100 == 0){message(sprintf("%s of %s",i,length(idx)))}
	#peaks co-accessible to peak with snp
	coi <- unique(c(which(connMat[,idx[i]]>0),which(connMat[idx[i],]>0)))
	if(length(coi) > 0){
		#create sub mat
		mati <- as(t(replicate(length(coi), matches[idx[i],])),"dgCMatrix")
		#add it to sub mat of connected peaks
		matches2[coi,,drop=FALSE] <- matches2[coi,,drop=FALSE] + mati
	}
}
diff <- Matrix::colSums(matches2) - Matrix::colSums(matches)
print(diff)

#Make Annotation SE
anno_ix <- SummarizedExperiment(assays=SimpleList(motifMatches=matches2), rowRanges=rowRanges(se))

#Compute Deviations
dev <- computeDeviations(se, anno_ix)

#compute variability
metadata(dev)$Variability <- computeVariability(dev)

#add matches
metadata(dev)$gwasMatches <- anno_ix

#save output
saveRDS(dev, "results/GWAS-Co-Accessibility-chromVAR-Summarized-Experiment.rds")



