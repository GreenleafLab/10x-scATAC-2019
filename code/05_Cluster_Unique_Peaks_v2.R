#Identifying Cluster Specific ATAC-seq peaks
#05/02/19
#Cite Satpathy*, Granja*, et al. 
#Massively parallel single-cell chromatin landscapes of human immune 
#cell development and intratumoral T cell exhaustion (2019)
#Created by Jeffrey Granja
library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(readr)
library(GenomicRanges)
library(magrittr)
library(data.table)
library(edgeR)
set.seed(1)

groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

groupSds <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE) {
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gs <- lapply(unique(groups), function(x) {
        if (sparse) {
            matrixStats::rowSds(as.matrix(mat[, which(groups == x), drop = F]), na.rm = na.rm)
        }
        else {
            matrixStats::rowSds(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gs) <- unique(groups)
    return(gs)
}

createPseudoBulk <- function(mat, groups, labels, minCells = 100, maxCells = 500, minReps = 3, ceiling = 1, prior.count = 3, nSim = 1000, distMethod = "vars", seed = 1){

  calcDiff <- function(mat, method = "vars"){
    if(tolower(method)=="vars"){
      sum(matrixStats::rowVars(mat))/ncol(mat)
    }else if(tolower(method)=="euclidean"){
      sum(dist(t(mat))/ncol(mat))
    }else{
      stop("Error method not found!")
    }
  }

  sumCells <- function(mat, groups, maxCells = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        idx <- which(groups == x)
        if(!is.null(maxCells)){
          idx <- sample(idx, size = min(maxCells, length(idx)), replace = FALSE)
        }
        if (sparse) {
          Matrix::rowSums(mat[, idx, drop = FALSE], na.rm = na.rm)
        }
        else{
          rowSums(mat[, idx, drop = FALSE], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .) %>% as.matrix
    colnames(gm) <- unique(groups)
    return(gm)
  }

  message(paste0("Setting Seed = ",seed))
  set.seed(seed)

  names(labels) <- colnames(mat)
  groupList <- split(labels,groups)

  if(minReps <= 1){
    stop("Minimum 2 replicates!")
  }

  if(is.numeric(ceiling)){
    message(paste0("Setting ceiling of input matrix to ", ceiling, "!"))
    mat@x[mat@x>ceiling]<-ceiling
  }

  #--------------------------------------------------------------------------
  # Constructing Bulk Pseudo ATAC Matrix v1.0
  #--------------------------------------------------------------------------

  pseudoAll <- lapply(seq_along(groupList), function(x){

    message(sprintf("####################################\n  Groups %s of %s : %s\n####################################",x,length(groupList), names(groupList)[x]))
    start <- Sys.time()
    groupx <- groupList[[x]]
    matx <- mat[,names(groupx)]
    bioReps <- names(table(groupx))[which(table(groupx)>minCells)]
    nToSim <- minReps - length(bioReps)

    if(length(bioReps) >= minReps){
      
      #--------------------------------------------------------------------------
      # If there is enough biological samples passing the minimal cells, great
      # good to go! Just merge into true pseudo bulk replicates!
      #--------------------------------------------------------------------------

      message(sprintf("Found %s Bio Reps which is more or equal than required (%s)", length(bioReps), minReps))
      nBio <- length(bioReps)
      groupBio <- groupx[which(groupx %in% bioReps)]
      pseudoBio <- sumCells(matx[,names(groupBio)],groups=groupBio,sparse=TRUE,na.rm=TRUE, maxCells = maxCells)
      nBio <- table(groupBio)[bioReps]
      if(!is.null(maxCells)){
        nBio[nBio > maxCells] <- maxCells
      }

      colnames(pseudoBio) <- lapply(seq_along(bioReps),function(k){
        paste0(names(groupList)[x],"._.BRep_",bioReps[k],".",nBio[which(names(nBio)==bioReps[k])],".FALSE")
      })
      pseudoMat <- pseudoBio

    }else if(length(bioReps) > 0 & ((length(groupx[groupx %ni% bioReps]) + 1) / min(nToSim, 2)) > minCells){
      
      #--------------------------------------------------------------------------
      # If there is at least 1 biological sample with the minimum cells but not
      # as many as required, we will make pseudo replicates with the true replicate
      # to attempt to capture real biological varation
      #--------------------------------------------------------------------------

      message("PSA : To ensure minimum replicates, simulation must be performed!")
      groupBio <- groupx[which(groupx %in% bioReps)]
      
      nBio <- table(groupBio)
      if(length(bioReps) == 1){
        pseudoBio <- Matrix::rowSums(matx[,names(groupBio)],na.rm=TRUE)
      }else{
        pseudoBio <- sumCells(matx[,names(groupBio)],groups=groupBio,sparse=TRUE,na.rm=TRUE, maxCells = maxCells)
      }
      pseudoBioLog <- edgeR::cpm(pseudoBio, log = TRUE, prior.count = prior.count)

      #Determine how sampling is to be performed, ideally we could have the minimum cells * non bio cells of cells!
      nSplit <- floor(length(groupx[groupx %ni% bioReps])/(nToSim))
      if(!is.null(maxCells)){
        nSplit <- min(nSplit, maxCells)
      }
      if(nSplit < minCells){
        message("Splitting cluster into overlapping cells using sampling with replacement BE CAREFUL OF LOW VARIANCE!!!!")
        replacement <- TRUE
        nSplit <- minCells
      }else{
        replacement <- FALSE
      }

      for(i in seq_len(nSim)){
        
        #Figure out how to split Matrix!
        randOrder <- sample(seq_len(ncol(matx)), ncol(matx))
        if(replacement){
          splitList <- lapply(seq_len(nToSim), function(x) sample(seq_len(ncol(matx)), size = nSplit, replace = replacement))
        }else{
          splitList <- split(randOrder, ceiling(seq_along(randOrder)/nSplit))[seq_len(nToSim)]
        }
        
        if(i == 1){

          pseudoMin <- lapply(seq_len(nToSim), function(j){
            jj <- splitList[[j]]
            Matrix::rowSums(matx[,jj,drop=FALSE],na.rm=TRUE)
          }) %>% Reduce("cbind",.)
          pseudoMinLog <- edgeR::cpm(pseudoMin, log = TRUE, prior.count = prior.count)
          diffMax <- calcDiff(cbind(pseudoBioLog, pseudoMinLog), method = distMethod)

        }else{

          pseudoI <- lapply(seq_len(nToSim), function(j){
            jj <- splitList[[j]]
            Matrix::rowSums(matx[,jj,drop=FALSE],na.rm=TRUE)
          }) %>% Reduce("cbind",.)

          pseudoILog <- edgeR::cpm(pseudoI, log = TRUE, prior.count = prior.count)
          diffI <- calcDiff(cbind(pseudoBioLog, pseudoILog), method = distMethod)
          message(sprintf("Trial %s, Current Distance = %s , Max Distance = %s", i, diffI, diffMax))
          
          if(diffI > diffMax){
            message("Found new maxima pseudo to be conservative...")
            pseudoMin <- pseudoI
            pseudoMinLog <- pseudoILog
            diffMax <- diffI
          }

        }

      }

      pseudoBio <- as.matrix(pseudoBio)
      pseudoMin <- as.matrix(pseudoMin)
      colnames(pseudoBio) <- lapply(seq_along(bioReps),function(k){
        paste0(names(groupList)[x],"._.BRep_",bioReps[k],".",nBio[which(names(nBio)==bioReps[k])],".FALSE")
      })
      colnames(pseudoMin) <- paste0(names(groupList)[x],"._.Rep",seq_len(nToSim),".",nSplit,".",replacement)
      pseudoMat <- cbind(pseudoBio, pseudoMin)

    }else{

      #--------------------------------------------------------------------------
      # If there is not at least 1 sample with enough cells we will bootstrap replicates.
      # This is not preferred as we will have a large underestimate of true biological
      # variation.
      #--------------------------------------------------------------------------

      message("PSA : No representation by at least one separate rep, please be cautious with result!")

      #Determine how sampling is to be performed, ideally we could have the minimum cells * non bio cells of cells!
      nToSim <- minReps
      nSplit <- floor(length(groupx[groupx %ni% bioReps])/(nToSim))
      if(floor(1.5 * minCells) > length(groupx)){
        nToSim <- 2
        nSplit <- floor(2 / 3 * length(groupx))
        message(sprintf("Warning! Group size (%s) is smaller than the 3/2 * minimal number of cells (%s)",length(groupx),floor(1.5 * minCells)))
        message(sprintf("To deal with this, we will sample %s replicates at 2/3 the group size (%s of %s)", nToSim, nSplit, length(groupx)))
      }
      if(!is.null(maxCells)){
        nSplit <- min(nSplit, maxCells)
      }
      if(nSplit < minCells){
        message("Splitting cluster into overlapping cells using sampling with replacement BE CAREFUL OF LOW VARIANCE!!!!")
        replacement <- TRUE
        nSplit <- minCells
      }else{
        replacement <- FALSE
      }

      for(i in seq_len(nSim)){

        #Figure out how to split Matrix!
        randOrder <- sample(seq_len(ncol(matx)), ncol(matx))
        if(replacement){
          splitList <- lapply(seq_len(nToSim), function(x) sample(seq_len(ncol(matx)), size = minCells, replace = replacement))
        }else{
          splitList <- split(randOrder, ceiling(seq_along(randOrder)/nSplit))[seq_len(nToSim)]
        }

        if(i == 1){

          pseudoMin <- lapply(seq_len(nToSim), function(j){
            jj <- splitList[[j]]
            Matrix::rowSums(matx[,jj,drop=FALSE],na.rm=TRUE)
          }) %>% Reduce("cbind",.)
          pseudoMinLog <- edgeR::cpm(pseudoMin, log = TRUE, prior.count = prior.count)
          diffMax <- calcDiff(pseudoMinLog, method = distMethod)

        }else{
          
          pseudoI <- lapply(seq_len(nToSim), function(j){
            jj <- splitList[[j]]
            Matrix::rowSums(matx[,jj,drop=FALSE],na.rm=TRUE)
          }) %>% Reduce("cbind",.)
          pseudoILog <- edgeR::cpm(pseudoI, log = TRUE, prior.count = prior.count)
          diffI <- calcDiff(pseudoILog, method = distMethod)
          
          message(sprintf("Trial %s, Current Distance = %s , Max Distance = %s", i, diffI, diffMax))
          if(diffI > diffMax){
            message("Found new maxima pseudo to be conservative...")
            pseudoMin <- pseudoI
            pseudoMinLog <- pseudoILog
            diffMax <- diffI
          }

        }
      }

      pseudoMin <- as.matrix(pseudoMin)
      colnames(pseudoMin) <- paste0(names(groupList)[x],"._.Rep_",seq_len(nToSim),".",nSplit,".",replacement)
      pseudoMat <- pseudoMin

    }

    print(Sys.time() - start)

    pseudoMat 

  }) %>% Reduce("cbind",.)

  out <- list(pseudoMat = pseudoAll, groupMat = sumCells(mat, groups=groups, sparse=TRUE, na.rm=TRUE))

  return(out)

}


uniqueFeatures <- function(
  mat, groups, padj = 0.01, minSdRatio = 0.001, 
  minLFC = 0.25, zCutoff = 1, breakPt = "last",
  padjMethod = "fdr", clusterCols = TRUE, 
  sparse = FALSE, twoWay = FALSE, groupMin = 10, 
  minGroupSize = 1, maxGroupSize = NULL){
  
  #----------------------
  # Functions
  #----------------------
  binarizeMatrix <- function(matMean, matSd, cutoff, method, minSdRatio, minLFC){
    
    binarizeVector <- function(vMean, vSd, cutoff, method, minSdRatio, minLFC){
      
      #Method
      if(method == "mean"){
        vTest = vMean
      }else{
        vTest <- vMean - cutoff * vSd    
      }

      #Order lowest to highest
      idx <- order(vTest)
      vTest <- vTest[idx]
      vMean <- vMean[idx]
      vSd <- vSd[idx]

      #Set which are too low for evaluation ie low Sd that is probably due to 0s or weird offsets
      sdToLow <- which(vSd < minSdRatio*vMean | vSd < 10^-5)
      vBinarySd <- rep(1,length(vSd))
      if(length(sdToLow) > 0){
        vBinarySd[sdToLow] <- 0
      }

      #Initialize
      vMeanCutSd <- vMean + cutoff*vSd;
      maxValue <- vMeanCutSd[1]
      maxBSd <- vBinarySd[1]
      maxMean <- vMean[1]

      #Create out vector and initialize breakpoint
      n <- length(vMean)
      out <- rep(0, n)
      breakPoint <- 0
      out[1] <- breakPoint

      #Evaluate
      for(i in seq(2,n)){
        
        #Check if break point assuming log space
        if((vTest[i] - maxValue) > 0 & (vMean[i] - maxMean) >= minLFC & maxBSd != 0){
          breakPoint <- breakPoint + 1
        }

        #Set current value of break point
        out[i] <- breakPoint

        #Keep Max value observed
        if(vMeanCutSd[i] > maxValue){
          maxValue <- vMeanCutSd[i]
          maxMean <- vMean[i]
          maxBSd <- vBinarySd[i]
        }
      }

      out <- out[order(idx)]
      return(out)

    }

    #Create binary matrix
    bMat <- matrix(NA,nrow=nrow(matMean),ncol=ncol(matMean))
    for(i in seq_len(nrow(matMean))){
      if(i%%5000==0){message(sprintf("%s of %s (percent = %s)", i, nrow(bMat), round(100*i/nrow(bMat),2)))}
      bMat[i,] <- binarizeVector(matMean[i,],matSd[i,],cutoff, method, minSdRatio, minLFC) 
    }

    #Add names
    colnames(bMat) <- colnames(matMean)
    rownames(bMat) <- rownames(matMean)
    return(bMat)

  }

  idxRow <- seq_len(nrow(mat))
  #-----------------------------------------------
  #Within Group Statistics
  #-----------------------------------------------
  message("Getting Within Group Stats...")
  intraMean <- groupMeans(mat, groups = groups, sparse = sparse, na.rm = TRUE)
  rownames(intraMean) <- rownames(mat)

  intraSd <- groupSds(mat, groups = groups, sparse = sparse, na.rm = TRUE)
  rownames(intraSd) <- rownames(mat)
  
  #-----------------------------------------------
  #Binarize Rows of Matrix
  #-----------------------------------------------
  message("Binarizing Features...")
  if (twoWay) {
    binarizedMat <- binarizeMatrix(intraMean, intraSd, zCutoff, "meanSd", minSdRatio, minLFC)
  }else {
    binarizedMat <- binarizeMatrix(intraMean, intraSd, zCutoff, "mean", minSdRatio, minLFC)
  }
  colnames(binarizedMat) <- colnames(intraMean)
  rownames(binarizedMat) <- rownames(intraMean)
  for(i in seq_len(nrow(binarizedMat))){
    bvi <- binarizedMat[i,]
    if(tolower(breakPt) == "last"){
      bvi[bvi!=max(bvi)] <- 0
      bvi[bvi>0] <- 1
    }else{
      bvi[bvi<1] <- 0
      bvi[bvi>0] <- 1
    }
    binarizedMat[i,] <- bvi
  }
  message(sprintf("Successful Binarization of %s Features...", sum(rowSums(binarizedMat) > 0)))

  #-----------------------------------------------
  #Get Test Statistics
  #-----------------------------------------------
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  groupSplit <- split(colnames(mat), groups)
  bpval <- unlist(lapply(seq_len(nrow(binarizedMat)), function(i){
    setTxtProgressBar(pb,round(i*100/nrow(binarizedMat),0))
    if(any(binarizedMat[i,]>0)){
      cu <- as.character(unlist(groupSplit[names(which(binarizedMat[i,]==max(binarizedMat[i,])))]))
      cd <- as.character(unlist(groupSplit[names(which(binarizedMat[i,]!=max(binarizedMat[i,])))]))
      mati <- mat[i,,drop=F]
      pval <- t.test(mati[,cu],mati[,cd])$p.value
    }else{
      pval <- 1
    }
    pval
  }))
  bpadj <- p.adjust(bpval, method=padjMethod) #sum(binarizedMat + noSdMat < 1))

  message(sprintf("\nFiltering by signficance %s of %s...", sum(bpadj < padj), sum(bpval < padj)))
  idxKeep <- which(bpadj < padj)
  bpadj <- bpadj[idxKeep]
  idxRow <- idxRow[idxKeep]
  binarizedMat <- binarizedMat[idxKeep, ]
  mat <- mat[idxKeep, ]
  intraMean <- intraMean[idxKeep, ]
  intraSd <- intraSd[idxKeep, ]

  #-----------------------------------------------
  #Determine which rows are above min group size and max group size to keep
  #-----------------------------------------------
  if (!is.null(maxGroupSize)) {
    idxKeep <- which(rowSums(binarizedMat) < (maxGroupSize + 1) & rowSums(binarizedMat) > (minGroupSize - 1))
  }else {
    idxKeep <- which(rowSums(binarizedMat) > (minGroupSize - 1))
  }
  message(sprintf("Filtering Features that are within Group Size %s of %s...", length(idxKeep), nrow(binarizedMat)))
  idxRow <- idxRow[idxKeep]
  bpadj <- bpadj[idxKeep]
  binarizedMat <- binarizedMat[idxKeep, ]
  mat <- mat[idxKeep, ]
  intraMean <- intraMean[idxKeep, ]
  intraSd <- intraSd[idxKeep, ]
  
  #-----------------------------------------------
  #Determine Pattern Occurences using data.table
  #-----------------------------------------------
  binarizedTable <- data.table::as.data.table(binarizedMat)[,.N,by=c(colnames(binarizedMat))]
  binarizedTable <- data.frame(binarizedTable[which(binarizedTable$N > groupMin),])[,which(colnames(binarizedTable) %ni% "N")]
  idxKeep <- unlist(lapply(seq_len(nrow(binarizedTable)), function(x){
    idx1  <- which(binarizedTable[x,,drop=TRUE] > 0)
    rs1   <- which(rowSums(binarizedMat[,idx1,drop=FALSE]) == length(idx1))
    rs2   <- which(rowSums(binarizedMat[,-idx1,drop=FALSE]) == 0)
    idxBM <- intersect(rs1,rs2)
    idxBM
  }))
  message(sprintf("Filtering Features Pattern Appearances %s of %s...", length(idxKeep), nrow(binarizedMat)))
  idxRow <- idxRow[idxKeep]
  bpadj <- bpadj[idxKeep]
  binarizedMat <- binarizedMat[idxKeep, ]
  mat <- mat[idxKeep, ]
  intraMean <- intraMean[idxKeep, ]
  intraSd <- intraSd[idxKeep, ]
  
  #-----------------------------------------------
  #Organize the output for maximum interpretation
  #-----------------------------------------------
  message(sprintf("Found %s unique elements!", length(idxKeep)))
  message("Finalizing Output...")
  colClust <- hclust(as.dist(1 - cor(intraMean)))
  colOrder <- unique(groups)[colClust$order]
  binarizedMat <- binarizedMat[, colClust$order]
  idxOrdered <- do.call("order", c(as.data.frame(binarizedMat)[seq_len(ncol(binarizedMat))], list(decreasing = TRUE)))
  binarizedMat <- binarizedMat[idxOrdered, ]
  idxRow <- idxRow[idxOrdered]
  intraMean <- intraMean[idxOrdered, colClust$order]
  mat <- mat[idxOrdered, order(match(groups, colOrder))]
  bpadj <- bpadj[idxOrdered]
  
  #-----------------------------------------------
  #Time to label each row
  #-----------------------------------------------
  binarizedTable <- data.frame(data.table::as.data.table(binarizedMat)[,.N,by=c(colnames(binarizedMat))])
  rownames(binarizedTable) <- paste0("Feature_",seq_len(nrow(binarizedTable)))
  dfUniuqe <- lapply(seq_len(nrow(binarizedTable)), function(x){
   idx1 <- which(binarizedTable[x,-ncol(binarizedTable),drop=TRUE] > 0)
   rs1 <- which(rowSums(binarizedMat[,idx1,drop=FALSE]) == length(idx1))
   rs2 <- which(rowSums(binarizedMat[,-idx1,drop=FALSE])==0)
   idxBM <- intersect(rs1,rs2)
   data.frame(feature = rownames(binarizedTable)[x], rows = idxBM)
  }) %>% Reduce("rbind",.)

  #-----------------------------------------------
  #Return Output
  #-----------------------------------------------
  out <- list(mat = mat, binaryMat = binarizedMat, groupMat = intraMean, binarizedTable = binarizedTable, dfFeature = dfUniuqe, rowOrder = idxRow, padj = bpadj)
  return(out)

}  

'%ni%' <- Negate('%in%')

#----------------------------
# Get Inputs
#----------------------------
se <- readRDS("results/scATAC-Summarized-Experiment.rds")

#----------------------------
# PseudoBulk
#----------------------------
colData(se)$Group <- "PBMC" #Replicate info
objPB <- createPseudoBulk(
    mat = assay(se), 
    groups=paste0(colData(se)$Clusters), 
    labels=colData(se)$Group, 
    ceiling = 1, 
    minCells = 100, 
    maxCells = 500, 
    minReps = 2, 
    prior.count = 3, 
    nSim = 250
  )

cdPB <- DataFrame(
    row.names=colnames(objPB[[1]]), 
    Group = stringr::str_split(colnames(objPB[[1]]), "._.", simplify = TRUE)[,1],
    Type  = stringr::str_split(stringr::str_split(colnames(objPB[[1]]), "._\\.", simplify = TRUE)[,2], "\\.", simplify = TRUE)[,1],
    Cells = as.integer(stringr::str_split(stringr::str_split(colnames(objPB[[1]]), "._\\.", simplify = TRUE)[,2], "\\.", simplify = TRUE)[,2]),
    Replace = stringr::str_split(stringr::str_split(colnames(objPB[[1]]), "._\\.", simplify = TRUE)[,2], "\\.", simplify = TRUE)[,3]
  )

sePB <- SummarizedExperiment(assays = SimpleList(counts = objPB[[1]]), rowRanges = rowRanges(se), colData = cdPB)
metadata(sePB)$matCS <- objPB[[2]] 

saveRDS(sePB, "results/Cluster_PseudoBulk-Summarized-Experiment.rds")

#----------------------------
# Unique Features
#----------------------------
params <- list(
    padj = 0.01,
    minSdRatio = 0.001,
    minLFC = 0.25,
    zCutoff = 1,
    breakPt = "last",
    groupMin = 25,
    maxGroupSize = max(floor(length(unique(colData(se)$Group))/3), 2),
    date = Sys.Date()
  )

#Unique Features
uf <- uniqueFeatures(
    edgeR::cpm(assay(sePB),log=TRUE,prior.count=3),
    groups = colData(sePB)$Group,
    padj = params$padj,
    minSdRatio = params$minSdRatio,
    minLFC = params$minLFC,
    zCutoff = params$zCutoff,
    breakPt = params$breakPt,
    groupMin = params$groupMin,
    maxGroupSize = params$maxGroupSize
  )

library(pheatmap)
pdf("results/Unique_Peaks.pdf",width=4,height=8)
m <- uf$mat
z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
z[z > 2] <- 2
z[z < -2] <- -2
pheatmap(z, cluster_cols=FALSE, cluster_rows=FALSE, show_rownames = FALSE)
dev.off()

saveRDS(uf, "results/Unique_Peaks.rds")
saveRDS(params, "results/Unique_Peaks_Params.rds")




