#Constructing a Trajectory in UMAP projected sub-space
#05/02/19
#Cite Satpathy*, Granja*, et al. 
#Massively parallel single-cell chromatin landscapes of human immune 
#cell development and intratumoral T cell exhaustion (2019)
#Created by Jeffrey Granja
library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(magrittr)
library(edgeR)
set.seed(1)

getQuantiles <- function(x){
    trunc(rank(x))/length(x)
}

alignTrajectory <- function(df, trajectory, filter = 0.05, dof = 250, spar = 1){
    findClosest <- function(x, y, fitx, fity){
        distxy <- sqrt(rowSums(cbind((fitx - x)^2 + (fity - y)^2)))
        idxMin <- which.min(distxy)
        if(idxMin==1){
            idxMin <- idxMin + 1
        }else if(idxMin==length(fitx)){
            idxMin <- idxMin - 1
        }
        if(distxy[idxMin + 1]  < distxy[idxMin - 1]){
            diff <- 1
        }else{
            diff <- -1
        }
        data.frame(idx = idxMin, dist = distxy[idxMin], diff = diff)
    }
    dfAll <- data.frame()
    for(x in seq_along(trajectory)){
        #Subset
        dfx <- df[df$Group==trajectory[x],]
        #Mean Diff Filter
        xmean <- colMeans(dfx[,c(1,2)])
        diffx <- sqrt(colSums((t(dfx[,1:2]) - xmean)^2))
        dfx <- dfx[which(diffx <= quantile(diffx,1 - filter)),]
        #Get diff
        if(x!=length(trajectory)){
            xmean1 <- colMeans(df[df$Group==trajectory[x+1],c(1,2)])
            diffx1 <- sqrt(colSums((t(dfx[,1:2]) - xmean1)^2))
            dfx$time <- (1 - getQuantiles(diffx1)) + x
        }else{
            xmean1 <- colMeans(df[df$Group==trajectory[x-1],c(1,2)])
            diffx1 <- sqrt(colSums((t(dfx[,1:2]) - xmean1)^2))
            dfx$time <- getQuantiles(diffx1) + x
        }
        dfAll <- rbind(dfAll , dfx)
    }
    sx <- smooth.spline(dfAll$time, dfAll$x, df = dof, spar = spar)
    sy <- smooth.spline(dfAll$time, dfAll$y, df = dof, spar = spar)
    dfFit <- data.frame(x = sx[[2]], y = sy[[2]], t = seq_along(sy[[2]]))
    dfTrajectory <- df[df$Group %in% trajectory,]
    dfTime <- lapply(seq_len(nrow(dfTrajectory)), function(x){
        findClosest(dfTrajectory[x,1],dfTrajectory[x,2], dfFit[,1],dfFit[,2])
    }) %>% Reduce("rbind",.)
    dfTime$distQ <- getQuantiles(dfTime$dist)
    dfTrajectory$pseudotime <- 100*getQuantiles(dfTime$idx + matrixStats::rowProds(as.matrix(dfTime[,c("diff","distQ")])))

    out <- list(trajectory=dfTrajectory, fitTrajectory = dfFit)
}

trajectoryStats <- function(mat, trajectory, nSim = 1000, nFeatures = 10000){
    #--------------
    # Functions
    #--------------
    rankTrajectory <- function(mat, trajectory, n = 10000, method = "vec"){
        if(method=="vec"){
            vecRank <- c()
            for(i in seq_along(trajectory)[-length(trajectory)]){
                if(i == 1){
                    rem <- c(trajectory[i])
                }else{
                    rem <- c(trajectory[i], trajectory[seq(1,i-1)])
                }
                trajectoryI <- trajectory[i]
                peaksI <- head(order(mat[,trajectoryI],decreasing=TRUE), n = 10000)
                distI <- rank(sqrt(colSums((mat[peaksI,] - mat[peaksI,trajectoryI])^2))[colnames(mat) %ni% rem])
                vecRank <- c(vecRank, distI[trajectory[i+1]])
            }   
            vecRank
        }else{
            matRank <- matrix(ncol=length(trajectory)-1,nrow=ncol(mat))
            rownames(matRank) <- colnames(mat)
            colnames(matRank) <- paste0("T",seq_along(trajectory[-1]))
            for(i in seq_along(trajectory)[-length(trajectory)]){
                if(i == 1){
                    rem <- c(trajectory[i])
                }else{
                    rem <- c(trajectory[i], trajectory[seq(1,i-1)])
                }
                trajectoryI <- trajectory[i]
                peaksI <- head(order(mat[,trajectoryI],decreasing=TRUE), n = 10000)
                distI <- sqrt(colSums((mat[peaksI,] - mat[peaksI,trajectoryI])^2))[colnames(mat) %ni% rem]
                matRank[names(distI),i] <- distI
            }   
            matRank
        }
    }
    nullTracjectory <- function(trajectory, n = 1000){
        nullTracjectory <- list()
        while(length(nullTracjectory) < n){
            trajx <- sample(trajectory, length(trajectory))
            if(!identical(trajx, trajectory)){
                nullTracjectory[[length(nullTracjectory) + 1]] <- trajx
            }

        }
        nullTracjectory
    }
    #Reverse Trajectory
    trajectory <- rev(trajectory)
    # Compute Null
    message("Computing Null Trajectories...")
    nullT <- nullTracjectory(trajectory, nSim)
    pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
    nullRanks <- Reduce("rbind",lapply(seq_along(nullT), function(x){
        setTxtProgressBar(pb, round(x * 100/length(nullT), 0))
        rankX <- rankTrajectory(mat, nullT[[x]], n = nFeatures)
        data.frame(mean = mean(rankX), median = median(rankX))
    }))
    # Compute Actual
    rankT <- rankTrajectory(mat, trajectory, n = nFeatures)
    rankTMat <- rankTrajectory(mat, trajectory, n = nFeatures, method = "mat")
    pvalMean <- sum(nullRanks$mean < mean(rankT)) / nrow(nullRanks)
    pvalMedian <- sum(nullRanks$median < median(rankT)) / nrow(nullRanks)
    out <- list(
        pvalMean = pvalMean,
        pvalMedian = pvalMedian,
        rankT = rankT,
        rankNull = nullRanks,
        rankTMat = rankTMat
    )
    return(out)
}

groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

'%ni%' <- Negate('%in%')

#----------------------------
# Get Inputs
#----------------------------
se <- readRDS("results/scATAC-Summarized-Experiment.rds")
trajectory <- paste0("Cluster",c(2,3,4,1,6,5))

#Align single cells to Trajectory
df <- data.frame(row.names = colnames(se), x = colData(se)$UMAP1, y = colData(se)$UMAP2, Group = colData(se)$Clusters)
trajAligned <- alignTrajectory(df, trajectory)

df2 <- trajAligned[[1]]
dfT <- trajAligned[[2]]

pdf("results/UMAP-PseudoTime.pdf")
library(ggplot2)
ggplot(df2, aes(x,y,color=pseudotime)) + geom_point() +
    theme_bw() + viridis::scale_color_viridis() +
    geom_path(data=data.frame(dfT), aes(x,y,color=NULL), size= 1, 
        arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")))
dev.off()
saveRDS(trajAligned, "results/Aligned-Trajectory.rds")

#Significance of Trajectory Order
clustMeans <- edgeR::cpm(groupSums(assay(se), colData(se)$Clusters, sparse = TRUE),log=TRUE,prior.count=3)
trajStats <- trajectoryStats(clustMeans, trajectory, nSim = 500, nFeatures = 10000)
rankMat <- trajStats$rankTMat
rankDF <- lapply(seq_len(ncol(rankMat)), function(x){
    data.frame(names=names(rankMat[,x]),ranks=rankMat[,x],t=x)
}) %>% Reduce("rbind",.)
rankV <- lapply(seq_len(ncol(rankMat)), function(x){
    data.frame(ranks=rankMat[rev(trajectory)[x+1],x],t=x)
}) %>% Reduce("rbind",.)

pdf("results/Trajectory_Plot_Distance.pdf", useDingbats = FALSE, width = 12, height = 12)
ggplot(rankDF, aes(t,ranks,color=names)) + 
    geom_point(size = 8) + 
    geom_point(data = rankV, aes(color=NULL), color = "black", size = 2) + 
    geom_line(data = rankV, aes(color=NULL), color = "black", size = 1, 
        lty = "dashed", arrow = arrow(type = "open", length = unit(0.2, "inches"))) + 
    theme_bw() + 
    ylab("Trajectory Distance (Dist logCPM)") + xlab(paste0("Trajectory 5000 Sim Pval = ",trajStats$pvalMean)) + 
    ggtitle(paste0(rev(trajectory),collapse=" -> "))
dev.off()

out <- list(trajStats = trajStats, rankDF = rankDF, rankV = rankV)
saveRDS(out, "results/Trajectory-Stats.rds")



