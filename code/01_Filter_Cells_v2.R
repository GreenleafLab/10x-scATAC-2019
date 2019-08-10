#Filtering Cells based on TSS enrichment and unique fragments
#05/02/19
#Cite Satpathy*, Granja*, et al. 
#Massively parallel single-cell chromatin landscapes of human immune 
#cell development and intratumoral T cell exhaustion (2019)
#Created by Jeffrey Granja
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(magrittr)
library(ggplot2)
library(Rcpp)
library(viridis)

#--------------------------------------------
# Functions
#--------------------------------------------

sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // [[Rcpp::export]]
  IntegerMatrix tabulate2dCpp(IntegerVector x1, int xmin, int xmax, IntegerVector y1, int ymin, int ymax){
    if(x1.size() != y1.size()){
      stop("width must equal size!");
    }
    IntegerVector x = clone(x1);
    IntegerVector y = clone(y1);
    int n = x.size();
    IntegerVector rx = seq(xmin,xmax);
    IntegerVector ry = seq(ymin,ymax);
    IntegerMatrix mat( ry.size() , rx.size() );
    int xi,yi;
    for(int i = 0; i < n; i++){
      xi = (x[i] - xmin);
      yi = (y[i] - ymin);
      if(yi >= 0 && yi < ry.size()){
        if(xi >= 0 && xi < rx.size()){
          mat( yi , xi ) = mat( yi , xi ) + 1; 
        }
      }
    }
    return mat;
  }'
)

insertionProfileSingles <- function(feature, fragments, by = "RG", getInsertions = TRUE, fix = "center", flank = 2000, norm = 100, smooth = 51, range = 100, batchSize = 100){
  
  insertionProfileSingles_helper <- function(feature, fragments, by = "RG", getInsertions = TRUE, fix = "center", flank = 2000, norm = 100, smooth = 51, range = 100, batchSize = 100){
    #Convert To Insertion Sites
    if(getInsertions){
        insertions <- c(
          GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
          GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
        )
        by <- "RG"
    }else{
      insertions <- fragments
    }
    remove(fragments)
    gc()

    #center the feature
    center <- unique(resize(feature, width = 1, fix = fix, ignore.strand = FALSE))
    
    #get overlaps between the feature and insertions only up to flank bp
    overlap <- DataFrame(findOverlaps(query = center, subject = insertions, maxgap = flank, ignore.strand = TRUE))
    overlap$strand <- strand(center)[overlap[,1]]
    overlap$name <- mcols(insertions)[overlap[,2],by]
    overlap <- transform(overlap, id=match(name, unique(name)))
    ids <- length(unique(overlap$name))
    
    #distance
    overlap$dist <- NA
    minus <- which(overlap$strand == "-")
    other <- which(overlap$strand != "-")
    overlap$dist[minus] <- start(center[overlap[minus,1]]) - start(insertions[overlap[minus,2]])
    overlap$dist[other] <- start(insertions[overlap[other,2]]) - start(center[overlap[other,1]])

    #Insertion Mat
    profile_mat <- tabulate2dCpp(x1 = overlap$id, y1 = overlap$dist, xmin = 1, xmax = ids, ymin = -flank, ymax = flank)
    colnames(profile_mat) <- unique(overlap$name)
    profile <- rowSums(profile_mat)

    #normalize
    profile_mat_norm <- apply(profile_mat, 2, function(x) x/max(mean(x[c(1:norm,(flank*2-norm+1):(flank*2+1))]), 0.5)) #Handles low depth cells
    profile_norm <- profile/mean(profile[c(1:norm,(flank*2-norm+1):(flank*2+1))])

    #smooth
    profile_mat_norm_smooth <- apply(profile_mat_norm, 2, function(x) zoo::rollmean(x, smooth, fill = 1))
    profile_norm_smooth <- zoo::rollmean(profile_norm, smooth, fill = 1)

    #enrichment
    max_finite <- function(x){
      suppressWarnings(max(x[is.finite(x)], na.rm=TRUE))
    }
    e_mat <- apply(profile_mat_norm_smooth, 2, function(x) max_finite(x[(flank-range):(flank+range)]))
    names(e_mat) <- colnames(profile_mat_norm_smooth)
    e <- max_finite(profile_norm_smooth[(flank-range):(flank+range)])

    #Summary
    df_mat <- data.frame(
      enrichment = e_mat,
      insertions = as.vector(table(mcols(insertions)[,by])[names(e_mat)]),
      insertionsWindow = as.vector(table(overlap$name)[names(e_mat)])
      )
    df_sum <- data.frame(bp = (-flank):flank, profile = profile, norm_profile = profile_norm, smooth_norm_profile = profile_norm_smooth, enrichment = e)
    rownames(df_sum) <-  NULL

    return(list(df = df_sum, dfall = df_mat, profileMat = profile_mat_norm, profileMatSmooth = profile_mat_norm_smooth))
  }
  
  uniqueTags <- as.character(unique(mcols(fragments)[,by]))
  splitTags <- split(uniqueTags, ceiling(seq_along(uniqueTags)/batchSize))
  
  pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  batchTSS <- lapply(seq_along(splitTags), function(x){
    setTxtProgressBar(pb, round(x * 100/length(splitTags), 0))
    profilex <- insertionProfileSingles_helper(
        feature=feature, 
        fragments=fragments[which(mcols(fragments)[,by] %in% splitTags[[x]])], 
        by = by, 
        getInsertions = getInsertions,
        fix = fix, 
        flank = flank, 
        norm = norm, 
        smooth = smooth, 
        range = range
      )

    return(profilex)
  })
  df <- lapply(batchTSS, function(x) x$df) %>% Reduce("rbind",.)
  dfall <- lapply(batchTSS, function(x) x$dfall) %>% Reduce("rbind",.)
  profileMat <- lapply(batchTSS, function(x) x$profileMat) %>% Reduce("cbind",.)
  profileMatSmooth <- lapply(batchTSS, function(x) x$profileMatSmooth) %>% Reduce("cbind",.)
  return(list(df = df, dfall = dfall, profileMat = profileMat, profileMatSmooth = profileMatSmooth))
}

#--------------------------------------------
# Input
#--------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
minFrags <- 100
filterFrags <- 1000
filterTSS <- 8
file_fragments <- "data/PBMC_10x-Sub25M-fragments.tsv.gz"
out_fragments <- "data/PBMC_10x-Sub25M-fragments.gr.rds"
name <- "PBMC"

#-----------------
# Reading Fragment Files
#-----------------
message("Reading in fragment files...")
fragments <- data.frame(readr::read_tsv(file_fragments, col_names=FALSE))
#fragmentSub <- fragments[sample(seq_len(nrow(fragments)),25*10^6),]
#write.table(fragmentSub, "data/PBMC_10x-Sub25M-fragments.tsv.gz", col.names=FALSE, row.names =FALSE, sep = "\t", quote = FALSE)

fragments <- GRanges(
  seqnames = fragments[,1], 
  IRanges(fragments[,2]+1, fragments[,3]), 
  RG = fragments[,4], 
  N = fragments[,5]
  )

message("Filtering Lowly Represented Cells...")
tabRG <- table(fragments$RG)
keep <- names(tabRG)[which(tabRG >= minFrags)]
fragments <- fragments[fragments$RG %in% keep,]
fragments <- sort(sortSeqlevels(fragments))

#-----------------
# TSS Profile
#-----------------
feature <- txdb %>% transcripts(.) %>% resize(., width = 1, fix = "start") %>% unique
tssProfile <- insertionProfileSingles(feature = feature, fragments = fragments, 
  getInsertions = TRUE, batchSize = 1000)
tssSingles <- tssProfile$dfall
tssSingles$uniqueFrags <- 0
tssSingles[names(tabRG),"uniqueFrags"] <- tabRG
tssSingles$cellCall <- 0
tssSingles$cellCall[tssSingles$uniqueFrags >= filterFrags & tssSingles$enrichment >= filterTSS] <- 1

#-----------------
# Plot Stats
#-----------------
tssSingles <- tssSingles[complete.cases(tssSingles),]
nPass  <- sum(tssSingles$cellCall==1)
nTotal <- sum(tssSingles$uniqueFrags >= filterFrags)

pdf("results/Filter-Cells.pdf")
ggplot(tssSingles[tssSingles$uniqueFrags > 500,], aes(x = log10(uniqueFrags), y = enrichment)) +
  geom_hex(bins = 100) +
  theme_bw() + scale_fill_viridis() +
  xlab("log10 Unique Fragments") +
  ylab("TSS Enrichment") +
  geom_hline(yintercept = filterTSS, lty = "dashed") +
  geom_vline(xintercept = log10(filterFrags), lty = "dashed") +
  ggtitle(sprintf("Pass Rate : %s of %s (%s)", nPass, nTotal, round(100*nPass/nTotal,2)))
dev.off()

write.table(tssSingles, "results/Filter-Cells.txt")

#Filter
fragments <- fragments[mcols(fragments)$RG %in% rownames(tssSingles)[tssSingles$cellCall==1]]
fragments$RG <- paste0(name,"#",fragments$RG)

#Save
saveRDS(fragments, out_fragments)







