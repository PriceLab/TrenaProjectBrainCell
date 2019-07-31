# convert the counts into an RData file

library(tidyverse)

count.files <- list.files(path = "../expression", pattern = "normalized.tsv", full.names=TRUE)

convertFile <- function(count.file){
        counts <- read.table(file = count.file, header=T, stringsAsFactors=F, sep="\t")
        rownames(counts) <- counts$ensembl_id
        counts$ensembl_id <- NULL
        counts$mean <- rowMeans(counts)
        counts <- counts[(counts$mean > 0.05),]
        counts$mean <- NULL
        mtx <- as.matrix(counts)
        mtx.t <- t(mtx)   
        mtx.t <- scale(mtx.t)
        mtx <- t(mtx.t)
        return(mtx)
}

#converted.counts <- lapply(count.files, convertFile)

# I didn't actually call the function, but just did the following to save the file:
#save(file = "../expression/TCX_counts_scaled.RData", mtx)



countsToScaledLog2p1 <- function(tsvFile) {
  counts <- read.table(tsvFile, header = T, row.names=1, stringsAsFactors = F)
  counts$mean <- rowMeans(counts)
  # the choice of .50 is based on the number of genes it leaves me with. I would expect it to differ based on the dataset.
  cutoff <- quantile(counts$mean, c(.50))
  counts <- counts[(counts$mean > cutoff),]
  counts$mean <- NULL
  mtx <- as.matrix(counts)
  pseudo <- min(1,min(mtx[mtx > 0]))
  mtx.t <- t(mtx)
  mtx.t <- scale(log(pseudo+mtx.t,base=2))
  mtx <- t(mtx.t)
  return(mtx)
}

converted.counts <- lapply(count.files, countsToScaledLog2p1)

#lMp1s <- countsToScaledLog2p1("MayoRNAseq_RNAseq_TCX_geneCounts_normalized.tsv")
#save(lMp1s,file="MayoRNAseq_RNAseq_TCX_geneCounts_normalized-log1p2-scaled.RData")

