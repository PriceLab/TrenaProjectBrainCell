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
save(file = "../expression/TCX_counts_scaled.RData", mtx)

