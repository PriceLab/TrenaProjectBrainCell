# convert the counts into an RData file

library(tidyverse)
library(biomaRt)
ensemble_ids <- rownames(mtx)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
human_mapping <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), mart = ensembl)

count.files <- list.files(path = "../expression", pattern = "normalized.tsv", full.names=TRUE)

countsToScaledLog2p1 <- function(tsvFile) {
  counts <- read.table(tsvFile, header = T, row.names=1, stringsAsFactors = F)
  counts$mean <- rowMeans(counts)
  
  # the choice of .50 is based on the number of genes it leaves me with. I would expect it to differ based on the dataset.
  cutoff <- quantile(counts$mean, c(.50))
  counts <- counts[(counts$mean > cutoff),]
  
  # convert ensemble ids to geneSymbols and get rid of duplicates that have a lower mean value
  ensemble_ids <- rownames(counts)
  counts$ensembl_gene_id <- rownames(counts)
  rownames(counts) <- NULL
  counts.temp <- left_join(counts, human_mapping)
  counts.temp <- counts.temp %>% dplyr::select(hgnc_symbol, ensembl_gene_id, everything())
  counts.temp$ensembl_gene_id <- NULL
  #count.temp$mean <- rowMeans(counts.temp)
  counts.temp <- counts.temp[order(counts.temp$hgnc_symbol, -abs(counts.temp$mean) ), ]
  counts.temp <- counts.temp[ !duplicated(counts.temp$hgnc_symbol), ]
  counts.temp$mean <- NULL
  present <- which(!(counts.temp$hgnc_symbol == ""))
  counts.temp <- counts.temp[present,]
  rownames(counts.temp) <- counts.temp$hgnc_symbol
  counts.temp$hgnc_symbol <- NULL

  # log transform the data and scale it
  mtx <- as.matrix(counts.temp)
  pseudo <- min(1,min(mtx[mtx > 0]))
  mtx.t <- t(mtx)
  mtx.t <- scale(log(pseudo+mtx.t,base=2))
  mtx <- t(mtx.t)
  return(mtx)
}

converted.counts <- lapply(count.files, countsToScaledLog2p1)

#lMp1s <- countsToScaledLog2p1("MayoRNAseq_RNAseq_TCX_geneCounts_normalized.tsv")
#save(lMp1s,file="MayoRNAseq_RNAseq_TCX_geneCounts_normalized-log1p2-scaled.RData")
