# convert the counts into an RData file

library(tidyverse)

count.files <- list.files(path = "../expression", pattern = ".tsv", full.names=TRUE)

convertFile <- function(count.file){
	counts <- read.table(file = count.file, header=T, stringsAsFactors=F, sep="\t")
	rownames(counts) <- counts$GeneSymbol
	counts$GeneSymbol <- NULL
	counts$mean <- rowMeans(counts)
	counts <- counts[(counts$mean > 0.05),]
	counts$mean <- NULL
	mtx <- as.matrix(counts)
	mtx <- asinh(mtx)
	return(mtx)	
}

converted.counts <- lapply(count.files, convertFile)

temp <- str_split(count.files, ".tsv") %>% map_chr(`[`, 1)
prefix <-str_split(count.files, "_") %>% map_chr(`[`, 3)
sufix <- str_split(temp, "_") %>% map_chr(`[`, 4)
names(converted.counts) <- paste0(prefix, "_", sufix)

lapply(names(converted.counts), function(x) {
       mtx <- converted.counts[[x]]
       save(mtx, file=paste0('../expression/', x, '.RData'))
  })
