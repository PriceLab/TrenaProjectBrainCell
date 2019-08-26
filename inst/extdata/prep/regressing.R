library(tidyverse)
library(Biobase)
library(ggplot2)

print(load("cibersort_bootstrap_summary.rda"))
resMean <- as.data.frame(resMean)
resMean$Symbol <- rownames(resMean)

counts <- read.table(file = "../expression/celltype_expr_Micro_TYROBP.tsv", header = T, row.names = 1, stringsAsFactors = F)
counts <- counts[complete.cases(counts), ]
counts <- as.matrix(counts)
counts.m <- rowMedians(counts)
counts <- cbind(counts.m, counts)
counts <- as.data.frame(counts)

# trying to figure out the best way to get rid of lowely expressed genes. Cutoff differs a lot between datasets
counts <- subset(counts, counts.m > 1)

#cutoff <- quantile(counts$counts.m, c(.72))
#counts <- counts[(counts$median > cutoff),]

names(counts) <- gsub(x = names(counts), pattern = "^X", replacement = "")
counts$counts.m <- NULL
counts.mtx <- as.matrix(counts)
counts.mtx <- t(counts.mtx)
counts <- as.data.frame(counts.mtx)
counts$Symbol <- rownames(counts)
count_names <- colnames(counts)

# get only the samples for which there are cell percentages
temp1 <- dplyr::inner_join(counts, resMean)
temp1 <- temp1 %>% dplyr::select(Symbol, everything())
ct_names <- colnames(resMean)
ct_sub <- temp1[, ct_names]
ct_sub <- ct_sub %>% dplyr::select(Symbol, everything())
rownames(ct_sub) <- ct_sub$Symbol
ct_sub$Symbol <- NULL
ct_sub <- as.matrix(ct_sub)

# get only the counts for which there are cell percentages
counts_sub <- temp1[,count_names]
counts_sub <- counts_sub %>% dplyr::select(Symbol, everything())
rownames(counts_sub) <- counts_sub$Symbol
counts_sub$Symbol <- NULL
counts_sub <- t(counts_sub)
counts_sub <- as.matrix(counts_sub)

# make a function out of this to go through every row in counts_sub
#test1 <- lm(counts_sub[1,]~ct_sub)
#test2 <- residuals(test1)

# because the lm functions works by row but the output is by column, it must be transposed
adj.counts <- t(apply(counts_sub, 1, function(x) lm(x~ct_sub) %>% residuals(.)))
