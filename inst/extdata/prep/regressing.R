library(tidyverse)

print(load("cibersort_bootstrap_summary.rda"))
resMean <- as.data.frame(resMean)
resMean$Symbol <- rownames(resMean)

counts <- read.table(file = "../expression/celltype_expr_Micro_TYROBP.tsv", header = T, row.names = 1, stringsAsFactors = F)
counts <- counts[complete.cases(counts), ]
counts$mean <- rowMeans(counts)

# the choice of .50 is based on the number of genes it leaves me with. I would expect it to differ based on the dataset.
cutoff <- quantile(counts$mean, c(.72))
counts <- counts[(counts$mean > cutoff),]
names(counts) <- gsub(x = names(counts), pattern = "^X", replacement = "")
counts$mean <- NULL
counts.mtx <- as.matrix(counts)
counts.mtx <- t(counts.mtx)
counts <- as.data.frame(counts.mtx)
counts$Symbol <- rownames(counts)
count_names <- colnames(counts)

temp1 <- dplyr::inner_join(counts, resMean)
temp1 <- temp1 %>% dplyr::select(Symbol, everything())
ct_names <- colnames(resMean)
ct_sub <- temp1[, ct_names]
ct_sub <- ct_sub %>% dplyr::select(Symbol, everything())
rownames(ct_sub) <- ct_sub$Symbol
ct_sub$Symbol <- NULL
ct_sub <- as.matrix(ct_sub)

counts_sub <- temp1[,count_names]
counts_sub <- counts_sub %>% dplyr::select(Symbol, everything())
rownames(counts_sub) <- counts_sub$Symbol
counts_sub$Symbol <- NULL
counts_sub <- t(counts_sub)
counts_sub <- as.matrix(counts_sub)

# make a function out of this to go through every row in counts_sub
test1 <- lm(counts_sub[1,]~ct_sub)
test1 <- residuals(test1)

