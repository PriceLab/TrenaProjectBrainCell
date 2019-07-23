mtx.vip <- get(load("Inh_VIP.RData"))
mtx.pvalb <- get(load("Inh_PVALB.RData"))
mtx.lamp5 <- get(load("Inh_LAMP5.RData"))

dim(mtx.vip); dim(mtx.pvalb); dim(mtx.lamp5)

length(which(is.na(mtx.vip)))  # [1] 19536

mtx.vip[which(is.na(mtx.vip))] <- 0
mtx.pvalb[which(is.na(mtx.pvalb))] <- 0
mtx.lamp5[which(is.na(mtx.lamp5))] <- 0

gene.names <- sort(unique(c(rownames(mtx.vip), rownames(mtx.pvalb), rownames(mtx.lamp5))))
nrow(mtx.pvalb); nrow(mtx.vip); nrow(mtx.lamp5); length(gene.names)  # 14726 14737 14738

mtx.out <- matrix(0, nrow=length(gene.names), ncol=ncol(mtx.pvalb), dimnames=list(gene.names, colnames(mtx.pvalb)))

mtx.out[rownames(mtx.pvalb), colnames(mtx.pvalb)] <- mtx.out[rownames(mtx.pvalb), colnames(mtx.pvalb)] + mtx.pvalb
mtx.out[rownames(mtx.vip), colnames(mtx.vip)] <- mtx.out[rownames(mtx.vip), colnames(mtx.vip)] + mtx.vip
mtx.out[rownames(mtx.lamp5), colnames(mtx.lamp5)] <- mtx.out[rownames(mtx.lamp5), colnames(mtx.lamp5)] + mtx.lamp5


# do a spot check - please check more spots yourself!

mtx.vip[1:10, 1:10]; mtx.pvalb[1:10, 1:10]; mtx.lamp5[1:10, 1:10]; mtx.out[1:10, 1:10]

# random checking

 for(i in 1:300){   # will fail only if “row” is the not in the smaller matrix. 
    row <- rownames(mtx.out)[sample(1:nrow(mtx.out),1)]
    col <- colnames(mtx.out)[sample(1:ncol(mtx.out),1)]
    stopifnot(mtx.out[row, col] == mtx.vip[row, col] + mtx.pvalb[row, col] + mtx.lamp5[row, col])
    }

# save the matrix
mtx <- mtx.out
save(file = "Inh_ALL.RData", mtx)


