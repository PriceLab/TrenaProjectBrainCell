mtx.linc <- get(load("Exc_LINC00507.RData"))
mtx.themis <- get(load("Exc_THEMIS.RData"))
mtx.rorb <- get(load("Exc_RORB.RData"))

dim(mtx.linc); dim(mtx.themis); dim(mtx.rorb)

length(which(is.na(mtx.linc)))  # [1] 19536

mtx.linc[which(is.na(mtx.linc))] <- 0
mtx.themis[which(is.na(mtx.themis))] <- 0
mtx.rorb[which(is.na(mtx.rorb))] <- 0

gene.names <- sort(unique(c(rownames(mtx.linc), rownames(mtx.themis), rownames(mtx.rorb))))
nrow(mtx.themis); nrow(mtx.linc); nrow(mtx.rorb); length(gene.names)  # 14726 14737 14738

mtx.out <- matrix(0, nrow=length(gene.names), ncol=ncol(mtx.themis), dimnames=list(gene.names, colnames(mtx.themis)))

mtx.out[rownames(mtx.themis), colnames(mtx.themis)] <- mtx.out[rownames(mtx.themis), colnames(mtx.themis)] + mtx.themis
mtx.out[rownames(mtx.linc), colnames(mtx.linc)] <- mtx.out[rownames(mtx.linc), colnames(mtx.linc)] + mtx.linc
mtx.out[rownames(mtx.rorb), colnames(mtx.rorb)] <- mtx.out[rownames(mtx.rorb), colnames(mtx.rorb)] + mtx.rorb


# do a spot check - please check more spots yourself!

mtx.linc[1:10, 1:10]; mtx.themis[1:10, 1:10]; mtx.rorb[1:10, 1:10]; mtx.out[1:10, 1:10]

# random checking

 for(i in 1:300){   # will fail only if “row” is the not in the smaller matrix. 
    row <- rownames(mtx.out)[sample(1:nrow(mtx.out),1)]
    col <- colnames(mtx.out)[sample(1:ncol(mtx.out),1)]
    stopifnot(mtx.out[row, col] == mtx.linc[row, col] + mtx.themis[row, col] + mtx.rorb[row,col])
    }

# save the matrix
mtx <- mtx.out
save(file = "Exc_ALL.RData", mtx)


