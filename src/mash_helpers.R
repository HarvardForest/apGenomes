as.mashdist <- function(x){
    lab <- unique(as.character(unlist(x[,1:2])))
    mat <- array(NA,dim = rep(length(lab),2))
    for (i in 1:nrow(x)){
        mat[lab == x[i,1],lab == x[i,2]] <- x[i,3]
    }
    rownames(mat) <- colnames(mat) <- lab
    mat
}
