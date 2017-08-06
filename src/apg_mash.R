### Analyze the mash distnaces



mash <- read.table("data/storage/apg/mash_dist.txt",sep = "\t")
mash <- as.mashdist(mash)
diag(mash) <- 0
mash <- as.dist(mash)
plot(hclust(mash))
