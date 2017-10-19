as.mashdist <- function(x){
}
geno.info <- read.csv("data/storage/apg/gen_seq_info.csv")
mash.txt <- read.table("data/storage/apg/mash_dist.txt",sep = "\t")
mash <- as.mashdist(mash.txt) 
rownames(mash) <- colnames(mash) <- paste0(as.character(geno.info[sapply(geno.info[,1],grep,x = rownames(mash)),2]),c("","","",1,2,rep("",nrow(geno.info) - 5)))
mash <- mash[grep("Aphaenogaster",rownames(mash)),grep("Aphaenogaster",rownames(mash))]
apg.geo  <- do.call(rbind,list('arudis1' = c(-78.9830464,36.0200847),
                 'flo1' = c(-82.031176,29.785325)))
colnames(apg.geo) <- c('Longitude','Latitude')
for (i in 1:nrow(ap.ctr)){
}
geo.cd <- as.dist(apg.gcd)
mash <- mash[order(apg.gcd[,"pic1"],decreasing = F),
             order(apg.gcd[,"pic1"], decreasing = F)]
diag(mash) <- NA
mystack <- ls_prism_data() %>%  prism_stack()  
mycrs <- mystack@crs@projargs
mypoints <- data.frame(id = rownames(apg.geo),
)
coordinates(mypoints) <- c('long', 'lat')
proj4string(mypoints) <- CRS(mycrs)
data <- data.frame(coordinates(mypoints), mypoints$id,
                   extract(mystack, mypoints))
colnames(data)[4:6] <- do.call(rbind,strsplit(colnames(data)[4:6],"_"))[,2]
rownames(data)[rownames(data) == "arudis1"] <- "rud1"
data[,"mypoints.id"] <- as.character(data[,"mypoints.id"])
data[rownames(data) == "rud1","mypoints.id"] <- "rud1"
clim.data <- data
clim.d <- dist(apply(clim.data[,c("ppt","tmax","tmin")],2,function(x) (x - mean(x)) / sd(x)))
temp.d <- dist(apply(clim.data[,c("tmax","tmin")],2,function(x) (x - mean(x)) / sd(x)))
ppt.d <- as.matrix(dist(clim.data[,c("ppt")]))
rownames(ppt.d) <- colnames(ppt.d) <- rownames(as.matrix(clim.d))
ppt.d <- as.dist(ppt.d)
clim.d <- as.matrix(clim.d)
clim.d <- as.dist(clim.d)
mash.d <- as.dist(mash)
mash.d <- as.matrix(mash.d)
mash.d <- mash.d[c(1,2,7,5,3,6,4),c(1,2,7,5,3,6,4)]
mash.d <- as.dist(mash.d)
temp.d <- as.matrix(temp.d)
temp.d <- as.dist(temp.d)
ppt.d <- as.matrix(ppt.d)
ppt.d <- as.dist(ppt.d)
geo.cd <- as.matrix(geo.cd)
geo.cd <- geo.cd[c(1,2,7,5,3,6,4),c(1,2,7,5,3,6,4)]
geo.cd <- as.dist(geo.cd)
mantel.tab <- list(clim.geog = ecodist::mantel(clim.d~geo.cd, nperm = 10000),
                  ppt = ecodist::mantel(mash.d~ppt.d, nperm = 10000))
mantel.tab <- do.call(rbind,mantel.tab)
write.csv(mantel.tab,"results/mantel_tab.csv")
