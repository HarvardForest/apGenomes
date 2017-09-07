### Analyze the mash distnaces
source("src/mash_helpers.R")
broad.info <- read.csv("data/storage/apg/broad_sample_key.csv")
broad.info[,"Collaborator.Sample.ID"] <- as.character(broad.info[,"Collaborator.Sample.ID"])
broad.info[broad.info[,"Collaborator.Sample.ID"] == "arudis1","Collaborator.Sample.ID"] <- "rud1"
broad.info[broad.info[,"Collaborator.Sample.ID"] == "rud6","Collaborator.Sample.ID"] <- "rud2"

## mash organziation
geno.info <- read.csv("data/storage/apg/gen_seq_info.csv")
mash.txt <- read.table("data/storage/apg/mash_dist.txt",sep = "\t")
mash <- as.mashdist(mash.txt) 
rownames(mash) <- colnames(mash) <- as.character(geno.info[sapply(geno.info[,1],grep,x = rownames(mash)),2])
ncbi.gen <- mash
ncbi.rv <- c(11,19,14,15,9,5,7,10,6,4,8,12,1,2,24,3,22,23,20,25,26,18,21,16,13,17)
mash <- mash[grep("Aphaenogaster",rownames(mash)),grep("Aphaenogaster",rownames(mash))]



## distances
mash.d <- as.dist(mash)
geo.d <- as.dist(apg.gd)
geo.cd <- as.dist(apg.gcd)
geod.pic <- apg.gcd[rownames(apg.gcd) != "pic1","pic1"]
dist.pic <- apg.geo[,'Latitude'] - apg.geo[,'Latitude']['pic1']
mash.dpic <- mash[rownames(mash) != "pic1","pic1"]

mash <- mash[order(apg.gcd[,"pic1"],decreasing = F),
             order(apg.gcd[,"pic1"], decreasing = F)]
gcd.pic <- apg.gcd["pic1",rownames(apg.gcd) != "pic1"]
apg.gcd <- apg.gcd[order(apg.gcd[,"pic1"]),order(apg.gcd[,"pic1"])]

diag(mash) <- NA
diag(apg.gcd) <- NA

### Reorg geo.ctr for latitude
geo.lat <- geo.ctr[c(6,6,5,4,3,1,2),"Lat"]

### Size similarity
size.d <- as.matrix(dist(as.matrix(gaemr.tab[,"TotalScaffoldLength"])))
size.dpic <- size.d[rownames(mash) == "pic1",rownames(mash) != "pic1"]

### GC similiarty
gc.d <- as.matrix(dist(as.matrix(gaemr.tab[,"AssemblyGC"])))
gc.dpic <- size.d[rownames(mash) == "pic1",rownames(mash) != "pic1"]

### Lat distance

lat.d <- as.matrix(dist(as.matrix(geo.lat)))


