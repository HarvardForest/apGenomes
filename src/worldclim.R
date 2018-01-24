library(raster)
library(sp)
library(ggmap)

### Functions
get_names <- function(x){
    x <- strsplit(x, split = " ")[[1]]
    paste(substr(x[1],1,1), x[2], sep = "_")
}

### Get locations for gaga_genome_info.csv
ncbi_info <- read.csv("../data/gaga_genome_info.csv")
locs <- ncbi_info[,"Location"]
locs <- sapply(as.character(locs),function(x) 
    paste(strsplit(x,",")[[1]][c(1,3)], collapse = ","))
ncbi_gps <- list()
for (i in 1:length(locs)){
    ncbi_gps[[i]] <- geocode(locs[i])
}
for (i in 1:length(locs)){
    if (any(is.na(ncbi_gps[[i]]))){
        ncbi_gps[[i]] <- geocode(locs[i])
    }
}
ncbi_gps <- do.call(rbind,ncbi_gps)
ncbi_gps[grepl("Unkown",rownames(ncbi_gps)),] <- c(NA,NA)
out <- data.frame(ncbi_info,ncbi_gps)
write.csv(out,"../data/gaga_genome_info_test.csv", row.names = FALSE)

ncbi.spc <- SpatialPoints(ncbi_gps, proj4string = r@crs)
ncbi.df <- cbind.data.frame(coordinates(ncbi_gps),ncbi.spc)

pdf(file = "tmp.pdf")
plot(r[[1]])
plot(ncbi.spc, add = TRUE, pch = "")
text(ncbi.spc, cex = 0.5, labels = substr(
                   gsub(",, ", " ", 
                        gsub("Unknown", "", rownames(ncbi_gps))),1,15))
dev.off()
system("scp tmp.pdf matthewklau@fas.harvard.edu:public_html")


ncbi.gps <- read.csv("data/gaga_genome_info.csv")
ncbi.gps[,"Species.name"] <- as.character(ncbi.gps[,"Species.name"])

### Merge ncbi and apg geo
apg.merging <- data.frame(Species.name = rownames(as.matrix(mash.d)), 
                          lon = apg.geo[,'Longitude'],
                          lat = apg.geo[,'Latitude'])
all.geo <- rbind(ncbi.gps[,c("Species.name","lon","lat")], 
                 apg.merging)
all.geo <- all.geo[match(rownames(ncbi.gen),all.geo[,"Species.name"]),]
all.geo <- na.omit(all.geo)
all.mash <- ncbi.gen
all.mash <- all.mash[sapply(rownames(all.mash), function(x) x %in% all.geo[,"Species.name"]),
                     sapply(colnames(all.mash), function(x) x %in% all.geo[,"Species.name"])]
all.mash <- all.mash[na.omit(match(all.geo[,"Species.name"],rownames(all.mash))),
         na.omit(match(all.geo[,"Species.name"],rownames(all.mash)))]
all.geo <- all.geo[sapply(all.geo[,"Species.name"], function(x) x %in% rownames(all.mash)),]

all(all(rownames(all.mash) == colnames(all.mash)) & 
        all(rownames(all.mash) == all.geo[,"Species.name"]))

### Get climate for ncbi locs
r <- getData("worldclim",var="bio",res=2.5)
wc <- r 
names(wc) <- c("MAT", "MDR", "Iso", "TS", "Tmax", "Tmin", "ATR", "MTWeQ", "MTDQ","MTWaQ", "MTCQ", "PA", "PWM", "PDM", "PS", "PWeQ", "PDQ", "PWaQ", "PCQ")
## BIO1 = Annual Mean Temperature "MAT"
## BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp)) "MDR"
## BIO3 = Isothermality (BIO2/BIO7) (* 100) "Iso"
## BIO4 = Temperature Seasonality (standard deviation *100) "TS"
## BIO5 = Max Temperature of Warmest Month "Tmax"
## BIO6 = Min Temperature of Coldest Month "Tmin"
## BIO7 = Temperature Annual Range (BIO5-BIO6) "ATR"
## BIO8 = Mean Temperature of Wettest Quarter "MTWeQ"
## BIO9 = Mean Temperature of Driest Quarter "MTDQ
## BIO10 = Mean Temperature of Warmest Quarter "MTWaQ"
## BIO11 = Mean Temperature of Coldest Quarter "MTCQ"
## BIO12 = Annual Precipitation "PA"
## BIO13 = Precipitation of Wettest Month "PWM"
## BIO14 = Precipitation of Driest Month "PDM"
## BIO15 = Precipitation Seasonality (Coefficient of Variation) "PS"
## BIO16 = Precipitation of Wettest Quarter "PWeQ"
## BIO17 = Precipitation of Driest Quarter "PDQ"
## BIO18 = Precipitation of Warmest Quarter "PWaQ"
## BIO19 = Precipitation of Coldest Quarter "PCQ"

coords <- data.frame(x=all.geo[,"lon"],y=all.geo[,"lat"])
points <- SpatialPoints(coords, proj4string = wc@crs)
values <- extract(wc,points)
df <- cbind.data.frame(coordinates(points),values)
rownames(df) <- all.geo[,"Species.name"]
colnames(df)[1:2] <- c("Lon","Lat")
### Bioclim imports temps as integers by default
### They just multiply by 10, so we divide to get the true temps
df[,grep("T",colnames(df))] <- df[,grep("T",colnames(df))] / 10

### Just Ap
apg.bio <- df[grep("Aphaenogaster",rownames(df)),]
apg.mash <- ncbi.gen[grep("Aphaenogaster",rownames(df)),grep("Aphaenogaster",rownames(df))]
all(rownames(apg.bio) == rownames(apg.mash))

bio.i <- 2
pdf(file = "tmp.pdf")
plot(r[[bio.i]],main = names(wc)[bio.i])
plot(points, add = TRUE, pch = 19, cex = 0.5)
dev.off()
system("scp tmp.pdf matthewklau@fas.harvard.edu:public_html")

ord <- nmds(as.dist(all.mash),3,3, nits = 500)
nms <- nmds.min(ord,3)
apg.ord <- nmds(as.dist(apg.mash),2,2, nits = 500)
apg.nms <- nmds.min(apg.ord, 2)
write.csv(apg.nms, "data/storage/apg/nmds_apg.csv", row.names = FALSE)
write.csv(nms, "data/storage/apg/nmds_all.csv", row.names = FALSE)

vec <- envfit(nms, df)


pdf(file = "tmp.pdf")
plot(nms[,1:2], pch = "")
text(nms[,1:2], labels = sapply(rownames(all.mash), get_names), cex = 0.5)
plot(vec, col = "darkgrey", cex = 0.5)
dev.off()
system("scp tmp.pdf matthewklau@fas.harvard.edu:public_html")


apg.vec <- envfit(apg.nms, apg.bio)
apg.vec

vec.out <- cbind(r = (vec[["vectors"]][["r"]]), p = vec[["vectors"]][["pvals"]])
apg.vec.out <- cbind(r = (apg.vec[["vectors"]][["r"]]), p = apg.vec[["vectors"]][["pvals"]])

