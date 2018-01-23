library(raster)
library(sp)
library(ggmap)

r <- getData("worldclim",var="bio",res=10)
r <- r[[c(1,12)]]
names(r) <- c("Temp","Prec")
lats <- c(9.093028 , 9.396111, 9.161417)
lons <- c(-11.7235, -11.72975, -11.709417)
coords <- data.frame(x=lons,y=lats)
points <- SpatialPoints(coords, proj4string = r@crs)

values <- extract(r,points)
df <- cbind.data.frame(coordinates(points),values)

pdf(file = "tmp.pdf")
plot(r[[1]])
plot(points, add = TRUE, pch = ".")
dev.off()
system("scp tmp.pdf matthewklau@fas.harvard.edu:public_html")


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

