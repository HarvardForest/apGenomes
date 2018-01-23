library(raster)
library(sp)
library(ggmap)

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
lats <- clim.data[,"lat"]
lons <- clim.data[,"long"]
coords <- data.frame(x=lons,y=lats)
points <- SpatialPoints(coords, proj4string = wc@crs)
values <- extract(wc,points)
df <- cbind.data.frame(coordinates(points),values)
rownames(df) <- rownames(clim.data)
### Bioclim imports temps as integers by default
### They just multiply by 10, so we divide to get the true temps
df[,grep("T",colnames(df))] <- df[,grep("T",colnames(df))] / 10

values <- extract(r,points)
df <- cbind.data.frame(coordinates(points),values)

pdf(file = "tmp.pdf")
heatmap(abs(cor(cbind(clim.data[,-3],df))))
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
