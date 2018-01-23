library(raster)
library(sp)

"data/gaga_genome_info.csv"

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

txtplot(clim.data[,"tmax"],df[,"Tmax"])
clim.data[,"tmax"] - df[,"Tmax"]

range(cor(cbind(clim.data[,-3],df)))
range(abs(cor(cbind(clim.data[,-3],df))))

pdf(file = "tmp.pdf")
heatmap(abs(cor(cbind(clim.data[,-3],df))))
dev.off()
system("scp tmp.pdf matthewklau@fas.harvard.edu:public_html")

bc.d <- dist((df[,-1:-2]))
mantel(bc.d~geo.d)
mantel(mash.d~bc.d)
mantel(mash.d~geo.d)
