library(raster)
library(sp)

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
library(ggmap)
ncbi_info <- read.csv("../data/gaga_genome_info.csv")
geocode("Boston, MA")
