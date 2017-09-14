### Get prism data for coordinates
### This is borrowed code:
### http://eremrah.com/articles/How-to-extract-data-from-PRISM-raster/

library(tidyr)
library(stringr)
library(prism)
library(raster)
library(magrittr)

## SEE: https://github.com/ropensci/prism
## Get PRISM data
if (2 > 3){system("rm -rf ~/prismtmp")}
options(prism.path = "~/prismtmp")
# get_prism_dailys(type="tmin", minDate = "1991-01-01", maxDate = "1991-01-31", keepZip = FALSE)
# get_prism_dailys(type="tmax", minDate = "1991-01-01", maxDate = "1991-01-31", keepZip = FALSE)

get_prism_normals(type="ppt", "800m", annual = TRUE)
get_prism_normals(type="tmin", "800m", annual = TRUE)
get_prism_normals(type="tmax", "800m", annual = TRUE)

## Stack files
mystack <- ls_prism_data() %>%  prism_stack()  

## Get proj from raster stack
mycrs <- mystack@crs@projargs

## My points
mypoints <- data.frame(id = rownames(apg.geo),
                       lat = apg.geo[,"Latitude"],
                       long = apg.geo[,"Longitude"]
)

## Convert points to spatial points data frame
coordinates(mypoints) <- c('long', 'lat')
proj4string(mypoints) <- CRS(mycrs)

## Extract data from raster
data <- data.frame(coordinates(mypoints), mypoints$id,
                   extract(mystack, mypoints))

## Rename column headers
colnames(data)[4:6] <- do.call(rbind,strsplit(colnames(data)[4:6],"_"))[,2]
rownames(data)[rownames(data) == "arudis1"] <- "rud1"
data[,"mypoints.id"] <- as.character(data[,"mypoints.id"])
data[rownames(data) == "rud1","mypoints.id"] <- "rud1"

## Rename object
clim.data <- data[match(c("pic1", "rud1", "rud6", "ful1", "flo1", "mia1", "ash1"),rownames(data)),]

## create a climate and distance objects
clim.d <- dist(apply(clim.data[,c("ppt","tmax","tmin")],2,function(x) (x - mean(x)) / sd(x)))
temp.d <- dist(apply(clim.data[,c("tmax","tmin")],2,function(x) (x - mean(x)) / sd(x)))
ppt.d <- as.matrix(dist(clim.data[,c("ppt")]))
rownames(ppt.d) <- colnames(ppt.d) <- rownames(clim.data)
ppt.d <- as.dist(ppt.d)

tmax.d <- as.matrix(dist(clim.data[,c("tmax")]))
rownames(tmax.d) <- colnames(tmax.d) <- rownames(clim.data)
tmax.d <- as.dist(tmax.d)

tmin.d <- as.matrix(dist(clim.data[,c("tmin")]))
rownames(tmin.d) <- colnames(tmin.d) <- rownames(clim.data)
tmin.d <- as.dist(tmin.d)
