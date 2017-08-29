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

