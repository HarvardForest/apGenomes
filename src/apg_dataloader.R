pkg <- c("gdata","prism","ggplot2","raster","AntWeb")
sapply(pkg,function(pkg) 
    if(!pkg %in% installed.packages()[,1]){
        install.packages(pkg)
    }
)

library(gdata)
library(prism)
library(ggplot2)
library(raster)
library(AntWeb)
library(geosphere)
library(rnoaa)

get.prism <- FALSE

## sample information
broad.info <- read.csv('data/storage/apg/SSF-1728_KeyforCollaborator.csv')
sample.info <- read.csv('data/storage/apg/colony_locations.csv')
ant.info <- read.csv('data/storage/apg/RADseq_mastersheet_2014.csv')
ant.info <- ant.info[ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID. %in% na.omit(sample.info$spec_epithet),]
ant.info <- ant.info[!ant.info$State == "",]
ant.geo <- read.csv('data/storage/apg/ant_sites.csv')
ant.geo[,c("Lon","Lat")] <- apply(ant.geo[,c("Lon","Lat")],2,as.numeric)
ant.info <- data.frame(ant.info,
                       ant.geo[match(as.character(ant.info$Locale),ant.geo$Site),c("Lon","Lat")])
ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID. <- 
    as.character(ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID.)
ant.info <- na.omit(ant.info)
ant.color <- ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID.
ant.factor <- factor(ant.color)



## geographic information
geo.ctr <- split(ant.info[,c("Lon","Lat")],ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID.)
geo.ctr <- lapply(geo.ctr,function(x) apply(x,2,mean))
geo.ctr <- do.call(rbind,geo.ctr)

ap.ctr <- do.call(rbind,list(rud1 = geo.ctr[6,],
                             rud6 = geo.ctr[6,],
                             pic1 = geo.ctr[5,],
                             mia1 = geo.ctr[4,],
                             ful1 = geo.ctr[3,],
                             ash1 = geo.ctr[1,],
                             flo1 = geo.ctr[2,]))

### apg sample geographic info
### arudis and picea are from google earth
### All other coords are from ant_sites.csv
apg.geo  <- do.call(rbind,list('arudis1' = c(-78.9830464,36.0200847),
                 'rud6' = c(-78.9830464,36.0200847),
                 'pic1' = c(-72.5847494,42.6004513),
                 'mia1' = c(-82.301773,29.657955),
                 'ful1' = c(-82.514575,32.692384),
                 'ash1' = c(-82.031176,29.785325),
                 'flo1' = c(-82.031176,29.785325)))
colnames(apg.geo) <- c('Longitude','Latitude')
apg.geo.labs <- paste0(substr(rownames(apg.geo),1,3),
                       substr(rownames(apg.geo),nchar(rownames(apg.geo)),
                              nchar(rownames(apg.geo))))
apg.geo.labs[1] <- 'rud1'

### geographic distance for samples
apg.gd <- array(NA,dim = rep(nrow(apg.geo),2))
rownames(apg.gd) <- colnames(apg.gd) <- rownames(apg.geo)
for (i in 1:nrow(apg.geo)){
    for (j in 1:nrow(apg.geo)){
        apg.gd[i,j] <- distm (apg.geo[i,], apg.geo[j,], 
                               fun = distHaversine)
    }
}

apg.gcd <- array(NA,dim = rep(nrow(ap.ctr),2))
rownames(apg.gcd) <- colnames(apg.gcd) <- rownames(ap.ctr)
for (i in 1:nrow(ap.ctr)){
    for (j in 1:nrow(ap.ctr)){
        apg.gcd[i,j] <- distm (ap.ctr[i,], ap.ctr[j,], 
                               fun = distHaversine)
    }
}


### 
if (get.prism){
    options(prism.path = "~/prismtmp")
    if (!'jnorm' %in% ls()){
        get_prism_monthlys(type="ppt", year = 2000:2016, mon = 1:12, keepZip=T)
        get_prism_monthlys(type="tmax", year = 2000:2016, mon = 1:12, keepZip=T)
        get_prism_monthlys(type="tmin", year = 2000:2016, mon = 1:12, keepZip=T)
    }
    abs.path <- ls_prism_data(absPath = T)[,2]
    
    jul.max <- grep('tmax_30yr_normal',abs.path,value = TRUE)
    jul.max <- grep('07_bil',jul.max,value = TRUE)
    jul.max <- raster(jul.max)
    jan.min <- grep('tmin_30yr_normal',abs.path,value = TRUE)
    jan.min <- grep('01_bil',jan.min,value = TRUE)
    jan.min <- raster(jan.min)

### crop rasters
### extent: -125.0208, -66.47917, 24.0625, 49.9375  (xmin, xmax, ymin, ymax)
    crop.ext <- extent(-85,-65,20,49)
    jan.min <- crop(jan.min,crop.ext)
    jul.max <- crop(jul.max,crop.ext)
}


