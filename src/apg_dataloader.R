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

## sample information
broad.info <- read.xls('/Volumes/ellison_lab/ap_genomes/SSF-1728_KeyforCollaborator.xlsx',2)
sample.info <- read.csv('../docs/colony_locations.csv')
ant.info <- read.csv('../docs/RADseq_mastersheet_2014.csv')
ant.info <- ant.info[ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID. %in% na.omit(sample.info$spec_epithet),]
ant.info <- ant.info[!ant.info$State == "",]
ant.geo <- read.csv('../docs/ant_sites.csv')
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

### 
geo.apg <- 

options(prism.path = "~/prismtmp")
if (!'jnorm' %in% ls()){
    get_prism_monthlys(type="tmin", year = 2000:2016, mon = 7, keepZip=F)
    get_prism_normals(type = "tmin",mon = 1,res = '4km',keepZip = F)
    get_prism_normals(type = "tmax",mon = 7,res = '4km',keepZip = F)
}
abs.path <- ls_prism_data(absPath = T)[,2]

jul.max <- grep('tmax_30yr_normal',abs.path,value = TRUE)
jul.max <- grep('07_bil',jul.max,value = TRUE)
jul.max <- raster(jul.max)

jan.min <- grep('tmin_30yr_normal',abs.path,value = TRUE)
jan.min <- grep('01_bil',jan.min,value = TRUE)
jan.min <- raster(jan.min)












