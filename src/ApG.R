### Analytical script for the Apaenogaster genome
### MKLau 07February2017
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

### Analysis Outline

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

geo.ctr <- split(ant.info[,c("Lon","Lat")],ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID.)
geo.ctr <- lapply(geo.ctr,function(x) apply(x,2,mean))
geo.ctr <- do.call(rbind,geo.ctr)

site <- ant.info[1,c("Lon","Lat")]
if (!'jnorm' %in% ls()){
    get_prism_monthlys(type="tmean", year = 2000:2016, mon = 1, keepZip=F)
}
to_slice <- grep("_[0-9]{4}[0][1]",ls_prism_data()[,1],value=T)
to_slice <- grep("tmean",to_slice, value = T)
jnorm <- raster(ls_prism_data(absPath=T)[1,2])
p <- prism_slice(as.numeric(site[1,]),to_slice)
p + stat_smooth(method="lm",se=F) + theme_bw() +
    ggtitle(paste0("Avg Jan Temp",ant.geo[1,1]))

par(mfrow = c(1,2))
plot(jnorm)
points(geo.ctr,pch=19,cex=0.25,col=1:6)
plot(geo.ctr,pch=19,cex=2,col=1:6,xlim = c(-85,-70), ylim = c(30,45))
text(geo.ctr,labels = rownames(geo.ctr),pos=4,cex=0.75)

### GAEMR Info
source('assembly_info.R')

### Make the table
print(xtable::xtable(stats),
      type = "latex",
      file = "../docs/manuscript/assembly_stats.tex")

### NCBI Genome Info
source('ncbi_genome_info.R')
colnames(ncbi.ant)[1] <- 'Organism'

### AntWeb Info
aw.apg <- list()
for (i in 1:nrow(sample.info)){
    aw.apg[[i]] <- aw_data(scientific_name = paste('Aphaenogaster',sample.info[i,'spec_epithet']),georeferenced = TRUE)
}
aw.ncbi <- list()
for (i in 1:nrow(sample.info)){
    aw.ncbi[[i]] <- aw_data(scientific_name = ncbi.ant[i,1],georeferenced = TRUE)
}


### Inter-species comparisons
size.xlim <- range(as.numeric(c((stats[,'TotalScaffoldLength'] / (1000000)),ncbi.ant[,'Size (Mb)'])))
size.xlim <- c(floor(size.xlim[1]),ceiling(size.xlim[2]))
gc.xlim <- range(as.numeric(c(stats[,'AssemblyGC'],ncbi.ant[,'GC%'])))
gc.xlim <- c(floor(gc.xlim[1]),ceiling(gc.xlim[2]))

stats.col <- c(6,6,5,4,3,1,2)
par(mfrow = c(1,2))
hist(as.numeric(ncbi.ant[,'Size (Mb)']),main = '',xlab = 'Size (Mb)')
abline(v = stats[,'TotalScaffoldLength'] / (1000000),col = stats.col)
hist(as.numeric(ncbi.ant[,'GC%']),main = '',xlab = 'GC%',xlim = gc.xlim)
abline(v = stats[,'AssemblyGC'],col = stats.col)


plot(log(as.numeric(ncbi.ant[,'Genes']))~log(as.numeric(ncbi.ant[,'Size (Mb)'])))
abline(lm(log(as.numeric(ncbi.ant[,'Genes']))~log(as.numeric(ncbi.ant[,'Size (Mb)']))))

ncbi.ant[ncbi.ant[,'Genes'] == max(ncbi.ant[,'Genes']),]

## cross reference site-collection with location
phyto.info <- read.xls('/Users/hermes/Dropbox/WarmAntDimensions/Phytotron\ 2013/Phytotron\ colonies\ 2013\ Transcriptome.xlsx',1)

## phytotron information
radseq.info <- read.xls('~/Dropbox/WarmAntDimensions/Genomics/Data_allocation_mastersheet_12-23-16.xlsx')
phyt.info <- read.xls('~/Dropbox/WarmAntDimensions/Phytotron 2013/Phytotron colonies 2013 Transcriptome.xlsx',3)
## phyt.info <- read.xls('~/Dropbox/WarmAntDimensions/Phytotron 2013/Aphaeno thermal data_field and phyto.xls')

## geographic information


## 0. genome quality (from the gaemer results)

## GC content

## apg.dat

## par(mfrow=c(2,2))
## plot(Assembly.GC~Lat,data=apg.dat)
## abline(lm(Assembly.GC~Lat,data=apg.dat))
## plot(Genome.size.estimate~Lat,data=apg.dat)
## abline(lm(Genome.size.estimate~Lat,data=apg.dat))
## plot(Scaffold.N50~Lat,data=apg.dat)
## abline(lm(Scaffold.N50~Lat,data=apg.dat))
## plot(SNP.rate....................~Lat,data=apg.dat)
## abline(lm(SNP.rate....................~Lat,data=apg.dat))

## summary(lm(Assembly.GC~Lat,data=apg.dat))
## summary(lm(Genome.size.estimate~Lat,data=apg.dat))
## summary(lm(Contig.N50~Lat,data=apg.dat))
## summary(lm(Scaffold.N50~Lat,data=apg.dat))
## summary(lm(SNP.rate....................~Lat,data=apg.dat))

## round(na.omit(apply(apply(apg.dat,2,as.numeric),2,mean)),3)

## apg.dat
## percent contaminants
## 

## 1. Heat stress response genes @Stanton-Geddes2016


## 2. Cold tolerance genes @Stanton-Geddes2016


## 3. Map transcripts onto genome and look for inter-specific variance @Stanton-Geddes2016


## 4. Immunity @Roux2014


## 5. Olfaction @Roux2014


## 6. Mitochondrial function @Roux2014


## 7. Social chromosome? @Brelsford2014
