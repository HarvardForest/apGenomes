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


### Nat History

##############################
#### Outline of results ######
##############################

## Get contaminent percentages


## Table. GEAMR output
## Figure. Sequenced ants and genomes
## Figure. BLAST hits to other ants 
## Figure. Genome sizes of other ants
## Figure. Bivariate plot of climate variables vs Size and GC

source('apg_dataloader.R')
library("ggmap")
library(maptools)
library(maps)

pdf('../docs/manuscript/map_apg.pdf',width = 5,height = )
map("state",
    fill=TRUE, col="white", bg="white",xlim = c(-85,-65),mar = c(0,0,0,0),asp = 1.25)
points(apg.geo,cex = 2,pch = 19,col = 'darkgrey')
###text(apg.geo,labels = rownames(apg.geo),cex = 0.45)
dev.off()



### QC, Composition, Structure
### GAEMR Info
if (!'table.stats' %in% ls()){source('assembly_info.R')}

xtable::xtable(table.stats)
names(xtable

print(,
      type = "latex",
      file = "../docs/manuscript/assembly_stats.tex",
      include.rownames = FALSE
      )

### Compare to other ant genomes and bees
ref.sizes <- c(as.numeric(ncbi.ant[,'Size (Mb)']),ant.gen.size[,'1C Genome Size (Mb)'])
if (!ncbi.ant %in% ls()){source('apg_genome_compare.R')}

stats.col <- c(6,6,5,4,3,1,2)
par(mfrow = c(1,2))
hist(ref.sizes,main = '',xlab = 'Size (Mb)')
abline(v = stats[,'TotalScaffoldLength'] / (1000000),col = stats.col)
hist(as.numeric(ncbi.ant[,'GC%']),main = '',xlab = 'GC%',xlim = gc.xlim)
abline(v = stats[,'AssemblyGC'],col = stats.col)


par(mfrow=c(1,1),mar = c(5.1,5.1,1.1,1.1),cex.lab = 2)
plot(log10(as.numeric(ncbi.insect[,'Genes']))~log10(as.numeric(ncbi.insect[,'Size (Mb)'])),
     xlab = expression('log'[10]*' Genome Size (Mb)'),ylab = expression('log'[10]*' Gene Size (bp)'),
     pch = 19, col = 'darkgrey')
abline(lm(log10(as.numeric(ncbi.insect[,'Genes']))~log10(as.numeric(ncbi.insect[,'Size (Mb)']))))
points(log10(stats[,'TotalScaffoldLength'] / (1000000)),log10(stats[,'MaxGap']),col = stats.col,pch = 1,cex = 1.5,lwd = 3)

legend('topleft',legend = c('Gene Size','Scaffold Size'), pch = c(19,1), 
       col = c('darkgrey','black'))

mtext(side = 4, line = 3, 'log10(Scaffold Size (bp))')

### Genomic biogeography = lat/lon, distance, climate = temperature, climatic similarities

### Blast Search for Transcriptome targets
## 1. Heat stress response genes @Stanton-Geddes2016
## 2. Cold tolerance genes @Stanton-Geddes2016
## 3. Map transcripts onto genome and look for inter-specific variance @Stanton-Geddes2016

if (!'ant.gen.size' %in% ls()){source('ant_genome_size.R')}

cor.test(stats[,'Removed.Scaffold'],apg.geo[,'Latitude'])
cor.test(stats[,'TotalScaffoldLength'],apg.geo[,'Latitude'])
cor.test(stats[,'Removed.Scaffold'],ap.ctr[,'Lat'])
cor.test(stats[,'TotalScaffoldLength'],ap.ctr[,'Lat'])

size <- as.numeric(c((stats[,'TotalScaffoldLength']/1000000),ncbi.ant[match(rownames(ncbi.geo),ncbi.ant[,1]),'Size (Mb)']))
lat <- c(apg.geo[,'Latitude'],ncbi.geo[,'Latitude'])





cor.test(size,lat)
plot(lat,size)

par(mfrow = c(2,2))
plot(stats[,'Removed.Scaffold']~apg.geo[,'Latitude'],
     xlab = 'Latitude',ylab = 'Contaminant Scaffold')
abline(lm(stats[,'Removed.Scaffold']~apg.geo[,'Latitude']))
plot(stats[,'TotalScaffoldLength']~apg.geo[,'Latitude'],
     xlab = 'Latitude',ylab = 'Total Scaffold Length')
plot(stats[,'Removed.Scaffold']~ap.ctr[,'Lat'],
     xlab = 'Lat',ylab = 'Contaminant Scaffold')
plot(stats[,'TotalScaffoldLength']~ap.ctr[,'Lat'],
     xlab = 'Latitude',ylab = 'Total Scaffold Length')



## cross reference site-collection with location
phyto.info <- read.xls('/Users/hermes/Dropbox/WarmAntDimensions/Phytotron\ 2013/Phytotron\ colonies\ 2013\ Transcriptome.xlsx',1)

## phytotron information
radseq.info <- read.xls('~/Dropbox/WarmAntDimensions/Genomics/Data_allocation_mastersheet_12-23-16.xlsx')
phyt.info <- read.xls('~/Dropbox/WarmAntDimensions/Phytotron 2013/Phytotron colonies 2013 Transcriptome.xlsx',3)
## phyt.info <- read.xls('~/Dropbox/WarmAntDimensions/Phytotron 2013/Aphaeno thermal data_field and phyto.xls')


