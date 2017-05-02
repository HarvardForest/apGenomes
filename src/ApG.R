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

pdf('../docs/manuscript/map_janmin.pdf',width = 5,height = )
plot(jan.min,main = 'January 30 Year Minimum (C)',xlim = c(-85,-66),asp = 1.25)
points(jitter(apg.geo,10),cex = 1.45,pch = 1)
dev.off()

plot(jul.max,main = 'July 30 year maximum temp')

### Add sample locations

### QC, Composition, Structure
### GAEMR Info
if (!'stats' %in% ls()){source('assembly_info.R')}

(table.stats - table.uf.stats)[,'Total Scaffold Length']
print(xtable::xtable(table.stats),
      type = "latex",
      file = "../docs/manuscript/assembly_stats.tex")

### Compare to other ant genomes and bees
if (!ncbi.ant %in% ls()){source('apg_genome_compare.R')}
stats.col <- c(6,6,5,4,3,1,2)
par(mfrow = c(1,2))
hist(as.numeric(ncbi.ant[,'Size (Mb)']),main = '',xlab = 'Size (Mb)')
abline(v = stats[,'TotalScaffoldLength'] / (1000000),col = stats.col)
hist(as.numeric(ncbi.ant[,'GC%']),main = '',xlab = 'GC%',xlim = gc.xlim)
abline(v = stats[,'AssemblyGC'],col = stats.col)

par(mfrow=c(1,1))
plot(log(as.numeric(ncbi.ant[,'Genes']))~log(as.numeric(ncbi.ant[,'Size (Mb)'])),
     xlab = 'log(Genome Size (Mb))',ylab = 'log(Gene Number)')
abline(lm(log(as.numeric(ncbi.ant[,'Genes']))~log(as.numeric(ncbi.ant[,'Size (Mb)']))),lty =2)
abline(v = log(stats[,'TotalScaffoldLength'] / (1000000)),col = stats.col)


### Genomic biogeography = lat/lon, distance, climate = temperature, climatic similarities

### Blast Search for Transcriptome targets
## 1. Heat stress response genes @Stanton-Geddes2016
## 2. Cold tolerance genes @Stanton-Geddes2016
## 3. Map transcripts onto genome and look for inter-specific variance @Stanton-Geddes2016

cor.test(stats[,'Removed.Scaffold'],apg.geo[,'Latitude'])
cor.test(stats[,'TotalScaffoldLength'],apg.geo[,'Latitude'])
cor.test(stats[,'Removed.Scaffold'],ap.ctr[,'Lat'])
cor.test(stats[,'TotalScaffoldLength'],ap.ctr[,'Lat'])


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


