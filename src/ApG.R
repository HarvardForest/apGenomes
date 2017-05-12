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
source('xtable_helpers.R')
source('mumplot.R')
source('mash_helpers.R')

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

## Figure: map_apg.jpg
pdf('../docs/manuscript/map_apg.pdf',width = 5,height = )
map("state",
    fill=TRUE, col="white", bg="white",xlim = c(-85,-65),mar = c(0,0,0,0),asp = 1.25)
points(apg.geo,cex = 2,pch = 19,col = 'darkgrey')
###text(apg.geo,labels = rownames(apg.geo),cex = 0.45)
dev.off()


## %%%%%%%%%%%%%%%%%
## % Antweb photos %
## %%%%%%%%%%%%%%%%%

## % Citing Antweb -> https://www.antweb.org/citing_antweb.jsp
## % ashmeadi https://www.antweb.org/bigPicture.do?name=casent0103563&shot=p&number=1
## % fordidiana https://www.antweb.org/bigPicture.do?name=casent0103583&shot=p&number=1
## % fulva https://www.antweb.org/bigPicture.do?name=casent0103585&shot=p&number=1
## % miamiana https://www.antweb.org/bigPicture.do?name=casent0103595&shot=p&number=1
## % picea https://www.antweb.org/bigPicture.do?name=casent0104844&shot=p&number=1
## % rudis https://www.antweb.org/bigPicture.do?name=casent0104843&shot=p&number=1


### QC, Composition, Structure
### GAEMR Info
if (!'table.stats' %in% ls()){source('assembly_info.R')}

### Formatting
apg.xtab <- xtable::xtable(t(table.stats[,-1]))

## Distibuish rudis colonies
names(apg.xtab) <- paste0(table.stats[,1],c(1,2,rep('',5)))

### Table: assembly_stats
print(apg.xtab,
      type = "latex",
      file = "../docs/manuscript/assembly_stats.tex",
      sanitize.colnames.function = italic,
      include.rownames = TRUE,
      include.colnames = TRUE
      )

### Mash dist
msh.d <- read.table('~/storage/ap_genomes/apg_mashdist.txt')
msh.d[,1:2] <- apply(msh.d[,1:2],2,substr,start = 4,stop = 8)
msh.d <- as.mashdist(msh.d)
rownames(msh.d) <- colnames(msh.d) <- rownames(apg.gd)

mash.mantel <- mantel(apg.gd,msh.d)
mash.mantel

mash.mantel <- mantel(apg.gd[-2,-2],msh.d[-2,-2])
mash.mantel

### Figure: clustering distance and mash
png('../docs/manuscript/clustdist.png',width = 900)
par(mfrow = c(1,2))
plot(hclust(as.dist(apg.gd)),main = '',xlab = 'Geographic Distance')
plot(hclust(as.dist(msh.d)),main = '',xlab = 'Genomic Distance (mash)')
dev.off()

### Mummer 
### system(paste("awk 'f{print;f=0} />/{f=1}'",file,"> .tmp"))
mum.thresh <- 15
sd.thresh <- 0
files <- grep('dists',
              dir('~/storage/ap_genomes/mums',
                  full = TRUE),val = TRUE)
if (!"pos.l" %in% ls()){
    pos.l <- lapply(files,function(x) 
        do.call(rbind,strsplit(readLines(x),' ')))
    pos.l <- lapply(pos.l,function(x)
        apply(x,2,as.numeric))
}
mum.ld <- lapply(pos.l,function(x) sum(abs(apply(x[,1:2],1,diff))))
mum.d <- unlist(mum.ld)
dist.pic <- apg.geo[,'Latitude'] - apg.geo[,'Latitude']['pic1']

### Figure: mum dissimilarity from picea 
par(mfrow=c(1,1))
png('../docs/manuscript/mumdispic.png')
plot(mum.d~dist.pic,pch = '',cex = 1,
     xlab = 'Latitudinal Distance from A. picea',
     ylab = 'Cumulative Change in Position of MUMs')
abline(lm(mum.d~dist.pic))
text(mum.d~dist.pic,labels=names(dist.pic))
dev.off()

### Figure: mum plots
png('../docs/manuscript/mumspark.png',
    width = 1000,height = 1500,
    pointsize = 18)
par(mfrow = c(3,2),cex.lab = 2, cex.axis = 2.0,mar = c(4.5,5.1,4.1,2.1))
for (i in ((1:nrow(table.stats))[-3])){
    print(files[i])
    mumspark(files[i],
        xlab = paste(as.character(table.stats[3,'Species']),
                     'MUM Position (bp)'),
        ylab = paste(as.character(table.stats[,'Species'])[i],
                     'MUM Position Change (bp)'),main = '',cex = 1)
}
dev.off()

### Compare mum and mash
png('../docs/manuscript/mum_mash.png',width = 500, height = 500)
plot(msh.d['pic1',],mum.d,pch = '',
     xlab = 'Mash Distance',
     ylab = 'MUM Cumulative Change')
text(msh.d['pic1',],mum.d,labels = rownames(msh.d))
dev.off()

### Compare to other ant genomes and bees
if (!ncbi.ant %in% ls()){source('apg_genome_compare.R')}
ref.sizes <- c(as.numeric(ncbi.ant[,'Size (Mb)']),ant.gen.size[,'1C Genome Size (Mb)'])

## Figure: hist_ncbi
stats.col <- c(6,6,5,4,3,1,2)
png('../docs/manuscript/hist_ncbi.png',width = 1000,height = 500)
par(mfrow = c(1,2))
hist(ref.sizes,main = '',xlab = 'Assembly Scaffold Size (Mb)')
abline(v = stats[,'TotalScaffoldLength'] / (1000000),col = stats.col,lwd = 2)
hist(as.numeric(ncbi.ant[,'GC%']),main = '',xlab = 'GC%',xlim = gc.xlim)
abline(v = stats[,'AssemblyGC'],col = stats.col,lwd = 2)
dev.off()

## Figure: genome_sizes
png('../docs/manuscript/genome_sizes.png',width = 500, height = 500)
par(mfrow=c(1,1),mar = c(5.1,5.1,1.1,1.1),cex.lab = 2)
plot(log10(as.numeric(ncbi.insect[,'Genes']))~log10(as.numeric(ncbi.insect[,'Size (Mb)'])),
     xlab = expression('log'[10]*' Genome Size (Mb)'),
     ylab = expression('log'[10]*' Gene or Scaffold Size (bp)'),
     pch = 19, col = 'darkgrey')
abline(lm(log10(as.numeric(ncbi.insect[,'Genes']))~log10(as.numeric(ncbi.insect[,'Size (Mb)']))))
points(log10(as.numeric(ncbi.ant[,'Size (Mb)'])),log10(as.numeric(ncbi.ant[,'Genes'])),
       col = 'black',pch = 19)
points(log10(stats[,'TotalScaffoldLength'] / (1000000)),log10(stats[,'MaxGap']),
       col = stats.col,pch = 1,cex = 1.5,lwd = 3)
legend('topleft',legend = c('Insect Gene Size','Ant Gene Size','Assembly Scaffold Size'), 
       pch = c(19,19,1),  col = c('black','darkgrey','black'))
dev.off()

## NCBI ant information
ncbi.xtab <- ncbi.ant[,c('BioProject Accession','Release Date')]
rownames(ncbi.xtab) <- ncbi.ant[,1]
ncbi.xtab <- xtable::xtable(ncbi.xtab)

## Table: create ncbi_ants 
print(ncbi.xtab,
      type = "latex",
      file = "../docs/manuscript/ncbi_ants.tex",
      sanitize.rownames.function = italic,
      include.rownames = TRUE,
      include.colnames = TRUE
      )


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


## Ant genome releases by year
plot(table(substr(ncbi.ant[,c('Release Date')],1,4)),xlab = 'Year',ylab = 'Genomes')


## cross reference site-collection with location
phyto.info <- read.xls('/Users/hermes/Dropbox/WarmAntDimensions/Phytotron\ 2013/Phytotron\ colonies\ 2013\ Transcriptome.xlsx',1)

## phytotron information
radseq.info <- read.xls('~/Dropbox/WarmAntDimensions/Genomics/Data_allocation_mastersheet_12-23-16.xlsx')
phyt.info <- read.xls('~/Dropbox/WarmAntDimensions/Phytotron 2013/Phytotron colonies 2013 Transcriptome.xlsx',3)
## phyt.info <- read.xls('~/Dropbox/WarmAntDimensions/Phytotron 2013/Aphaeno thermal data_field and phyto.xls')


