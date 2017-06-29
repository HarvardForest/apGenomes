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
library(gdata)

### Analysis Outline

## sample information
broad.info <- read.csv('data/storage/apg/broad_sample_key.csv')
sample.info <- read.csv('data/colony_locations.csv')
ant.info <- read.csv('data/RADseq_mastersheet_2014.csv')
ant.info <- ant.info[ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID. %in% na.omit(sample.info$spec_epithet),]
ant.info <- ant.info[!ant.info$State == "",]
ant.geo <- read.csv('data/ant_sites.csv')
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

## site <- ant.info[1,c("Lon","Lat")]
## get_prism_monthlys(type="tmean", year = 1982:2014, mon = 1, keepZip=F)
## to_slice <- grep("_[0-9]{4}[0][1]",ls_prism_data()[,1],value=T)
## to_slice <- grep("tmean",to_slice, value = T)
## p <- prism_slice(as.numeric(site[1,]),to_slice)
## p + stat_smooth(method="lm",se=F) + theme_bw() +
##     ggtitle(paste0("Avg Jan Temp",ant.geo[1,1]))
jnorm <- raster(ls_prism_data(absPath=T)[1,2])
## j2013 <- raster(ls_prism_data(absPath=T)[52,2])

plot(jnorm)
points(geo.ctr,pch=19,cex=0.25,col=1:5)

plot(geo.ctr,pch=19,cex=2,col=1:5,xlim = c(-85,-70), ylim = c(30,45))
text(geo.ctr,labels = rownames(geo.ctr),pos=4)

## cross reference site-collection with location
phyto.info <- read.csv('data/storage/apg/Phytotron\ colonies\ 2013\ Transcriptome.csv')

## gaemr info
gaemr.tab <- read.csv('data/gaemr-table.csv')
metrics <- c("Mean Total Aligned Depth","Total Contig Length","Genome size estimate","Contig N50","Scaffold N50","Assembly GC","SNP rate                   ~","Snps")
gaemr.tab <- gaemr.tab[gaemr.tab$Metric %in% metrics,]
gaemr.tab$Value <- sub(' bases','',gaemr.tab$Value)
for (i in 1:nrow(gaemr.tab)){
    if (grepl('SNP rate',gaemr.tab$Metric[i])){
        x <- as.numeric(strsplit(gaemr.tab$Value[i],'\\/')[[1]])
        gaemr.tab$Value[i] <- round((x[1] / x[2]), 7)
    }
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


