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

## Get contiminent percentages

## Table. GEAMR output
## Figure. Sequenced ants and genomes
## Figure. BLAST hits to other ants 
## Figure. Genome sizes of other ants
## Figure. Bivariate plot of climate variables vs Size and GC

source('apg_dataloader.R')

site <- ant.info[1,c("Lon","Lat")]

par(mfrow = c(2,2))
plot(jmin,main = 'January Minimum')
plot(jnorm,main = 'January Average')
plot(julmin,main = 'July Minimum')
plot(julmin2010,main = 'July Minimum')

par(mfrow = c(1,1))
plot(jan.min,main = 'January 30 year minimum temp')

plot(jul.max,main = 'July 30 year maximum temp')

### Add sample locations

### QC, Composition, Structure
### GAEMR Info
if (!'stats' %in% ls()){source('assembly_info.R')}
table.stats <- stats[,c(
    grep('GC',colnames(stats)),
    grep('gap',colnames(stats),ign = TRUE),
    grep('contig',colnames(stats),ign = TRUE),
    grep('Scaffold',colnames(stats))
    )]
tab.names <- c("AssemblyGC","TotalGapLength","CapturedGaps",
"Contigs","MaxContig","MeanContig","ContigN50","ContigN90","TotalContigLength",
"Scaffolds","MaxScaffold","MeanScaffold","ScaffoldN50","ScaffoldN90","TotalScaffoldLength")
table.stats <- table.stats[,match(tab.names,colnames(table.stats))]
colnames(table.stats) <- c("Assembly GC","Total Gap Length","Captured Gaps",
"Contigs","Max Contig","Mean Contig","Contig N50","Contig N90","Total Contig Length",
"Scaffolds","Max Scaffold","Mean Scaffold","Scaffold N50","Scaffold N90","Total Scaffold Length")

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

## cross reference site-collection with location
phyto.info <- read.xls('/Users/hermes/Dropbox/WarmAntDimensions/Phytotron\ 2013/Phytotron\ colonies\ 2013\ Transcriptome.xlsx',1)

## phytotron information
radseq.info <- read.xls('~/Dropbox/WarmAntDimensions/Genomics/Data_allocation_mastersheet_12-23-16.xlsx')
phyt.info <- read.xls('~/Dropbox/WarmAntDimensions/Phytotron 2013/Phytotron colonies 2013 Transcriptome.xlsx',3)
## phyt.info <- read.xls('~/Dropbox/WarmAntDimensions/Phytotron 2013/Aphaeno thermal data_field and phyto.xls')


