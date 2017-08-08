### Analytical script for the Apaenogaster genome
### MKLau 07February2017

if (substr(getwd(),(nchar(getwd()) - 2),nchar(getwd())) == "src"){setwd("..")}

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

geo.ctr <- split(ant.info[,c("Lon","Lat")],ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID.)
geo.ctr <- lapply(geo.ctr,function(x) apply(x,2,mean))
geo.ctr <- do.call(rbind,geo.ctr)


## gaemr info
gaemr.tab <- read.csv('data/storage/apg/gaemr-table.csv')
gaemr.tab <- gaemr.tab[!grepl("Name",gaemr.tab[,"Metric"]),]
gaemr.tab <- gaemr.tab[!grepl("Assembler",gaemr.tab[,"Metric"]),]
gaemr.tab <- split(gaemr.tab[,c("Metric","Value")],gaemr.tab[,"ID"])
gaemr.tab <- lapply(gaemr.tab,t)
metrics <- gaemr.tab[[1]]["Metric",]
gaemr.tab <- do.call(rbind,lapply(gaemr.tab,function(x) as.numeric(x["Value",])))
colnames(gaemr.tab) <- metrics


### Compare to other ant genomes and bees
## if (!ncbi.ant %in% ls()){source('src/apg_genome_compare.R')}
ref.sizes <- c(as.numeric(ncbi.ant[,'Size (Mb)']),ant.gen.size[,'1C Genome Size (Mb)'])

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

if (!'ant.gen.size' %in% ls()){source('src/ant_genome_size.R')}




