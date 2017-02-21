### Analytical script for the Apaenogaster genome
### MKLau 07February2017

library(gdata)
library(prism)
library(ggplot2)
library(raster)

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
phyto.info <- read.xls('/Users/hermes/Dropbox/WarmAntDimensions/Phytotron\ 2013/Phytotron\ colonies\ 2013\ Transcriptome.xlsx',1)

## gaemr info
gaemr.tab <- read.csv('../docs/abstract_esa2017/tables/gaemr-table.csv')
metrics <- c("Mean Total Aligned Depth","Total Contig Length","Genome size estimate","Contig N50","Scaffold N50","Assembly GC","SNP rate                   ~","Snps")
gaemr.tab <- gaemr.tab[gaemr.tab$Metric %in% metrics,]
gaemr.tab$Value <- sub(' bases','',gaemr.tab$Value)
for (i in 1:nrow(gaemr.tab)){
    if (grepl('SNP rate',gaemr.tab$Metric[i])){
        x <- as.numeric(strsplit(gaemr.tab$Value[i],'\\/')[[1]])
        gaemr.tab$Value[i] <- round((x[1] / x[2]), 7)
    }
}
gaemr.tab$Value <- as.numeric(gsub(',','',as.character(gaemr.tab$Value)))
gaemr <- split(gaemr.tab,gaemr.tab$ID)
gaemr
for (i in 1:length(gaemr)){
    tmp <- as.character(gaemr[[i]]$Metric)
    gaemr[[i]] <- t(gaemr[[i]][,4])
    colnames(gaemr[[i]]) <- tmp
}
tmp <- do.call(rbind,gaemr)
gaemr <- data.frame(id = names(gaemr),tmp)
sample.info <- data.frame(sample.info,geo.ctr[match(sample.info$spec_epithet,rownames(geo.ctr)),])
sample.info
sample.info <- data.frame(seqID = broad.info$Sample.ID[match(sample.info$broadID,broad.info$Collaborator.Sample.ID)] , sample.info)
gaemr <- gaemr[match(sample.info$seqID,gaemr$id),]
apg.dat <- data.frame(sample.info,gaemr)

## Average depth of coverage == "Mean Total Aligned Depth"
## Assembly size == "Total Contig Length"
## Total genome size ==  "Genome size estimate" "Genome size estimate CN = 1" "Genome size estimate CN > 1"
## Contig N50 == "Contig N50"
## Scaffold N50 == "Scaffold N50"
## % genomic G+C base composition == "Assembly GC"
## SNP == "SNP rate                   ~"      "Snps"

## Haploid chromosome number == ?
## Protein-coding genes == ?
## Manually curated genes == ?
## Genes with EST support == ?
## miRNA == ?
## Total repeat content == ?


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
