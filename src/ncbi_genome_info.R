### Gather genome size information from NCBI
getGeo <- function(x){
    if (any(grepl('latitude',tolower(x)))){
        x <- grep('latitude',tolower(x),val = T)
        x <- strsplit(x,'latitude')[[1]][2]
        x <- strsplit(x,'description')[[1]][1]
        x <- strsplit(x,'">')[[1]][2]
        x <- strsplit(x,'</a>')[[1]][1]
        x <- strsplit(x,' ')[[1]]
        x <- c(paste(x[2:1],collapse = ''),
               paste(x[4:3],collapse = ''))
        x <- sapply(x,gsub,pattern = '[s|w]',replacement = '-')
        x <- sapply(x,gsub,pattern = '[n|e]',replacement = '')
        x <- as.numeric(x)
        x
    }else{NA}
}

openBioSam <- function(x){
    system(paste0('open ','https://www.ncbi.nlm.nih.gov/biosample/',x))
}


### entrez_db_searchable('genome')
ants <- entrez_search(db = 'genome',term = 'Formicidae[ALL]')
ant.sum <- entrez_summary('genome',id = ants$ids)

### extract_from_esummary(ant.sum, "projectid") 

ncbi.euk <- readLines('ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt')
ncbi.euk <- do.call(rbind,strsplit(ncbi.euk,split = '\t'))
colnames(ncbi.euk) <- ncbi.euk[1,]
ncbi.euk <- ncbi.euk[-1,]
ncbi.insect <- ncbi.euk[ncbi.euk[,'SubGroup'] == 'Insects',]
ncbi.insect <- ncbi.insect[ncbi.insect[,'Genes'] != '-',]
ncbi.insect <- ncbi.insect[log10(as.numeric(ncbi.insect[,'Genes'])) > 2,]
ncbi.ant <- ncbi.euk[ncbi.euk[,'Assembly Accession'] %in% extract_from_esummary(ant.sum, "assembly_accession"),]

biosam.ant <- lapply(ncbi.ant[,'BioSample Accession'],function(accn) 
    readLines(paste0('https://www.ncbi.nlm.nih.gov/biosample/',accn))
)
names(biosam.ant) <- ncbi.ant[,1]

### lapply(ncbi.ant[,'BioSample Accession'],openBioSam)

ncbi.geo <- lapply(biosam.ant,getGeo)
ncbi.geo <- ncbi.geo[!unlist(lapply(ncbi.geo,function(x)any(is.na(x))))]
ncbi.geo <- ncbi.geo[!unlist(lapply(ncbi.geo,function(x) any(grepl('NA',x))))]
ncbi.geo <- do.call(rbind,ncbi.geo)
ncbi.geo <- ncbi.geo[,2:1]
colnames(ncbi.geo) <- c('Longitude','Latitude')

### write to disk
write.csv(ncbi.ant, file = "data/storage/apg/ncbi_ant.csv")
write.csv(ncbi.insect, file = "data/storage/apg/ncbi_insect.csv")
write.csv(ncbi.euk, file = "data/storage/apg/ncbi_euk.csv")
write.csv(ncbi.geo, file = "data/storage/apg/ncbi_geo.csv")
