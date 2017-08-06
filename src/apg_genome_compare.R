### NCBI Genome Info
source('src/ncbi_genome_info.R')
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

