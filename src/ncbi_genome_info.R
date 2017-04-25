### Gather genome size information from NCBI
library(rentrez)
### entrez_db_searchable('genome')
ants <- entrez_search(db = 'genome',term = 'Formicidae[ALL]')
ant.sum <- entrez_summary('genome',id = ants$ids)
extract_from_esummary(ant.sum, "projectid") 
ncbi.euk <- readLines('ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt')
ncbi.euk <- do.call(rbind,strsplit(ncbi.euk,split = '\t'))
colnames(ncbi.euk) <- ncbi.euk[1,]
ncbi.euk <- ncbi.euk[-1,]
ncbi.ant <- ncbi.euk[ncbi.euk[,'Assembly Accession'] %in% extract_from_esummary(ant.sum, "assembly_accession"),]
