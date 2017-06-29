### Table of information from gaemr results
### MKLau 20February2017

library(XML)
source('src/helpers.R')

## organize the results from broad
broad.info <- read.csv('data/storage/ap_genomes/broad_sample_key.csv')
gaemr.tab <- lapply(as.character(broad.info$Sample.ID),get.broad,'data/storage/ap_genomes/')
out <- do.call(rbind,gaemr.tab)
write.csv('data/gaemr-table.csv')

## 


