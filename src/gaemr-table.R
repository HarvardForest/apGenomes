### Table of information from gaemr results
### MKLau 20February2017

library(XML)
source('helpers.R')
library(gdata)

## organize the results from broad
broad.info <- read.xls('data/storage/ap_genomes/SSF-1728_KeyforCollaborator.xlsx',2)
gaemr.tab <- lapply(as.character(broad.info$Sample.ID),get.broad,'data/storage/ap_genomes/')
out <- do.call(rbind,gaemr.tab)
write.csv('data/gaemr-table.csv')

## 


