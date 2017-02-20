### Table of information from gaemr results
### MKLau 20February2017

library(XML)
library(magrittr)
source('helpers.R')

## organize the results from broad
gaemr.tab <- lapply(as.character(broad.info$Sample.ID),get.broad)
out <- do.call(rbind,gaemr.tab)
write.csv('./gaemr-table.csv')

## 


