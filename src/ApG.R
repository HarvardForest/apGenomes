### Analytical script for the Apaenogaster genome
### MKLau 07February2017

### Analysis Outline

## Basic analysis of genomic similarity across species
## Basic analysis of target genes

source("https://bioconductor.org/biocLite.R")
library(magrittr)

## sample information

info <- gdata::read.xls('/Volumes/ellison_lab/ap_genomes/SSF-1728_KeyforCollaborator.xlsx',sheet = 2)

fasta <- (lapply(dir('/Volumes/ellison_lab/ap_genomes/heads/',
                    'fasta',full = TRUE),readLines) %>% #read fasta files from Odyssey
              lapply(FUN = paste0,collapse = '') %>% #paste lines together
                  lapply(FUN = strsplit,split = '>') %>% #split based on contig signifier
                      lapply(FUN = function(x) x[[1]][c(-1,-length(x))])
          )

## 0. genome quality (from the gaemer results)

## GC content
## percent contaminants
## 

## 1. Heat stress response genes @Stanton-Geddes2016


## 2. Cold tolerance genes @Stanton-Geddes2016


## 3. Map transcripts onto genome and look for inter-specific variance @Stanton-Geddes2016


## 4. Immunity @Roux2014


## 5. Olfaction @Roux2014


## 6. Mitochondrial function @Roux2014


## 7. Social chromosome? @Brelsford2014
