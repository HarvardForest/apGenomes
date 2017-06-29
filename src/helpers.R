library(magrittr)

#' get.broad Produces a table of output from Broad Institute GAEMR results
#' @return Table of values for metrics produced by GAEMR

get.broad <- function(x = 'sample id (Broad Inst)',path = '/Volumes/ellison_lab/ap_genomes/'){
    print(x)
    ## Grab info from PILON
    print('assembly')
    assembly.tab <- readHTMLTable(
        paste0(path,x,'/gaemr/html/assembly_stats.html')
               )
    print('pilon')
    if (any(grepl('PILON',dir(paste0(path,x,'/gaemr/html/'))))){
        pilon.tab <- readHTMLTable(
            paste0(path,x,'/gaemr/html/PILON.html')
            ) 
    }else{print('PILON.html missing.')}
    print('allpaths')
    allpaths.tab <- readHTMLTable(
        paste0(path,x,'/gaemr/html/GAEMR-ALLPATHS.html')
        )
    print('blasts')
    blast.loc <- dir(paste0(path,x,'/gaemr/table'), full = T)[grep(
        'blast_hit_taxonomy.table.txt',
        dir(paste0(path,x,'/gaemr/table'))
        )]
    blast.lines <- readLines(
        blast.loc
        )
    blast.tab <- (
        blast.lines[-1] %>% 
            strsplit(split = '\\|') %>%
                do.call(what = rbind)
        )
    colnames(blast.tab) <- blast.tab[1,]
    blast.tab <- data.frame(blast.tab[-1,])
    ant.hits <- grepl('formicidae',
                      blast.tab$'X.BestScoringTaxonomy.',
                      ign = T)
    pct.ants <- c('Percent.Ants',(sum(ant.hits) / length(ant.hits)) * 100)
    print('exporting')
    out <- list(
        assembly.tab$'Basic Assembly Stats',
        pilon.tab$'Pilon Metrics',
        allpaths.tab$'Kmer Stats Table',
        data.frame(t(pct.ants))
        )
    for (i in 1:length(out)){
        colnames(out[[i]]) <- c('Metric','Value')
    }
    out <- do.call(rbind,out)
    ID <- rep(x,nrow(out))
    out <- data.frame(ID,out)
}
