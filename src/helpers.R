#' get.broad Produces a table of output from Broad Institute GAEMR results
#' @return Table of values for metrics produced by GAEMR

get.broad <- function(x = 'sample id (Broad Inst)',path = '/Volumes/ellison_lab/ap_genomes/'){
    print(x)
    print('assembly info')
    assembly.loc <- dir(paste0(path,x,'/gaemr/table'), full = T)[grep(
        'assembly.basic_assembly_stats.table.txt',
        dir(paste0(path,x,'/gaemr/table'))
        )]
    assembly.lines <- readLines(assembly.loc)
    assembly.lines <- gsub('\\,','',assembly.lines)
    assembly.lines <- gsub('\\ ','',assembly.lines)
    assembly.tab <- do.call(rbind,strsplit(assembly.lines[-1],split = "\\|"))
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
                      blast.tab$'X.TaxonomicString',
                      ign = T)
    pct.ants <- c('Percent.Ants',(sum(ant.hits) / length(ant.hits)) * 100)
    aph.hits <- grepl('Aphaenogaster',
                      blast.tab$'X.TaxonomicString',
                      ign = T)
    pct.aph <- c('Percent.Aphaenogaster',(sum(aph.hits) / length(aph.hits)) * 100)
    print('exporting')
    out <- list(
        assembly.tab,
        data.frame(t(pct.ants)),
        data.frame(t(pct.aph))
        )
    for (i in 1:length(out)){
        colnames(out[[i]]) <- c('Metric','Value')
    }
    out <- do.call(rbind,out)
    ID <- rep(x,nrow(out))
    out <- data.frame(ID,out)
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)

    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    if (is.null(layout)) {
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        for (i in 1:numPlots) {
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                  layout.pos.col = matchidx$col))
        }
    }
}

as.mashdist <- function(x){
    lab <- unique(as.character(unlist(x[,1:2])))
    mat <- array(NA,dim = rep(length(lab),2))
    for (i in 1:nrow(x)){
        mat[lab == x[i,1],lab == x[i,2]] <- x[i,3]
    }
    rownames(mat) <- colnames(mat) <- lab
    mat
}

### Functions for accessing information from NCBI
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

### https://cran.r-project.org/web/packages/xtable/vignettes/xtableGallery.pdf
italic <- function(x){paste0('{\\emph{',x,'}}')}

### Example.
## mtcars$cyl <- factor(mtcars$cyl, levels = c("four","six","eight"),
##                      labels = c("four",italic("six"),"eight"))
## tbl <- ftable(mtcars$cyl, mtcars$vs, mtcars$am, mtcars$gear,
##               row.vars = c(2, 4),
##               dnn = c("Cylinders", "V/S", "Transmission", "Gears"))
## xftbl <- xtableFtable(tbl, method = "row.compact")
## print.xtableFtable(xftbl,
##                    sanitize.rownames.function = large,
##                    sanitize.colnames.function = bold,
##                    rotate.colnames = TRUE,
##                    rotate.rownames = TRUE)
