library(magrittr)

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

