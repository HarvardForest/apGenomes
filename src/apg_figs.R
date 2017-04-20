### Collect gaemr data for making tables and figures

apgs <- dir('~/storage/ap_genomes',full = TRUE)
apgs <- apgs[grepl('SM-',apgs)]

### overview
df <- list()
for (i in 1:length(apgs)){
    x <- apgs[i]
    tab <- lapply(paste0(x,'/gaemr/table/',
                         c("assembly.basic_assembly_stats.table.txt",
                           "assembly.blast_hit_taxonomy.table.txt",
                           "assembly.contig_detail.table.txt",
                           "assembly.gap_analysis.table.txt",
                           "assembly.kmer_copy_number.table.txt",
                           "assembly.scaffold_detail.table.txt")
                         ),
                  readLines)
    names(tab) <- unlist(lapply(tab,function(x) x[1]))
    tab <- lapply(tab,function(x) x[-1])
    df[[i]] <- lapply(tab,function(x) 
        do.call(rbind,strsplit(x,split = "\\|"))
                      )
}


stats <- lapply(df,function(x) x[[1]])
stats <- lapply(stats,function(x) x[-1:-2,])
stats <- lapply(stats,sub, pattern = " ",replacement = "")
stats <- lapply(stats,sub, pattern = ",",replacement = "")
stats <- lapply(stats,sub, pattern = ",",replacement = "")
stats.lab <- stats[[1]][,1]
stats <- lapply(stats,function(x) x[,-1])
stats <- do.call(cbind,stats)
ap.lab <- substr(apgs,(nchar(apgs[1]) - 4),nchar(apgs))
colnames(stats) <- ap.lab
stats <- t(stats)
stats <- apply(stats,2,as.numeric)
colnames(stats) <- stats.lab
rownames(stats) <- ap.lab

pdf("apg_stats.pdf")
par(mfrow = c(3,6))
apply(stats,2,barplot,las = 2)
dev.off()
system("scp apg_stats.pdf matthewklau@fas.harvard.edu:public_html/")

library(xtable)
xtable(stats)
