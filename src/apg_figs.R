### Collect gaemr data for making tables and figures

apgs <- dir('~/storage/ap_genomes',full = TRUE)
apgs <- apgs[grepl('SM-',apgs)]

### overview
x <- apgs[2]



tab <- lapply(paste0(x,'/gaemr/table/',c("assembly.basic_assembly_stats.table.txt",
                "assembly.blast_hit_taxonomy.table.txt",
                "assembly.contig_detail.table.txt",
                "assembly.gap_analysis.table.txt",
                "assembly.gap_ss_analysis.table.txt",
                "assembly.kmer_copy_number.table.txt",
                "assembly.scaffold_detail.table.txt")
              ),
              readLines)
names(tab) <- unlist(lapply(tab,function(x) x[1]))
tab <- lapply(tab,function(x) x[-1])

df <- lapply(tab,function(x) 
    do.call(rbind,strsplit(x,split = "|"))
              )
