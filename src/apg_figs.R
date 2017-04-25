### Load other ant genome information
source('genome_size.R')

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
colnames(stats) <- sub(' ','',colnames(stats))
colnames(stats) <- sub(' ','',colnames(stats))
stats.ord <- c("Contigs","MaxContig","MeanContig","ContigN50","ContigN90","TotalContigLength","TotalGapLength","CapturedGaps","MaxGap","MeanGap","GapN50","Scaffolds","MaxScaffold","MeanScaffold","ScaffoldN50","ScaffoldN90","AssemblyGC","TotalScaffoldLength")

all(colnames(stats)[match(stats.ord,colnames(stats))] == stats.ord)
stats <- stats[,match(stats.ord,colnames(stats))]

### Make the table
print(xtable::xtable(stats),type = "latex",file = "../docs/manuscript/assembly_stats.tex")

### Inter-species comparisons
size.xlim <- range(as.numeric(c((stats[,'TotalScaffoldLength'] / (1000000)),ncbi.ant[,'Size (Mb)'])))
size.xlim <- c(floor(size.xlim[1]),ceiling(size.xlim[2]))
gc.xlim <- range(as.numeric(c(stats[,'AssemblyGC'],ncbi.ant[,'GC%'])))
gc.xlim <- c(floor(gc.xlim[1]),ceiling(gc.xlim[2]))

par(mfrow = c(1,2))
plot(density(as.numeric(ncbi.ant[,'Size (Mb)'])),main = '',xlab = 'Size (Mb)',xlim = size.xlim)
abline(v = stats[,'TotalScaffoldLength'] / (1000000),col = 'grey')
plot(density(as.numeric(ncbi.ant[,'GC%'])),main = '',xlab = 'GC%',xlim = gc.xlim)
abline(v = stats[,'AssemblyGC'],col = 'grey')


