### Load other ant genome information

### Collect gaemr data for making tables and figures
apgs <- dir('~/storage/apg',full = TRUE)
apgs <- apgs[grepl('SM-',apgs)]
### overview
df <- list()
for (i in 1:length(apgs)){
    x <- apgs[i]
    tab <- lapply(paste0(x,'/gaemr/table/',
                         c("assembly.basic_assembly_stats.table.txt"
#,
#                           "assembly.blast_hit_taxonomy.table.txt",
#                           "assembly.contig_detail.table.txt",
#                           "assembly.gap_analysis.table.txt",
#                           "assembly.kmer_copy_number.table.txt",
#                           "assembly.scaffold_detail.table.txt"
                           )
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
stats.ord <- c("AssemblyGC","Contigs","MaxContig","MeanContig","ContigN50","ContigN90","TotalContigLength","TotalGapLength","CapturedGaps","MaxGap","MeanGap","GapN50","Scaffolds","MaxScaffold","MeanScaffold","ScaffoldN50","ScaffoldN90","TotalScaffoldLength")
all(colnames(stats)[match(stats.ord,colnames(stats))] == stats.ord)
stats <- stats[,match(stats.ord,colnames(stats))]
### Repeat for contaminants
df <- list()
for (i in 1:length(apgs)){
    tab <- sapply(c(
        "basic_assembly_stats.table.txt"
#,
#        "blast_hit_taxonomy.table.txt",
#        "contig_detail.table.txt",
#        "gap_analysis.table.txt",
#        "kmer_copy_number.table.txt",
#        "scaffold_detail.table.txt"
        ),function(x,y) y[grep(x,y)],
                  y = dir(paste0(apgs[i],'/unfiltered_gaemr/table/'),
                      full = TRUE))
    tab <- sapply(tab,readLines)
    tab <- tab[-1]
    tab <- do.call(rbind,strsplit(tab,'\\|'))
    df[[i]] <- tab
}
uf.stats <- df
uf.stats <- lapply(uf.stats,function(x) x[-1:-2,])
uf.stats <- lapply(uf.stats,sub, pattern = " ",replacement = "")
uf.stats <- lapply(uf.stats,sub, pattern = ",",replacement = "")
uf.stats <- lapply(uf.stats,sub, pattern = ",",replacement = "")
uf.stats.lab <- uf.stats[[1]][,1]
uf.stats <- lapply(uf.stats,function(x) x[,-1])
uf.stats <- do.call(cbind,uf.stats)
ap.lab <- substr(apgs,(nchar(apgs[1]) - 4),nchar(apgs))
colnames(uf.stats) <- ap.lab
uf.stats <- t(uf.stats)
uf.stats <- apply(uf.stats,2,as.numeric)
colnames(uf.stats) <- uf.stats.lab
rownames(uf.stats) <- ap.lab
colnames(uf.stats) <- sub(' ','',colnames(uf.stats))
colnames(uf.stats) <- sub(' ','',colnames(uf.stats))
uf.stats.ord <- c("Contigs","MaxContig","MeanContig","ContigN50","ContigN90","TotalContigLength","TotalGapLength","CapturedGaps","MaxGap","MeanGap","GapN50","Scaffolds","MaxScaffold","MeanScaffold","ScaffoldN50","ScaffoldN90","AssemblyGC","TotalScaffoldLength")
all(colnames(uf.stats)[match(uf.stats.ord,colnames(uf.stats))] == uf.stats.ord)
uf.stats <- uf.stats[,match(uf.stats.ord,colnames(uf.stats))]
### table for filtered

table.stats <- stats[,c(
    grep('GC',colnames(stats)),
    grep('gap',colnames(stats),ign = TRUE),
    grep('contig',colnames(stats),ign = TRUE),
    grep('Scaffold',colnames(stats))
    )]
tab.names <- c("AssemblyGC","TotalGapLength","CapturedGaps",
"Contigs","MaxContig","MeanContig","ContigN50","ContigN90","TotalContigLength",
"Scaffolds","MaxScaffold","MeanScaffold","ScaffoldN50","ScaffoldN90","TotalScaffoldLength")
table.stats <- table.stats[,match(tab.names,colnames(table.stats))]
colnames(table.stats) <- c("Assembly GC","Total Gap Length","Captured Gaps",
"Contigs","Max Contig","Mean Contig","Contig N50","Contig N90","Total Contig Length",
"Scaffolds","Max Scaffold","Mean Scaffold","Scaffold N50","Scaffold N90","Total Scaffold Length")
### table for UN-filtered
table.uf.stats <- uf.stats[,c(
    grep('GC',colnames(uf.stats)),
    grep('gap',colnames(uf.stats),ign = TRUE),
    grep('contig',colnames(uf.stats),ign = TRUE),
    grep('Scaffold',colnames(uf.stats))
    )]
tab.names <- c("AssemblyGC","TotalGapLength","CapturedGaps",
"Contigs","MaxContig","MeanContig","ContigN50","ContigN90","TotalContigLength",
"Scaffolds","MaxScaffold","MeanScaffold","ScaffoldN50","ScaffoldN90","TotalScaffoldLength")
table.uf.stats <- table.uf.stats[,match(tab.names,colnames(table.uf.stats))]
colnames(table.uf.stats) <- c("Assembly GC","Total Gap Length","Captured Gaps",
"Contigs","Max Contig","Mean Contig","Contig N50","Contig N90","Total Contig Length",
"Scaffolds","Max Scaffold","Mean Scaffold","Scaffold N50","Scaffold N90","Total Scaffold Length")
dif.stats <- uf.stats - stats
pc.contam <- (dif.stats[,'TotalScaffoldLength'] / uf.stats[,'TotalScaffoldLength']) * 100
pc.coverage <- (1 - (uf.stats[,'TotalGapLength'] / uf.stats[,'TotalScaffoldLength'])) * 100
### Add contaimant and species columns
table.stats <- data.frame('PercentRemoved' = pc.contam,stats)
### Convert to Mb
table.stats[,grepl('Total',colnames(table.stats))] <- table.stats[,grepl('Total',colnames(table.stats))] / 1000000
### ROund
table.stats <- round(table.stats,2)
table.stats <- data.frame('Species' = paste('A.',sample.info[match(broad.info[,'Collaborator.Sample.ID'],sample.info[,'broadID']),'spec_epithet']),table.stats)
### Remove N90
table.stats <- table.stats[,!grepl('N90',colnames(table.stats))]
table.stats <- table.stats[,!grepl('MaxScaffold',colnames(table.stats))]
table.stats <- table.stats[,!grepl('MaxContig',colnames(table.stats))]
table.stats <- table.stats[,!grepl('MeanGap',colnames(table.stats))]
table.stats <- table.stats[,!grepl('GapN50',colnames(table.stats))]
table.stats <- table.stats[,!grepl('Mean',colnames(table.stats))]
## Add percent coverage
table.stats <- data.frame(Species = table.stats[,1],
                          'TotalScaffoldLength' = table.stats[,ncol(table.stats)],
                          'PercentCoverage' = pc.coverage,
                          table.stats[,((ncol(table.stats)-1):2)])
