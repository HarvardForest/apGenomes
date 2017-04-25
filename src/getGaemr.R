### Compile gaemr output

apg.dir <- tail(commandArgs(trailingOnly = TRUE),1)

apg.dir <- '/Volumes/ellison_lab/ap_genomes'

apg.files <- dir(apg.dir)
apg.samples <- apg.files[grepl("SM",apg.files)]
apg.files <- paste0(apg.dir,"/",apg.samples,"/gaemr")
apg.stats <- lapply(paste0(apg.files,"/","table/",
                           "assembly.basic_assembly_stats.table.txt"),
                    readLines)
### gep CpG content ????
apg.stats <- lapply(apg.stats,function(x) x[-1:-3])
apg.stats <- lapply(apg.stats,strsplit,split="\\|")
apg.stats <- lapply(apg.stats,function(x) do.call(rbind,x))
apg.stats <- lapply(apg.stats,t)
apg.stats[2:length(apg.stats)] <- lapply(apg.stats[2:length(apg.stats)],
                                         function(x) x[-1,])
apg.stats <- do.call(rbind,apg.stats)
colnames(apg.stats) <- sub(" ","\\.",apg.stats[1,])
colnames(apg.stats) <- sub(" ","\\.",colnames(apg.stats))
apg.stats <- sub(",","",apg.stats)
apg.stats <- sub(",","",apg.stats)
apg.stats <- sub(" ","",apg.stats)
apg.stats <- apg.stats[-1,]
apg.stats <- data.frame(apg.stats)
apg.stats <- apply(apg.stats,2,as.numeric)
rownames(apg.stats) <- apg.samples
write.csv(apg.stats,"../docs/filtered_gaemr.csv")

