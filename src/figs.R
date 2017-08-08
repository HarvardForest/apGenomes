source("src/ApG.R")
source("src/ant_genome_size.R")


### Ant genomes previously sequenced
png("results/gaga_world.png",width = 700, height = 700)
ggplot(gaga.loc, aes(Country)) + 
    geom_bar(stat = "count") + 
        xlab("") + ylab("Frequency") +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=20,face="bold"),
                  axis.text.x = element_text(angle = 25, hjust = 1))
dev.off()

png("results/gaga_usa.png",width = 700, height = 700)
ggplot(gaga.loc[gaga.loc[,"Country"] == "USA",], aes(State)) + 
    geom_bar(stat = "count") + 
        xlab("") + ylab("Frequency") +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=20,face="bold"),
                  axis.text.x = element_text(angle = 25, hjust = 1))
dev.off()

### Ant Genome Comparisons

## GC Content
png("results/GC.png",width = 700, height = 700)
ggplot(data.frame(GC = as.numeric(ncbi.ant[,"GC%"])), aes(GC)) + 
    geom_histogram(binwidth = 2) + 
        geom_vline(xintercept = gaemr.tab[,"AssemblyGC"]) + 
            xlab("GC %") + ylab("Frequency") +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=20,face="bold"))
dev.off()

## Scaffold N50
png("results/ScaffoldN50.png",width = 700, height = 700)
ggplot(gaga, aes(Scaffold.N50.length.kb.)) + 
    geom_histogram(binwidth = 1200) + 
        geom_vline(xintercept = gaemr.tab[,"ScaffoldN50"] / 100) + 
            xlab("Scaffold N50") + ylab("Frequency") +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=20,face="bold"))
dev.off()

## Contig N50
png("results/ContigN50.png",width = 700, height = 700)
ggplot(gaga, aes(Contig.N50.length.kb.)) + 
    geom_histogram(binwidth = 10) + 
        geom_vline(xintercept = gaemr.tab[,"ContigN50"] / 1000) + 
            xlab("Contig N50") + ylab("Frequency") +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=20,face="bold"))
dev.off()

## Size
png("results/Assembly_gaga.png",width = 700, height = 700)
ggplot(gaga, aes(Assembly.size.Mb.)) + 
    geom_histogram(binwidth = 50) + 
        geom_vline(xintercept = gaemr.tab[,"TotalScaffoldLength"] / 1000000) + 
            xlab("Assembly Size") + ylab("Frequency") +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=20,face="bold"))
dev.off()

png("results/Assembly.png",width = 700, height = 700)
ggplot(data.frame(size = ref.sizes), aes(size)) + 
    geom_histogram(binwidth = 50) + 
        geom_vline(xintercept = gaemr.tab[,"TotalScaffoldLength"] / 1000000) + 
            xlab("Assembly Size") + ylab("Frequency") +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=20,face="bold"))
dev.off()

### Hits to ants and Aphaenogaster
png("results/gaemr_pc_ant.png",width = 700, height = 700)
ggplot(data.frame(gaemr.tab,names = rownames(mash)), aes(names,Percent.Ants)) + 
    geom_bar(stat = "identity") + 
        xlab("") + ylab("Percent Ant Hits (NCBI BLAST)") +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=20,face="bold"),
                  axis.text.x = element_text(angle = 25, hjust = 1))
dev.off()

png("results/gaemr_pc_apg.png",width = 700, height = 700)
ggplot(data.frame(gaemr.tab,names = rownames(mash)), 
       aes(names,Percent.Aphaenogaster)) + 
    geom_bar(stat = "identity") + 
        xlab("") + ylab("Percent Aph. Hits (NCBI BLAST)") +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=20,face="bold"),
                  axis.text.x = element_text(angle = 25, hjust = 1))
dev.off()




### Distance analyses
png("results/geoVmash.png")
df <- data.frame(geo = apg.gcd[lower.tri(apg.gcd)] / 1000, 
                 mash = mash[lower.tri(mash)])
fit <- lm(mash~geo,data = df)
summary(lm(mash~geo,data = df))
coef <- coef(fit)
ggplot(df,aes(geo,mash)) + geom_point() + 
    geom_abline(intercept = coef[1],slope = coef[2]) +
        xlab("Geographic Distance (km)") + ylab("Genomic Distance (MASH)") +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=20,face="bold"),
                  axis.text.x = element_text(angle = 0, hjust = 1))
dev.off()

### Latitude analysis
png("results/latVmash.png")
df <- data.frame(geo = lat.d[lower.tri(lat.d)] / 1000, 
                 mash = mash[lower.tri(mash)])
fit <- lm(mash~geo,data = df)
summary(lm(mash~geo,data = df))
coef <- coef(fit)
ggplot(df,aes(geo,mash)) + geom_point() + 
    geom_abline(intercept = coef[1],slope = coef[2]) +
        xlab("Latitudinal Distance (km)") + ylab("Genomic Distance (MASH)") +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=20,face="bold"),
                  axis.text.x = element_text(angle = 0, hjust = 1))
dev.off()

### By species
png("results/geoVmashXspp.png")

i <- 1
gcd.spp <- apg.gcd[i,-i]
mash.spp <- mash[i,-i]
df <- data.frame(geo = gcd.spp, 
                 mash = mash.spp,
                 rep(rownames(apg.gcd)[i],length(gcd.spp)))
for (i in 2:nrow(apg.gcd)){
gcd.spp <- apg.gcd[i,-i]
mash.spp <- mash[i,-i]
df <- rbind(df,data.frame(geo = gcd.spp, 
                 mash = mash.spp,
                 rep(rownames(apg.gcd)[i],length(gcd.spp))))
}
colnames(df) <- c("geo","mash","Species")
df[,"geo"] <- df[,"geo"] / 1000
ggplot(df) +
    geom_jitter(aes(geo,mash, colour = Species),) + geom_smooth(aes(geo,mash, colour = Species), method = lm, se = FALSE) +
        facet_wrap(~Species, scales="free_x") + scale_colour_discrete(guide = FALSE) + 
              labs(x = "Latitudinal Distance", y = "Genomic Distance (MASH)") +
            theme(axis.text=element_text(size=10),
                  axis.title=element_text(size=20,face="bold"),
                  axis.text.x = element_text(angle = 0))

dev.off()

system("cp results/*.png docs/esa2017")
