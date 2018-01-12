if (substr(getwd(),(nchar(getwd()) - 2),nchar(getwd())) == "src"){setwd("..")}

### Check install of package dependencies
pkg <- c("gdata", "prism", "ggplot2", "raster", "AntWeb", "geosphere",
         "rnoaa", "gdata", "prism", "ggplot2", "raster", "vegan", "gdata",
         "tidyr", "stringr", "prism", "raster", "magrittr", "XML", "RCurl",
         "rlist", "rentrez","xtable","broom","ecodist","tibble","igraph","Rgraphviz")
## 
if (any(!pkg %in% installed.packages()[,1])){
    unlist(sapply(pkg,function(pkg) 
        if(!pkg %in% installed.packages()[,1]){
            install.packages(pkg)
        }
                  ))
}

### Load packages and user defined functions
all(unlist(lapply(pkg, require, character.only = TRUE)))

### user defined functions
reorder_size <- function(x,decreasing = TRUE) {
    factor(x, levels = names(sort(table(x), decreasing = decreasing)))
}

get.spp <- function(x){
    x <- grep("species", strsplit(x,split = ";")[[1]], value = T)
    return(gsub("species=","",x))
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

get.mash.p <- function(x){
    lab <- unique(as.character(unlist(x[,1:2])))
    mat <- array(NA,dim = rep(length(lab),2))
    for (i in 1:nrow(x)){
        mat[lab == x[i,1],lab == x[i,2]] <- x[i,4]
    }
    rownames(mat) <- colnames(mat) <- lab
    mat
}


### Restrict strings to letters
onlyAz <- function(x){
    paste(strsplit(x, split = "")[[1]][
        strsplit(x, split = "")[[1]] %in% 
            c(LETTERS,letters)], collapse = "")
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


### rebuild the stats tables
make.stats.table <- FALSE
broad.info <- read.csv('data/storage/apg/broad_sample_key.csv')
sample.info <- read.csv('data/storage/apg/colony_locations.csv')
ant.info <- read.csv('data/storage/apg/RADseq_mastersheet_2014.csv')
ant.info <- ant.info[ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID. %in% na.omit(sample.info$spec_epithet),]
ant.info <- ant.info[!ant.info$State == "",]
ant.geo <- read.csv('data/storage/apg/ant_sites.csv')
ant.geo[,c("Lon","Lat")] <- apply(ant.geo[,c("Lon","Lat")],2,as.numeric)
ant.info <- data.frame(ant.info,
                       ant.geo[match(as.character(ant.info$Locale),ant.geo$Site),c("Lon","Lat")])
ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID. <- 
    as.character(ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID.)
ant.info <- na.omit(ant.info)
ant.color <- ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID.
ant.factor <- factor(ant.color)

geo.ctr <- split(ant.info[,c("Lon","Lat")],ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID.)
geo.ctr <- lapply(geo.ctr,function(x) apply(x,2,mean))
geo.ctr <- do.call(rbind,geo.ctr)

## source("src/ant_genome_size.R")
### https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-8-64
### Genome sizes were analyzed by flow cytometry
theurl <- getURL("https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-8-64",.opts = list(ssl.verifypeer = FALSE) )
tables <- readHTMLTable(theurl)
tables <- list.clean(tables, fun = is.null, recursive = FALSE)
ant.gen.size <- tables[[1]]
ant.gen.size <- ant.gen.size[!apply(ant.gen.size,1,function(x) any(grepl('MEAN',x))),]

ant.gen.size <- ant.gen.size[ant.gen.size[,1] != '',]
ant.gen.size[,'1C Genome Size (Mb)'] <- as.numeric(as.character(ant.gen.size[,'1C Genome Size (Mb)']))

### Genome Size by Location
ant.gs.world <- table(substr(ant.gen.size[,"Collection Info"],1,3))
ant.gs.usa <- table(grep("USA",ant.gen.size[,"Collection Info"],value = TRUE))
names(ant.gs.usa) <- gsub("USA: ","",names(ant.gs.usa))
ags.w <- ant.gs.world[order(ant.gs.world,decreasing = TRUE)]
ags.w <- data.frame(region = names(ags.w),count = as.numeric(ags.w))
ags.u <- ant.gs.usa[order(ant.gs.usa,decreasing = TRUE)]
ags.u <- data.frame(region = names(ags.u),count = as.numeric(ags.u))
ant.count.spp <- strsplit(as.character(ant.gen.size[,"Species"],1,4), split = " ")
ant.count.spp <- lapply(ant.count.spp, function(x) paste(x[1], x[2]))
ant.spp <- unique(unlist(ant.count.spp))
ant.count.gen <- table(do.call(rbind, strsplit(unique(unlist(ant.count.spp)), split = " "))[,1])

### From the GAGA group
gaga <- na.omit(read.csv("data/gaga_genome_info.csv"))
gaga <- gaga[gaga[,"Contig.N50.length.kb."] != "No data ", ]
gaga[,"Scaffold.N50.length.kb."] <- as.numeric(gsub(",","",as.character(
    gaga[,"Scaffold.N50.length.kb."])))
gaga[,"Contig.N50.length.kb."] <- as.numeric(gsub(",","",as.character(
    gaga[,"Contig.N50.length.kb."])), warn = FALSE)
gaga.loc <- do.call(rbind,strsplit(as.character(gaga[,"Location"]),","))
gaga.loc <- apply(gaga.loc,2,function(x) sapply(x, onlyAz))
colnames(gaga.loc) <- c("City","State","Country")
## Change ordering
gaga.loc <- data.frame(gaga.loc)



## source("src/apg_mash.R")
### Analyze the mash distnaces

## gaemr info
gaemr.tab <- read.csv('data/storage/apg/gaemr-table.csv')
gaemr.tab <- gaemr.tab[!grepl("Name",gaemr.tab[,"Metric"]),]
gaemr.tab <- gaemr.tab[!grepl("Assembler",gaemr.tab[,"Metric"]),]
gaemr.tab <- split(gaemr.tab[,c("Metric","Value")],gaemr.tab[,"ID"])
gaemr.tab <- lapply(gaemr.tab,t)
metrics <- gaemr.tab[[1]]["Metric",]
gaemr.tab <- do.call(rbind,lapply(gaemr.tab,function(x) as.numeric(x["Value",])))
colnames(gaemr.tab) <- metrics

broad.info <- read.csv("data/storage/apg/broad_sample_key.csv")
broad.info[,"Collaborator.Sample.ID"] <- as.character(broad.info[,"Collaborator.Sample.ID"])
broad.info[broad.info[,"Collaborator.Sample.ID"] == "arudis1","Collaborator.Sample.ID"] <- "rud1"
broad.info[broad.info[,"Collaborator.Sample.ID"] == "rud6","Collaborator.Sample.ID"] <- "rud2"

## mash organziation
## MASH scripts are located in apGenomes/bin
geno.info <- read.csv("data/storage/apg/gen_seq_info.csv")
mash.txt <- read.table("data/storage/apg/mash_dist.txt",sep = "\t")
mash <- as.mashdist(mash.txt) 
rownames(mash) <- colnames(mash) <- paste0(as.character(geno.info[sapply(geno.info[,1],grep,x = rownames(mash)),2]),c("","","",1,2,rep("",nrow(geno.info) - 5)))
ncbi.gen <- mash
ncbi.rv <- c(11,19,14,15,9,5,7,10,6,4,8,12,1,2,24,3,22,23,20,25,26,18,21,16,13,17)
mash <- mash[grep("Aphaenogaster",rownames(mash)),grep("Aphaenogaster",rownames(mash))]

## mash network for ants
mashP.ncbi <- get.mash.p(mash.txt) 
mashD.ncbi <- as.mashdist(mash.txt) 
mash.net <- mashD.ncbi
mash.net[ mashD.ncbi > 0.05] <- 0
## gplot(mash.net)

## distances
## geographic information
geo.ctr <- split(ant.info[,c("Lon","Lat")],ant.info$Species..varying.ID.sources...Bernice.if.different.from.original.ID.)
geo.ctr <- lapply(geo.ctr,function(x) apply(x,2,mean))
geo.ctr <- do.call(rbind,geo.ctr)

ap.ctr <- do.call(rbind,list(rud1 = geo.ctr[6,],
                             rud6 = geo.ctr[6,],
                             pic1 = geo.ctr[5,],
                             mia1 = geo.ctr[4,],
                             ful1 = geo.ctr[3,],
                             ash1 = geo.ctr[1,],
                             flo1 = geo.ctr[2,]))

### apg sample geographic info
### arudis and picea are from google earth
### All other coords are from ant_sites.csv
apg.geo  <- do.call(rbind,list('arudis1' = c(-78.9830464,36.0200847),
                 'rud6' = c(-78.9830464,36.0200847),
                 'pic1' = c(-72.5847494,42.6004513),
                 'mia1' = c(-82.301773,29.657955),
                 'ful1' = c(-82.514575,32.692384),
                 'ash1' = c(-82.031176,29.785325),
                 'flo1' = c(-82.031176,29.785325)))
colnames(apg.geo) <- c('Longitude','Latitude')
apg.geo.labs <- paste0(substr(rownames(apg.geo),1,3),
                       substr(rownames(apg.geo),nchar(rownames(apg.geo)),
                              nchar(rownames(apg.geo))))
apg.geo.labs[1] <- 'rud1'

### geographic distance for samples
apg.gd <- array(NA,dim = rep(nrow(apg.geo),2))
rownames(apg.gd) <- colnames(apg.gd) <- rownames(apg.geo)
for (i in 1:nrow(apg.geo)){
    for (j in 1:nrow(apg.geo)){
        apg.gd[i,j] <- distm (apg.geo[i,], apg.geo[j,], 
                               fun = distHaversine)
    }
}

apg.gcd <- array(NA,dim = rep(nrow(ap.ctr),2))
rownames(apg.gcd) <- colnames(apg.gcd) <- rownames(ap.ctr)
for (i in 1:nrow(ap.ctr)){
    for (j in 1:nrow(ap.ctr)){
        apg.gcd[i,j] <- distm (ap.ctr[i,], ap.ctr[j,], 
                               fun = distHaversine)
    }
}

mash.d <- as.dist(mash)
geo.d <- as.dist(apg.gd)
geo.cd <- as.dist(apg.gcd)
geod.pic <- apg.gcd[rownames(apg.gcd) != "pic1","pic1"]
dist.pic <- apg.geo[,'Latitude'] - apg.geo[,'Latitude']['pic1']
mash.dpic <- mash[rownames(mash) != "Aphaenogaster picea","Aphaenogaster picea"]

mash <- mash[order(apg.gcd[,"pic1"],decreasing = F),
             order(apg.gcd[,"pic1"], decreasing = F)]
gcd.pic <- apg.gcd["pic1",rownames(apg.gcd) != "pic1"]
apg.gcd <- apg.gcd[order(apg.gcd[,"pic1"]),order(apg.gcd[,"pic1"])]

## diag(mash) <- NA
## diag(apg.gcd) <- NA

### Reorg geo.ctr for latitude
geo.lat <- geo.ctr[c(6,6,5,4,3,1,2),"Lat"]
names(geo.lat) <- c("rudis1","rudis6","picea","miamiana","fulva","ashmeadi","floridana")

### Size similarity
size.d <- as.matrix(dist(as.matrix(gaemr.tab[,"TotalScaffoldLength"])))
size.dpic <- size.d[rownames(mash) == "pic1",rownames(mash) != "pic1"]

### GC similiarty
gc.d <- as.matrix(dist(as.matrix(gaemr.tab[,"AssemblyGC"])))
gc.dpic <- size.d[rownames(mash) == "pic1",rownames(mash) != "pic1"]

### Lat distance
lat.d <- as.matrix(dist(as.matrix(geo.lat)))
lat.d.reorder <- match(c("picea","rudis1","rudis6","fulva","floridana","miamiana","ashmeadi"),rownames(lat.d))
lat.d <- lat.d[lat.d.reorder,lat.d.reorder]

## source("src/apg_prism.R")
### Get prism data for coordinates
### This is borrowed code:
### http://eremrah.com/articles/How-to-extract-data-from-PRISM-raster/

## SEE: https://github.com/ropensci/prism
## Get PRISM data
if (!dir.exists("~/prismtmpnormals")){
    options(prism.path = "~/prismtmpnormals")
    get_prism_normals(type="ppt", "800m", annual = TRUE)
    get_prism_normals(type="tmin", "800m", mon = 1, annual = TRUE)
    get_prism_normals(type="tmax", "800m", mon = 7, annual = TRUE)
}

options(prism.path = "~/prismtmpnormals")
## Stack files
mystack <- ls_prism_data() %>%  prism_stack()  
## Get proj from raster stack
mycrs <- mystack@crs@projargs
## My points
mypoints <- data.frame(id = rownames(apg.geo),
                       lat = apg.geo[,"Latitude"],
                       long = apg.geo[,"Longitude"]
)
## Convert points to spatial points data frame
coordinates(mypoints) <- c('long', 'lat')
proj4string(mypoints) <- CRS(mycrs)
## Extract data from raster
data <- data.frame(coordinates(mypoints), mypoints$id,
                   extract(mystack, mypoints))
## Rename column headers
colnames(data)[4:6] <- do.call(rbind,strsplit(colnames(data)[4:6],"_"))[,2]
rownames(data)[rownames(data) == "arudis1"] <- "rud1"
data[,"mypoints.id"] <- as.character(data[,"mypoints.id"])
data[rownames(data) == "rud1","mypoints.id"] <- "rud1"

## Rename object
clim.data <- data

## create a climate and distance objects
## Climate distances
clim.d <- dist(apply(clim.data[,c("ppt","tmax","tmin")],2,function(x) (x - mean(x)) / sd(x)))
## Temperature distances
temp.d <- dist(apply(clim.data[,c("tmax","tmin")],2,function(x) (x - mean(x)) / sd(x)))
### Precip distances
ppt.d <- as.matrix(dist(clim.data[,c("ppt")]))
rownames(ppt.d) <- colnames(ppt.d) <- rownames(as.matrix(clim.d))
ppt.d <- as.dist(ppt.d)
### Tmax distances
tmax.d <- as.matrix(dist(clim.data[,c("tmax")]))
rownames(tmax.d) <- colnames(tmax.d) <- rownames(as.matrix(clim.d))
tmax.d <- as.dist(tmax.d)
### Tmin distances
tmin.d <- as.matrix(dist(clim.data[,c("tmin")]))
rownames(tmin.d) <- colnames(tmin.d) <- rownames(as.matrix(clim.d))
tmin.d <- as.dist(tmin.d)

### Analysis Outline
## sample information

### Compare to other ant genomes and bees
## This is the data to make the genome to gene/contig fig
### NCBI Genome Info
if (!any(grepl("ncbi_ant.csv", dir("data/storage/apg")))){
    source('src/ncbi_genome_info.R')
}
ncbi.ant <- read.csv("data/storage/apg/ncbi_ant.csv")
colnames(ncbi.ant)[1] <- 'Organism'

### Parse the reference sizes
ref.sizes <- c(as.numeric(ncbi.ant[,'Size..Mb.']),ant.gen.size[,'1C Genome Size (Mb)'])

### AntWeb Info
## aw.apg <- list()
## for (i in 1:nrow(sample.info)){
##     aw.apg[[i]] <- aw_data(
##         scientific_name = 
##             paste('Aphaenogaster',
##                   sample.info[i,'spec_epithet']),
##         georeferenced = TRUE)
## }
## ### ant web info for ncbi
## aw.ncbi <- list()
## for (i in 1:nrow(sample.info)){
##     aw.ncbi[[i]] <- aw_data(
##         scientific_name = 
##             as.character(ncbi.ant[i,2]),
##         georeferenced = TRUE)
## }

### Inter-species comparisons
### Load other ant genome information
apgs <- dir('data/storage/apg',full = TRUE)
sample.info <- read.csv('data/storage/apg/colony_locations.csv')
broad.info <- read.csv('data/storage/apg/broad_sample_key.csv')

### Collect gaemr data for making tables and figures
if (make.stats.table){
    apgs <- dir('data/storage/apg',full = TRUE)
    apgs <- apgs[grepl('SM-',apgs)]
### overview
    df <- list()
    for (i in 1:length(apgs)){
        x <- apgs[i]
        tab <- lapply(paste0(x,'/gaemr/table/',
                             c("assembly.basic_assembly_stats.table.txt"
                                        ,
                                         "assembly.blast_hit_taxonomy.table.txt"
                                        #,
                                        # "assembly.contig_detail.table.txt",
                                        # "assembly.gap_analysis.table.txt",
                                        # "assembly.kmer_copy_number.table.txt",
                                        # "assembly.scaffold_detail.table.txt"
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
    apgs <- dir('data/storage/apg/20161122/',full = TRUE)
    apgs <- grep('SM-A',apgs,value = TRUE)
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
                   "Contigs","MaxContig","MeanContig","ContigN50","ContigN90",
                   "TotalContigLength",
                   "Scaffolds","MaxScaffold","MeanScaffold","ScaffoldN50",
                   "ScaffoldN90","TotalScaffoldLength")
    table.stats <- table.stats[,match(tab.names,colnames(table.stats))]
    colnames(table.stats) <- c("Assembly GC","Total Gap Length","Captured Gaps",
                               "Contigs","Max Contig","Mean Contig","Contig N50",
                               "Contig N90",
                               "Total Contig Length",
                               "Scaffolds","Max Scaffold","Mean Scaffold",
                               "Scaffold N50","Scaffold N90","Total Scaffold Length")
### table for UN-filtered
    table.uf.stats <- uf.stats[,c(
        grep('GC',colnames(uf.stats)),
        grep('gap',colnames(uf.stats),ign = TRUE),
        grep('contig',colnames(uf.stats),ign = TRUE),
        grep('Scaffold',colnames(uf.stats))
        )]
    tab.names <- c("AssemblyGC","TotalGapLength","CapturedGaps",
                   "Contigs","MaxContig","MeanContig","ContigN50","ContigN90",
                   "TotalContigLength",
                   "Scaffolds","MaxScaffold","MeanScaffold","ScaffoldN50",
                   "ScaffoldN90","TotalScaffoldLength")
    table.uf.stats <- table.uf.stats[,match(tab.names,colnames(table.uf.stats))]
    colnames(table.uf.stats) <- c("Assembly GC","Total Gap Length","Captured Gaps",
                                  "Contigs","Max Contig","Mean Contig","Contig N50",
                                  "Contig N90","Total Contig Length",
                                  "Scaffolds","Max Scaffold","Mean Scaffold",
                                  "Scaffold N50","Scaffold N90","Total Scaffold Length")
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
### Write to csv
    write.csv(table.stats,'data/apg_summary.csv')
### Write latex
    xtab <- xtable(table.stats, align = c("l",">{\\itshape}l",rep('l',(ncol(table.stats)-1))))
    names(xtab) <- c("Species" , "Percent Removed" , "GC Content" , "Contigs" , "ContigN50" , "Total Contig Length" , "Total Gap Length" , "Captured Gaps" , "Max Gap Length" , "Scaffolds" , "Scaffold N50" , "Total Scaffold Length")
    capture.output(xtab,file = "docs/manuscript/seq_info_tab.tex")
    size.xlim <- range(as.numeric(c((
        stats[,'TotalScaffoldLength'] / (1000000)),
                                    ncbi.ant[,'Size..Mb.)'])))
    size.xlim <- c(floor(size.xlim[1]),ceiling(size.xlim[2]))
    gc.xlim <- range(as.numeric(c(
        stats[,'AssemblyGC'],
        ncbi.ant[,'GC%'])))
    gc.xlim <- c(floor(gc.xlim[1]),ceiling(gc.xlim[2]))
    ref.sizes <- c(as.numeric(ncbi.ant[,'Size..Mb.)']),
                   ant.gen.size[,'1C Genome Size (Mb)'])
    ## NCBI ant information
    ncbi.xtab <- ncbi.ant[,c('X.Organism.Name','BioProject.Accession','BioSample.Accession')]
    colnames(ncbi.xtab) <- c('Ant Species','BioProject Accession','BioSample Accession')
    rownames(ncbi.xtab) <- ncbi.xtab[,"Ant Species"]
    ncbi.xtab <- ncbi.xtab[,-1]
    ncbi.xtab <- ncbi.xtab[order(rownames(ncbi.xtab)),]
    ncbi.xtab <- xtable::xtable(ncbi.xtab, caption = "NCBI genome database accession information for the previously sequenced ant genomes.")
    ## Table: create ncbi_ants 
    print(ncbi.xtab,
          type = "latex",
          file = "results/ncbi_ants.tex",
          sanitize.rownames.function = italic,
          include.rownames = TRUE,
          include.colnames = TRUE
          )
}

### Genomic biogeography = lat/lon, distance, climate = temperature, climatic similarities
### Tests of climate variable correlations
cor.tminmax <- cor.test(clim.data[,"tmax"],clim.data[,"tmin"])
cor.ppttmax <- cor.test(clim.data[,"ppt"],clim.data[,"tmax"])
cor.ppttmin <- cor.test(clim.data[,"ppt"],clim.data[,"tmin"])
clim.cor <- do.call(rbind,lapply(lst(cor.tminmax,cor.ppttmax,cor.ppttmin),tidy))
write.csv(clim.cor,"results/clim_cor.csv")

## Table: create climate correlation table
climcor.xtab <- xtable::xtable(clim.cor, digits = 5)
print(climcor.xtab,
      type = "latex",
      file = "results/clim_cor.tex",
      include.rownames = TRUE,
      include.colnames = TRUE
      )


## Mantel of MASH
## check ordering 
if (!all(
    colnames(as.matrix(mash.d)) == c("Aphaenogaster rudis1", 
                "Aphaenogaster rudis2", 
                "Aphaenogaster picea",
                "Aphaenogaster miamiana",
                "Aphaenogaster fulva",
                "Aphaenogaster ashmeadi", 
                "Aphaenogaster floridana")) & 
    all(colnames(as.matrix(temp.d)) == c(
                    "rud1", "rud6", "pic1", "mia1", "ful1", "ash1", "flo1")) & 
    (all(colnames(as.matrix(temp.d)) == colnames(as.matrix(clim.d))) & 
         all(colnames(as.matrix(clim.d)) == colnames(as.matrix(geo.cd))) & 
             all(colnames(as.matrix(clim.d)) == colnames(as.matrix(ppt.d))))){
    warning("Distance matrix ordering incorrect!")
}

set.seed(1649)
mantel.tab <- list(clim.geo = ecodist::mantel(clim.d~geo.d, nperm = 10000),
                  mash.clim = ecodist::mantel(mash.d~clim.d, nperm = 10000),
                  mash.temp = ecodist::mantel(mash.d~temp.d, nperm = 10000),
                  temp.geo = ecodist::mantel(temp.d~geo.d, nperm = 10000),
                  mash.temp_geo = ecodist::mantel(mash.d~temp.d + geo.d, nperm = 10000),
                  mash.tmin_geo = ecodist::mantel(mash.d~tmin.d+geo.d, nperm = 10000),
                  mash.tmax_geo = ecodist::mantel(mash.d~tmax.d+geo.d, nperm = 10000),
                  ppt = ecodist::mantel(mash.d~ppt.d, nperm = 10000),
                  ppt.geo = ecodist::mantel(ppt.d~geo.d, nperm = 10000),
                  mash.ppt_geo = ecodist::mantel(mash.d~ppt.d+geo.d, nperm = 10000)
)
mantel.tab <- do.call(rbind,mantel.tab)
mantel.tab <- mantel.tab[,c(1,2)]
mantel.tab
write.csv(mantel.tab,"results/mantel_tab.csv")
mantel.xtab <- xtable::xtable(mantel.tab, digits = 5)
## Table: create mantel
print(mantel.xtab,
      type = "latex",
      file = "results/mantel_tab.tex",
      include.rownames = TRUE,
      include.colnames = TRUE
      )

## Geography and climate
mantel.geog <- ecodist::mantel(clim.d~geo.cd, nperm = 10000)

### gaemr correlations
gc.dat <- data.frame(GC = gaemr.tab[,"AssemblyGC"], 
                     Latitude = apg.geo[,"Latitude"], 
                     Precipitation = clim.data[,"ppt"])
size.dat <- data.frame(GenomeSize = gaemr.tab[,"TotalScaffoldLength"], 
                     Latitude = apg.geo[,"Latitude"], 
                     Precipitation = clim.data[,"ppt"])

(summary(lm(GC~Latitude*Precipitation, data = gc.dat)))
(summary(lm(GenomeSize~Latitude*Precipitation, data = size.dat)))


### Figures
### Ant genomes previously sequenced
png("results/gaga_world.png",width = 700, height = 700)
ggplot(gaga.loc, aes(reorder_size(Country))) + 
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
ggplot(data.frame(GC = as.numeric(ncbi.ant[,"GC."])), aes(GC)) + 
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
png("results/geoVmash.png",width = 700, height = 700)
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
png("results/latVmash.png",width = 700, height = 700)
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

## Heatmap
png("results/ncbi_heat.png",width = 1200, height = 800, pointsize = 25)
heatmap(ncbi.gen, 
        RowSideColors=rainbow(nlevels(geno.info[,"subfamily"]))[as.numeric(geno.info[,"subfamily"])],
        symm = T, margins = c(1,10),labCol = "")
dev.off()
png("results/apg_heat.png",width = 1200, height = 800, pointsize = 25)
heatmap(mash, 
        symm = T, margins = c(1,10),labCol = "")
dev.off()

### BioGeographic plots
### By species

### Geographic distance
stats.mash.spp <- list()
png("results/geoVmashXspp.png")
i <- 1
gcd.spp <- apg.gcd[i,-i]
mash.spp <- mash[i,-i]
df <- data.frame(geo = gcd.spp, 
                 mash = mash.spp,
                 rep(rownames(apg.gcd)[i],length(gcd.spp)))
stats.mash.spp[[i]] <- cor.test(gcd.spp,mash.spp^2)
for (i in 2:nrow(apg.gcd)){
gcd.spp <- apg.gcd[i,-i]
mash.spp <- mash[i,-i]
stats.mash.spp[[i]] <- cor.test(gcd.spp,mash.spp^2)
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

### Climate distance
rm(df)
stats.mash.spp <- list()
png("results/climVmashXspp.png")
i <- 1
clim.spp <- as.matrix(clim.d)[i,-i]
mash.spp <- as.matrix(mash.d)[i,-i]
df <- data.frame(clim = clim.spp, 
                 mash = mash.spp,
                 rep(rownames(as.matrix(clim.d))[i],length(clim.spp)))
colnames(df) <- c("clim","mash","Species")
stats.mash.spp[[i]] <- cor.test(clim.spp,mash.spp^2)
for (i in 2:nrow(as.matrix(clim.d))){
    clim.spp <- as.matrix(clim.d)[i,-i]
    mash.spp <- as.matrix(mash.d)[i,-i]
    stats.mash.spp[[i]] <- cor.test(clim.spp,mash.spp^2)
    df <- rbind(df,data.frame(clim = clim.spp, 
                              mash = mash.spp,
                              Species = rep(rownames(as.matrix(clim.d))[i],length(clim.spp))))
}
ggplot(df) +
    geom_jitter(aes(clim,mash, colour = Species),) + geom_smooth(aes(clim,mash, colour = Species), method = lm, se = FALSE) +
        facet_wrap(~Species, scales="free_x") + scale_colour_discrete(guide = FALSE) + 
              labs(x = "Climate Distance", y = "Genomic Distance (MASH)") +
            theme(axis.text=element_text(size=10),
                  axis.title=element_text(size=20,face="bold"),
                  axis.text.x = element_text(angle = 0))
dev.off()

### Precip distance
rm(df)
stats.mash.spp <- list()
png("results/pptVmashXspp.png")
i <- 1
ppt.spp <- as.matrix(ppt.d)[i,-i]
mash.spp <- as.matrix(mash.d)[i,-i]
df <- data.frame(ppt = ppt.spp, 
                 mash = mash.spp,
                 rep(rownames(as.matrix(ppt.d))[i],length(ppt.spp)))
colnames(df) <- c("ppt","mash","Species")
stats.mash.spp[[i]] <- cor.test(ppt.spp,mash.spp^2)
for (i in 2:nrow(as.matrix(ppt.d))){
    ppt.spp <- as.matrix(ppt.d)[i,-i]
    mash.spp <- as.matrix(mash.d)[i,-i]
    stats.mash.spp[[i]] <- cor.test(ppt.spp,mash.spp^2)
    df <- rbind(df,data.frame(ppt = ppt.spp, 
                              mash = mash.spp,
                              Species = rep(rownames(as.matrix(ppt.d))[i],length(ppt.spp))))
}
ggplot(df) +
    geom_jitter(aes(ppt,mash, colour = Species),) + geom_smooth(aes(ppt,mash, colour = Species), method = lm, se = FALSE) +
        facet_wrap(~Species, scales="free_x") + scale_colour_discrete(guide = FALSE) + 
              labs(x = "Precipitation Difference", y = "Genomic Distance (MASH)") +
            theme(axis.text=element_text(size=10),
                  axis.title=element_text(size=20,face="bold"),
                  axis.text.x = element_text(angle = 0))
dev.off()

### Temp distance
rm(df)
stats.mash.spp <- list()
png("results/tempVmashXspp.png")
i <- 1
temp.spp <- as.matrix(temp.d)[i,-i]
mash.spp <- as.matrix(mash.d)[i,-i]
df <- data.frame(temp = temp.spp, 
                 mash = mash.spp,
                 rep(rownames(as.matrix(temp.d))[i],length(temp.spp)))
colnames(df) <- c("temp","mash","Species")
stats.mash.spp[[i]] <- cor.test(temp.spp,mash.spp^2)
for (i in 2:nrow(as.matrix(temp.d))){
    temp.spp <- as.matrix(temp.d)[i,-i]
    mash.spp <- as.matrix(mash.d)[i,-i]
    stats.mash.spp[[i]] <- cor.test(temp.spp,mash.spp^2)
    df <- rbind(df,data.frame(temp = temp.spp, 
                              mash = mash.spp,
                              Species = rep(rownames(as.matrix(temp.d))[i],length(temp.spp))))
}
ggplot(df) +
    geom_jitter(aes(temp,mash, colour = Species),) + geom_smooth(aes(temp,mash, colour = Species), method = lm, se = FALSE) +
        facet_wrap(~Species, scales="free_x") + scale_colour_discrete(guide = FALSE) + 
              labs(x = "Temperature Distance", y = "Genomic Distance (MASH)") +
            theme(axis.text=element_text(size=10),
                  axis.title=element_text(size=20,face="bold"),
                  axis.text.x = element_text(angle = 0))
dev.off()

### Ordination
nms <- nmds(mash.d)
ord <- nmds.min(nms)
vec <- envfit(ord, clim.data[,c("long","lat", "tmax", "tmin", "ppt")])
clim.data[clim.data[,"mypoints.id"] == "rud6" , "mypoints.id"] <- "rud2"

png("results/apg_ord.png")
plot(ord, xlab = "NMDS 1", ylab = "NMDS 2", pch = "")
text(ord, labels = clim.data[,"mypoints.id"], col = "black")
plot(vec, col = "darkgrey")
dev.off()

### Parsing the impact of geography + climate on MASH
### Simple correlation vs mantel
summary(lm(ppt~tmax,data = clim.data))
summary(lm(ppt~tmin,data = clim.data))
summary(lm(ppt~tmin*tmax,data = clim.data))
mantel(ppt.d~temp.d, nperm = 10000)

### Parsing a structural model using mantel
### geo.centroidd -> temp.D 
###        |           |    \_-->  MASH
###        |           V    _--->   D
###         \------> ppt.D /        ^
###                                 |  
###                               unkown
### geo -> climate 
set.seed(12345)
clim.d_geo.cd <- ecodist::mantel(clim.d~geo.cd)
temp.d_geo.cd <- ecodist::mantel(temp.d ~ geo.cd, nperm = 10000)
ppt.d_geo.cd <- ecodist::mantel(ppt.d ~ geo.cd, nperm = 10000)
mash.d_geo.cd_temp.d_ppt.d <- ecodist::mantel(
    mash.d ~ geo.cd + 
        temp.d + 
            ppt.d, 
    nperm = 10000) # unkown
### geo -> MASH
mash.d_geo.cd <- ecodist::mantel(mash.d ~ geo.cd, nperm = 10000)
### temp -> ppt
ppt.d_temp.d_geo.cd <- ecodist::mantel(ppt.d ~ temp.d + geo.cd, nperm = 10000)
###  ppt -> temp
temp.d_ppt.d_geo.cd <- ecodist::mantel(temp.d ~ ppt.d + geo.cd, nperm = 10000)
### clim -> MASH
mash.d_clim.d_geo.cd <- ecodist::mantel(mash.d ~ clim.d + geo.cd, nperm = 10000)
### temp.d -> MASH
mash.d_temp.d_ppt.d_geo.cd <- ecodist::mantel(
    mash.d ~ temp.d + ppt.d + geo.cd, nperm = 10000)
### ppt.d -> MASH

### rm rudis1
## backup <- FALSE
## if (backup){
##     mash.d.. <- mash.d
##     ppt.d.. <- ppt.d
##     temp.d.. <- temp.d
##     geo.cd.. <- geo.cd
## }


## mash.d <- mash.d..
## ppt.d <-ppt.d..
## temp.d <- temp.d..
## geo.cd <-geo.cd..

## rm.rud1 <- TRUE
## if (rm.rud1){
##     mash.d <- as.matrix(mash.d)
##     mash.d <- mash.d[rownames(mash.d) != 'Aphaenogaster rudis2',colnames(mash.d) != 'Aphaenogaster rudis2']
##     mash.d <- as.dist(mash.d)
##     ppt.d <- as.matrix(ppt.d)
##     ppt.d <- ppt.d[rownames(ppt.d) != 'rud6',colnames(ppt.d) != 'rud6']
##     ppt.d <- as.dist(ppt.d)
##     temp.d <- as.matrix(temp.d)
##     temp.d <- temp.d[rownames(temp.d) != 'rud6',colnames(temp.d) != 'rud6']
##     temp.d <- as.dist(temp.d)
##     geo.cd <- as.matrix(geo.cd)
##     geo.cd <- geo.cd[rownames(geo.cd) != 'rud6',colnames(geo.cd) != 'rud6']
##     geo.cd <- as.dist(geo.cd)
## }

mash.d_ppt.d_temp.d_geo.cd <- ecodist::mantel(mash.d ~ ppt.d + temp.d + geo.cd, nperm = 10000)
mantel.path <- list(clim.d_geo.cd = clim.d_geo.cd, 
                    temp.d_geo.cd = temp.d_geo.cd, 
                    ppt.d_geo.cd = ppt.d_geo.cd, 
                    mash.d_geo.cd_temp.d_ppt.d = mash.d_geo.cd_temp.d_ppt.d, 
                    mash.d_geo.cd = mash.d_geo.cd, 
                    ppt.d_temp.d_geo.cd = ppt.d_temp.d_geo.cd, 
                    temp.d_ppt.d_geo.cd = temp.d_ppt.d_geo.cd, 
                    mash.d_clim.d_geo.cd = mash.d_clim.d_geo.cd, 
                    mash.d_temp.d_ppt.d_geo.cd = mash.d_temp.d_ppt.d_geo.cd, 
                    mash.d_ppt.d_temp.d_geo.cd = mash.d_ppt.d_temp.d_geo.cd)
mantel.path <- do.call(rbind,mantel.path)[,1:4]
mash.path <- list(r = matrix(NA, nrow = 4, ncol = 4),
                  p = matrix(NA, nrow = 4, ncol = 4))
rownames(mash.path[[1]]) <- colnames(mash.path[[1]]) <- 
    rownames(mash.path[[2]]) <- colnames(mash.path[[2]]) <- c("Geographic Distance", 
                                                "Temperature Difference",
                                                "Precipitation Difference",
                                                "Genomic Distance (MASH)")
mash.path[["r"]]["Geographic Distance","Temperature Difference"] <- mantel.path["temp.d_geo.cd","mantelr"]
mash.path[["r"]]["Geographic Distance","Precipitation Difference"] <- mantel.path["ppt.d_geo.cd","mantelr"]
mash.path[["r"]]["Temperature Difference","Precipitation Difference"] <- mantel.path["ppt.d_temp.d_geo.cd",
                                                                                     "mantelr"]
mash.path[["r"]]["Precipitation Difference","Temperature Difference"] <- mantel.path["temp.d_ppt.d_geo.cd",
                                                                                     "mantelr"]
mash.path[["r"]]["Precipitation Difference","Genomic Distance (MASH)"] <- mantel.path["mash.d_ppt.d_temp.d_geo.cd",
                                                                                     "mantelr"]
mash.path[["r"]]["Temperature Difference","Genomic Distance (MASH)"] <- mantel.path["mash.d_temp.d_ppt.d_geo.cd",
                                                                                     "mantelr"]
mash.path[["p"]]["Geographic Distance","Temperature Difference"] <- mantel.path["temp.d_geo.cd","pval3"]
mash.path[["p"]]["Geographic Distance","Precipitation Difference"] <- mantel.path["ppt.d_geo.cd","pval3"]
mash.path[["p"]]["Temperature Difference","Precipitation Difference"] <- mantel.path["ppt.d_temp.d_geo.cd",
                                                                                     "pval3"]
mash.path[["p"]]["Precipitation Difference","Temperature Difference"] <- mantel.path["temp.d_ppt.d_geo.cd",
                                                                                     "pval3"]
mash.path[["p"]]["Precipitation Difference","Genomic Distance (MASH)"] <- mantel.path["mash.d_ppt.d_temp.d_geo.cd",
                                                                                     "pval3"]
mash.path[["p"]]["Temperature Difference","Genomic Distance (MASH)"] <- mantel.path["mash.d_temp.d_ppt.d_geo.cd",
                                                                                     "pval3"]
path.ig <- mash.path
mash.path <- lapply(mash.path, round, digits = 3)
for (i in 1:nrow(mash.path[["r"]])){
    for (j in 1:ncol(mash.path[["r"]])){
        if (is.na(mash.path[["r"]][i,j])){mash.path[["r"]][i,j] <- ""}else{
            mash.path[["r"]][i,j] <- paste0(mash.path[["r"]][i,j],", ",mash.path[["p"]][i,j])
        }
    }
}
mash.path.xtab <- xtable::xtable(mash.path[["r"]])
## Table: create mash mantel path analysis
print(mash.path.xtab,
      type = "latex",
      file = "results/mash_path.tex",
      include.rownames = TRUE,
      include.colnames = TRUE
      )


## Path diagram
mash.A <- round(path.ig[["r"]],3)
mash.A[is.na(mash.A)] <- 0
mash.Ap <- round(path.ig[["p"]],3)
mash.Ap[is.na(mash.Ap)] <- 0
ig <- graph_from_adjacency_matrix(mash.A, weighted = TRUE)
ig <- igraph.to.graphNEL(ig)
ig.p <- graph_from_adjacency_matrix(mash.Ap, weighted = TRUE)
ig.p <- igraph.to.graphNEL(ig.p)
attr <- list(node = list(shape = "box"))
ew.r <- as.character(unlist(edgeWeights(ig)))
ew.r <- ew.r[setdiff(seq(along = ew.r), removedEdges(ig))]
ec.p <- unlist(edgeWeights(ig.p))
ec.p[unlist(edgeWeights(ig.p)) < 0.1] <- "black"
ec.p[unlist(edgeWeights(ig.p)) >= 0.1] <- "grey"
ec.p <- as.character(ec.p)
names(ec.p) <- names(ew.r) <- edgeNames(ig)
attr.e <- list(label = ew.r, color = ec.p)

pdf("results/mash_path.pdf",height = 5, width = 5)
plot(ig, attrs = attr, edgeAttrs = attr.e)
dev.off()

### Climate table
clim.tab <- clim.data
rownames(clim.tab) <- clim.tab[,"mypoints.id"] <- rownames(as.matrix(mash.d))
clim.tab <- clim.tab[order(clim.tab[,"mypoints.id"]),]
clim.tab <- clim.tab[,c(2,1,6,5,4)]
colnames(clim.tab) <- c("Lat","Lon","Tmin (C)","Tmax (C)","Precip (mm)")

clim.xtab <- xtable::xtable(clim.tab, caption =
"Climate variables for colony sample sites. Climate are 30 year normal values (1976-2016) for January minimum temperature (Tmin), July maximum temperature (Tmax) and total precipitation (Precip).")
print(clim.xtab,
      type = "latex",
      file = "results/climate.tex",
      sanitize.rownames.function = italic,
      include.rownames = TRUE,
      include.colnames = TRUE
      )


### system("scp results/mash_path.pdf matthewklau@fas.harvard.edu:public_html")


### Update figures in presentations and manuscripts
## system("cp results/*.png docs/esa2017")
# system("cp results/*.png docs/manuscript")

print("Done!")
