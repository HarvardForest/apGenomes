### https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-8-64
### Genome sizes were analyzed by flow cytometry

library(XML)
library(RCurl)
library(rlist)
theurl <- getURL("https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-8-64",.opts = list(ssl.verifypeer = FALSE) )
tables <- readHTMLTable(theurl)
tables <- list.clean(tables, fun = is.null, recursive = FALSE)
ant.gen.size <- tables[[1]]
ant.gen.size <- ant.gen.size[!apply(ant.gen.size,1,function(x) any(grepl('MEAN',x))),]
ant.gen.size <- ant.gen.size[ant.gen.size[,1] == '',]
ant.gen.size[,'1C Genome Size (Mb)'] <- as.numeric(as.character(ant.gen.size[,'1C Genome Size (Mb)']))

### Genome Size by Location
ant.gs.world <- table(substr(ant.gen.size[,"Collection Info"],1,3))
ant.gs.usa <- table(grep("USA",ant.gen.size[,"Collection Info"],value = TRUE))
names(ant.gs.usa) <- gsub("USA: ","",names(ant.gs.usa))
ags.w <- ant.gs.world[order(ant.gs.world,decreasing = TRUE)]
ags.w <- data.frame(region = names(ags.w),count = as.numeric(ags.w))
ags.u <- ant.gs.usa[order(ant.gs.usa,decreasing = TRUE)]
ags.u <- data.frame(region = names(ags.u),count = as.numeric(ags.u))



### From the GAGA group
gaga <- na.omit(read.csv("data/gaga_genome_info.csv"))
gaga[,"Scaffold.N50.length.kb."] <- as.numeric(gsub(",","",as.character(
    gaga[,"Scaffold.N50.length.kb."])))
gaga[,"Contig.N50.length.kb."] <- as.numeric(gsub(",","",as.character(
    gaga[,"Contig.N50.length.kb."])))
gaga.loc <- do.call(rbind,strsplit(as.character(gaga[,"Location"]),","))
gaga.loc <- apply(gaga.loc,2,gsub,pattern = "\001",replacement = "")
gaga.loc <- apply(gaga.loc,2,gsub,pattern = " ",replacement = "")
colnames(gaga.loc) <- c("City","State","Country")
gaga.loc <- data.frame(gaga.loc)
