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
