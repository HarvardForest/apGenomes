library(vegan)
library(ecodist)


## Climate analyses
rownames(site.clim)[rownames(site.clim) == "arudis1"] <- "rud1"

clim.d <- dist(site.clim[,4:6])
temp.d <- dist(site.clim[,5:6])
max.d <- dist(site.clim[,"tmax"])
min.d <- dist(site.clim[,"tmin"])
ppt.d <- dist(site.clim[,"ppt"])

m.d <- mash
rownames(m.d)[rownames(m.d) == "rud2"] <- "rud6"
colnames(m.d)[colnames(m.d) == "rud2"] <- "rud6"
m.d <- m.d[match(rownames(as.matrix(clim.d)),rownames(m.d)),match(rownames(as.matrix(clim.d)),rownames(m.d))]
m.d <- as.dist(m.d)

## Mantels
vegan::mantel(clim.d,m.d,perm = 5000)
vegan::mantel(temp.d,m.d,perm = 5000)
vegan::mantel(ppt.d,m.d,perm = 5000)
vegan::mantel(max.d,m.d,perm = 5000)
vegan::mantel(min.d,m.d,perm = 5000)

## Ordination and Vectors
nms.ord <- ecodist::nmds(m.d,nits=500)
ord <- ecodist::nmds.min(nms.ord)

## Plot
ord.gg <- ggplot(data = ord,
                 map = aes(X1,X2,label = rownames(as.matrix(m.d)))) + 
    geom_text()



