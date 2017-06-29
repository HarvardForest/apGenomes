### Plot mummer alignments

mumplot <- function(file,main,xlab,ylab,pch=19,cex=0.25,sd.thresh = 2){
    if (missing(main)){main <- gsub('.delta.dists','',file)}
    if (missing(xlab) | missing(ylab)){
        .file <- gsub('.delta.dists','',file)
        if (missing(xlab)){xlab <- strsplit(.file,'_')[[1]][1]}
        if (missing(ylab)){ylab <- strsplit(.file,'_')[[1]][2]}        
    }
    pos <- readLines(file)
    pos <- do.call(rbind,strsplit(pos,' '))
    pos <- apply(pos,2,as.numeric)
    pos.dif <- abs(apply(pos[,1:2],1,diff))
    pos <- pos[pos.dif >= (mean(pos.dif) + (sd.thresh * sd(pos.dif))),]
    plot(pos[,3:4],xlab = xlab,ylab = ylab, pch = pch,cex = cex,main = main,
         xlim = c(0,max(pos[,3:4])),ylim = c(0,max(pos[,3:4])))
    lines(c(0,max(pos[,3:4])),c(0,max(pos[,3:4])),
          lwd = 0.25,col = 'darkgrey')
}


mumspark <- function(file,main,xlab,ylab,pch=19,cex=0.25,sd.thresh = 2){
    if (missing(main)){main <- gsub('.delta.dists','',file)}
    if (missing(xlab) | missing(ylab)){
        .file <- gsub('.delta.dists','',file)
        if (missing(xlab)){xlab <- strsplit(.file,'_')[[1]][1]}
        if (missing(ylab)){ylab <- strsplit(.file,'_')[[1]][2]}        
    }
    pos <- readLines(file)
    pos <- do.call(rbind,strsplit(pos,' '))
    pos <- apply(pos,2,as.numeric)
    pos.dif <- apply(pos[,1:2],1,diff)
    pos <- pos[abs(pos.dif) >= (mean(pos.dif) + (sd.thresh * sd(pos.dif))),]
    pos.dif <- apply(pos[,1:2],1,diff)
    plot(pos.dif~pos[,1],xlab = xlab,ylab = ylab, 
         pch = pch,cex = cex,main = main,
         ylim = c(-max(pos.dif),max(pos.dif)))
    abline(h = 0,col = 'darkgrey',lty = 2)
}

