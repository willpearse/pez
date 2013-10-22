plotLinesphylo<-function (x, y, ord, assoc.x, assoc.y, assoc.z, spb=10, bnk.cols="red",use.edge.length = TRUE, space = 0, length.line=1, gap=2, type="phylogram",
                          col.x="black",col.y="black",col.z="red",lwd.x=1,lwd.y=1,lwd.z=1,lty.x=1,lty.y=1,lty.z=1,return=TRUE, font=3,show.tip.label=FALSE)
{
  res <- list()
  left <- max(nchar(x$tip.label, type = "width")) + length.line
  right <- max(nchar(y$tip.label, type = "width")) + length.line
  space.min <- left + right + gap * 2
  if ((space <= 0) || (space < space.min))
    space <- space.min
  N.tip.x <- Ntip(x)
  N.tip.y <- Ntip(y)
  res$N.tip.x <- N.tip.x
  res$N.tip.y <- N.tip.y
  a <- plotPhyloCoor(x, use.edge.length = use.edge.length,
                     type = type)
  res$a <- a
  b <- plotPhyloCoor(y, use.edge.length = use.edge.length,
                     direction = "leftwards", type = type)
  a[, 2] <- a[, 2] - min(a[, 2])
  b[, 2] <- b[, 2] - min(b[, 2])
  res$b <- b
  b2 <- b
  b2[, 1] <- b[1:nrow(b), 1] * (max(a[, 1])/max(b[, 1])) +
    space + max(a[, 1])
  b2[, 2] <- b[1:nrow(b), 2] * (max(a[, 2])/max(b[, 2]))
  res$b2 <- b2
  c <- matrix(ncol = 2, nrow = nrow(a) + nrow(b))
  c[1:nrow(a), ] <- a[1:nrow(a), ]
  c[nrow(a) + 1:nrow(b), 1] <- b2[, 1]
  c[nrow(a) + 1:nrow(b), 2] <- b2[, 2]
  res$c <- c
  plot(rbind(c,apply(c,2,max)+spb,apply(c,2,min)-spb), type = "n", xlim = NULL, ylim = NULL, log = "", main = NULL,
       sub = NULL, xlab = NULL, ylab = NULL, ann = FALSE, axes = FALSE,
       frame.plot = FALSE)
  if (type == "cladogram") {
    for (i in 1:(nrow(a) - 1)) segments(a[x$edge[i, 1], 1],
                                        a[x$edge[i, 1], 2], a[x$edge[i, 2], 1], a[x$edge[i,
                                                                                         2], 2], col = "red")
    for (i in 1:(nrow(b) - 1)) segments(b2[y$edge[i, 1],
                                           1], b2[y$edge[i, 1], 2], b2[y$edge[i, 2], 1], b2[y$edge[i,
                                                                                                   2], 2])
  }
  if (type == "phylogram") {
    for (i in (N.tip.x + 1):nrow(a)) {
      l <- length(x$edge[x$edge[, 1] == i, ][, 1])
      for (j in 1:l) {
        segments(a[x$edge[x$edge[, 1] == i, ][1, 1],
                   1], a[x$edge[x$edge[, 1] == i, 2], 2][1], a[x$edge[x$edge[,
                                                                             1] == i, ][1, 1], 1], a[x$edge[x$edge[, 1] ==
                                                                                                              i, 2], 2][j])
        segments(a[x$edge[x$edge[, 1] == i, ][1, 1],
                   1],   a[x$edge[x$edge[, 1] == i, 2], 2][j], a[x$edge[x$edge[,
                                                                               1] == i, 2], 1][j], a[x$edge[x$edge[, 1] ==
                                                                                                              i, 2], 2][j])
      }
    }
    for (i in (N.tip.y + 1):nrow(b)) {
      l <- length(y$edge[y$edge[, 1] == i, ][, 1])
      for (j in 1:l) {
        segments(b2[y$edge[y$edge[, 1] == i, ][1, 1],
                    1], b2[y$edge[y$edge[, 1] == i, 2], 2][1],
                 b2[y$edge[y$edge[, 1] == i, ][1, 1], 1], b2[y$edge[y$edge[,
                                                                           1] == i, 2], 2][j])
        segments(b2[y$edge[y$edge[, 1] == i, ][1, 1],
                    1], b2[y$edge[y$edge[, 1] == i, 2], 2][j],
                 b2[y$edge[y$edge[, 1] == i, 2], 1][j], b2[y$edge[y$edge[,
                                                                         1] == i, 2], 2][j])
      }
    }
  }
  if (show.tip.label) {
    text(a[1:N.tip.x, ], cex = 0, font = font, pos = 4, labels = x$tip.label)
    text(b2[1:N.tip.y, ], cex = 1, font = font, pos = 2,
         labels = y$tip.label)
  }
  
  lsa <- 1:N.tip.x
  lsb <- 1:N.tip.y
  
  a.x<-(a[1,1])
  a.y<-rep(a[lsa,2])
  
  b.x<-(b2[1,1])
  b.y<-rep(b2[lsb,2])
  
  min.y<-min(a[,2])-spb
  max.y<-max(a[,2])+spb
  
  r.y<-seq(from=min.y,to=max.y,by=(max.y-min.y)/(length(ord)-1))
  r.y<-sort(r.y,decreasing=T)
  r.x<-(a[1,1]+space/2)
  
  
  if (length(col.x) == 1)
    colors.a <- c(rep(col.x, sum(assoc.x)))
  else
    colors.a <- col.x
  if (length(col.y) == 1)
    colors.b <- c(rep(col.y, sum(assoc.y)))
  else
    colors.a <- col.y
  colors.z <- c(rep(col.z, sum(assoc.z)))
  
  if (length(lwd.x) == 1)
    lwidths.a <- c(rep(lwd.x, sum(assoc.x)))
  else
    lwidths.a <- lwd.x
  if (length(lwd.y) == 1)
    lwidths.b <- c(rep(lwd.y, sum(assoc.y)))
  else
    lwidths.b <- lwd.y
  
  lwidths.z <- lwd.z
  
  if (length(lty.x) == 1)
    ltype.a <- c(rep(lty.x, sum(assoc.x)))
  else
    ltype.a <- lty.x
  if (length(lty.y) == 1)
    ltype.b <- c(rep(lty.y, sum(assoc.y)))
  else
    ltype.b <- lty.y
  
  ltype.z <- lty.z
  
  if(length(bnk.cols)==1)
    bnk.cols<-rep(bnk.cols,length(ord))
  
  #plot the segments X past
  for (i in colnames(assoc.x)) {
    if (show.tip.label) {
      decx <- strwidth(i)
    }
    else {
      decx <- 0
    }
    hold<-assoc.x[,i]
    for (j in names(hold)[hold==1])
    {
      segments(a.x + decx + gap, a.y[x$tip.label == i], r.x, r.y[ord == j], col = bnk.cols[ord == j], lwd = lwidths.a,lty = ltype.a)
    }
  }
  
  #plot the segments Y PRESENT
  for (i in colnames(assoc.y)) {
    if (show.tip.label) {
      decy <- strwidth(i)
    }
    else {
      decy <- 0
    }
    hold<-assoc.y[,i]
    for (j in names(hold)[hold==1])
    {
      segments(b.x - (decy + gap), b.y[y$tip.label == i], r.x, r.y[ord == j], col = colors.b, lwd = lwidths.b,lty = ltype.b)
    }
  }
  
  #plot the segments Z PRESENT
  for (i in colnames(assoc.z)) {
    if (show.tip.label) {
      decy <- strwidth(i)
    }
    else {
      decy <- 0
    }
    hold<-assoc.z[,i]
    for (j in names(hold)[hold==1])
    {
      segments(b.x - (decy + gap), b.y[y$tip.label == i], r.x, r.y[ord == j], col = bnk.cols[ord == j], lwd = lwidths.z,lty = ltype.z)
    }
  }
  
  
  
  hilight <- function(x,y,s, bg="white",border="black",col="red", cex=1) {
    text.width <- strwidth(s)*(cex+cex/10)
    text.height <- strheight(s)*(cex+cex/10)
    rect(x-text.width/2,y-text.height/2,x+text.width/2,y+text.height/2,col=bg,border=border,cex=cex)
    text(x,y,s,cex=cex,col=col)
  }
  
  
  hilight(r.x,r.y+(spb/4),ord,cex=.7,col="black",bg=rgb(1,1,1,.9),border=NA)
  points(rep(r.x,length(r.y)),r.y,pch=21,bg=bnk.cols,col="black",cex=1)
  
  x.past<-a.x
  y.past<-r.y[2]
  x.pres<-b.x
  text(x.past,y.past,"Native",pos=2,cex=1.1)
  text(x.pres,y.past,"Exotic",pos=4,cex=1.1)
  #legend(x.pres,r.y[length(r.y)-2]-.5,legend=c("native","established exotic"),lwd=c(lwd.y,lwd.z),lty=c(lty.y,lty.z),col=c(col.y,col.z),cex=0.7,bty="n")
  
  s.bar<-a[1,1]
  e.bar<-s.bar-s.bar/10
  y.bar<-r.y[length(r.y)-1]
  arrows(s.bar,y.bar,e.bar,y.bar,angle=90,length=0.02,code=3)
  text(e.bar,y.bar,"10 MYA",pos=2,cex=0.7)
  
  #hilight(b.x+13,b.y[y$tip.label == "sagrei"],"sagrei",cex=.7,col="black",bg=rgb(1,1,1,1),border="black")
  
  #text(,pch=21, col="black",bg="red")
  #points(b.x,b.y[y$tip.label == "sagrei"],pch=21, col="black",bg="red")
  if (return == TRUE)
    return(res)
}


####################################################################################################################################################################
################ Plots two trees, with lines conecting to a central point
##
cophyloplotB<-function (x, y, ord, assoc.x, assoc.y, assoc.z=assoc.z, bnk.cols="black", spb=10, col.z="red", lwd.z=2,lty.z=1, use.edge.length = FALSE, space = 0,
                        length.line = 1, gap = 2, type = "phylogram", rotate = TRUE, col.x="black",col.y="black",lwd.x=1,lwd.y=1,lty.x=1,lty.y=1,
                        col = par("fg"), lwd = par("lwd"), lty = par("lty"), show.tip.label = FALSE, font = 3, ...)
{
  if (rotate == TRUE) {
    cat("\n    Click on a node to rotate (right click to exit)\n\n")
    repeat {
      res <- plotLinesphylo(x=x, y=y, ord=ord, assoc.x=assoc.x, assoc.y=assoc.y, assoc.z=assoc.z, bnk.cols=bnk.cols, spb=spb,col.z=col.z, lwd.z=lwd.z,lty.z=lty.z, use.edge.length = use.edge.length,
                            space = space, length.line=length.line, gap=gap, type=type, show.tip.label=show.tip.label,col.x=col.x,col.y=col.y,
                            lwd.x=lwd.x,lwd.y=lwd.y,lty.x=lty.x,lty.y=lty.y,return=TRUE, font=font)
      
      click <- identify(res$c[, 1], res$c[, 2], n = 1)
      x <- rotate(x, click)
      y <- rotate(y, click)
      plotLinesphylo(x=x, y=y, ord=ord, assoc.x=assoc.x, assoc.y=assoc.y, assoc.z=assoc.z,bnk.cols=bnk.cols, spb=spb, col.z=col.z, lwd.z=lwd.z,lty.z=lty.z, use.edge.length = use.edge.length, space = space,
                     length.line=length.line, gap=gap, type=type, show.tip.label=show.tip.label,col.x=col.x,col.y=col.y,lwd.x=lwd.x,lwd.y=lwd.y,
                     lty.x=lty.x,lty.y=lty.y,return=TRUE,font=font)
    }
    on.exit(print("done"))
  } else {plotLinesphylo(x=x, y=y, ord=ord, assoc.x=assoc.x, assoc.y=assoc.y, assoc.z=assoc.z,bnk.cols=bnk.cols, spb=spb, col.z=col.z, lwd.z=lwd.z,lty.z=lty.z, use.edge.length = use.edge.length,
                         space = space, length.line=length.line, gap=gap, type=type, show.tip.label=show.tip.label,col.x=col.x,col.y=col.y,
                         lwd.x=lwd.x,lwd.y=lwd.y,lty.x=lty.x,lty.y=lty.y,return=TRUE,font=font)}
}
