library(ggplot2)
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
perf.roc <- function(theta.hat, theta.star) {
  
  roc <- function(theta) {
    nzero <- which(theta != 0)
    zero  <- which(theta == 0)
    
    true.nzero <- which(theta.star != 0)
    true.zero  <- which(theta.star == 0)
    
    TP <- sum(nzero %in% true.nzero)
    TN <- sum(zero %in%  true.zero)
    FP <- sum(nzero %in% true.zero)
    FN <- sum(zero %in%  true.nzero)
    
    recall    <- TP/(TP+FN) ## also recall and sensitivity
    fallout   <- FP/(FP+TN) ## also 1 - specificit
    
    res <-  round(c(fallout,recall),3)
    res[is.nan(res)] <- 0
    names(res) <- c("fallout","recall")
    return(res)
  }
  
  if (is.list(theta.hat)) {
    return(as.data.frame(do.call(rbind, lapply(theta.hat, roc))))
  } else {
    return(roc(theta.hat))
  }
}

perf.auc <- function(roc) {
  fallout <- c(0,roc$fallout,1)
  recall  <- c(0,roc$recall, 1)
  dx <- diff(fallout)
  return(sum(c(recall[-1]*dx, recall[-length(recall)]*dx))/2)
}
