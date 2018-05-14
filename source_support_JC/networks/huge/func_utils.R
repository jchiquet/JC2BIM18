mse <- function(theta, theta.star) {
  return(lapply(1:length(theta),function(i){ mean((theta[[i]]-theta.star)^2, na.rm=TRUE) } ))
}

ham <- function(theta, theta.star) {
  return( sum( ((as.vector(theta)!=0)*1) != ((as.vector(theta.star)!=0)*1) , na.rm=TRUE) )
}

perf.mse <- function(theta.hat, theta.star) {

  if (is.list(theta.hat)) {
    return(sapply(theta.hat, mse, theta.star))
  } else {
    return(mse(theta.hat, theta.star))
  }
}

perf.ham <- function(theta.hat, theta.star) {

  if (is.list(theta.hat)) {
    return(sapply(theta.hat, ham, theta.star))
  } else {
    return(ham(theta.hat, theta.star))
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

perf.auc <- function(roc,ind=1) {
  fallout <- c(0,roc$fallout,ind)
  recall  <- c(0,roc$recall, ind)
  dx <- diff(fallout)
  return(sum(c(recall[-1]*dx, recall[-length(recall)]*dx))/2)
}
