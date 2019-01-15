estim_ncpMultilevel <-function (X,  ifac=1, ncpW.min = 1, ncpW.max = 5, ncpB.min = 1, ncpB.max = 5,
         scale = TRUE, nbsim = 100, pNA = 0.05, threshold = 1e-04, nb.cores =NULL, verbose = TRUE) 
{
  
  tab.disjonctif.NA <- function(tab) {
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) {
      moda <- tab[, i]
      nom <- names(tab)[i]
      n <- length(moda)
      moda <- as.factor(moda)
      x <- matrix(0, n, length(levels(moda)))
      ind <- (1:n) + n * (unclass(moda) - 1)
      indNA <- which(is.na(ind))
      x[(1:n) + n * (unclass(moda) - 1)] <- 1
      x[indNA, ] <- NA
      if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), "n", "N", "y", "Y"))) 
        dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), sep = "."))
      else dimnames(x) <- list(row.names(tab), levels(moda))
      return(x)
    }
    if (ncol(tab) == 1) 
      res <- modalite.disjonctif(1)
    else {
      res <- lapply(1:ncol(tab), modalite.disjonctif)
      res <- as.matrix(data.frame(res, check.names = FALSE))
    }
    return(res)
  }
  prodna<-function (x, noNA){
    n <- nrow(x)
    p <- ncol(x)
    NAloc <- rep(FALSE, n * p)
    NAloc[sample(n * p, floor(n * p * noNA))] <- TRUE
    x[matrix(NAloc, nrow = n, ncol = p)] <- NA
    return(x)
  }
  if (is.null(nb.cores)) cl <- makeCluster(detectCores()-1)
  else cl <- makeCluster(min(nb.cores,detectCores()-1))
  registerDoParallel(cl)
  X <- as.data.frame(X)
  ncpW.max <- min(ncpW.max, ncol(X) - 2, nrow(X) - nlevels(X[, ifac])-1)
  ncpB.max <- min(ncpB.max, nlevels(X[, ifac])- 2, ncol(X) - 2)
  jeu<-X[,c(which(sapply(X,is.numeric)),setdiff(which(sapply(X,is.factor)),ifac)),drop=FALSE]
  nbquanti<-sum(sapply(X,is.numeric))
  tab.jeu <- NULL
  if (nbquanti>0) tab.jeu <- sapply(jeu[,1:nbquanti,drop=FALSE],as.double)
  if (nbquanti<ncol(jeu)) tab.jeu <- cbind.data.frame(tab.jeu,tab.disjonctif.NA(jeu[,(nbquanti+1):ncol(jeu),drop=FALSE]))

  opts=NULL
  if (verbose){
    pb <- txtProgressBar(min = 1/nbsim * 100, max = 100,  style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }
  res <- foreach (sim = 1:nbsim, .options.snow = opts, .export = c("imputePCA"), .packages = c("missMDA","FactoMineR"), .combine = rbind) %dopar%  { 
    
    continue<-TRUE
    while(continue){
      jeuNA <- prodna(jeu, pNA)
      continue<-    continue<- (sum(unlist(sapply(as.data.frame(jeuNA[,-c(1:nbquanti),drop=F]),nlevels)))!=sum(unlist(sapply(jeu,nlevels))))
    }
    auxi <- matrix(NA,ncpW.max - ncpW.min + 1,  ncpB.max - ncpB.min + 1)
    for (nbaxesB in ncpB.min:ncpB.max) {
      for (nbaxesW in ncpW.min:ncpW.max) {
        # if (nbaxesW == 0) {
        # Xhat <- XNA
        # for (j in 1:ncol(X)) Xhat[, j] <- replace(XNA[,
        #                                             j], is.na(XNA[, j]), mean(XNA[, j], na.rm = TRUE))
        # }
        Xhat <-imputeMultilevel(cbind.data.frame(X[,ifac],jeuNA), ifac = 1, ncpB = nbaxesB, ncpW = nbaxesW, scale = scale,  threshold = 1e-04, maxiter = 1000)$Xhat
        # res[nbaxesW - ncpW.min + 1, nbaxesB - ncpB.min + 1, sim] <- sum((Xhat[,-ifac] - X[,-ifac])^2, na.rm = TRUE)
        auxi[nbaxesW - ncpW.min + 1, nbaxesB - ncpB.min + 1]<- sum((Xhat - tab.jeu)^2, na.rm = TRUE)
      }
    }
    as.vector(auxi)
  }
  stopCluster(cl)
  if (verbose) close(pb)
  
  resu <- apply(res, 2, mean)
  resu <- matrix(resu,nrow=ncpW.max-ncpW.min+1)
  rownames(resu) <- paste("W", ncpW.min:ncpW.max)
  colnames(resu) <- paste("B",ncpB.min:ncpB.max)
  mini<- which(resu==min(resu),arr.ind=TRUE)+c(ncpW.min-1, ncpB.min-1)
  result <- list(criterion = resu, res = res, ncpW=mini[1],ncpB=mini[2])
  return(result)
}
