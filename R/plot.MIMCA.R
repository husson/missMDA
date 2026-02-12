plot.MIMCA<-function (x, choice = "all", axes = c(1, 2), new.plot = TRUE, 
                      main = NULL, level.conf = 0.95, graph.type=c("ggplot","classic"),...){
 procrustes <- function(amat, target, orthogonal = FALSE, translate = FALSE,
        magnify = FALSE) {
#        for (i in nrow(amat):1) {
#            if (any(is.na(amat)[i, ]) | any(is.na(target)[i,
#                ])) {
#                amat <- amat[-i, ]
#                target <- target[-i, ]
#            }
#        }
#        dA <- dim(amat)
        dX <- dim(target)
#       if (length(dA) != 2 || length(dX) != 2)
#            stop("arguments amat and target must be matrices")
#        if (any(dA != dX))
#            stop("dimensions of amat and target must match")
#        if (length(attr(amat, "tmat")))
#            stop("oblique loadings matrix not allowed for amat")       
if (orthogonal) {
            if (translate) {
                p <- dX[1]
                target.m <- colMeans(target)
                amat.m <- colMeans(amat)
                target.c <- scale(target, center = target.m,scale = FALSE)
                amat.c <- scale(amat, center = amat.m, scale = FALSE)
                j <- svd(crossprod(target.c, amat.c))
            }
#			else {
#                amat.c <- amat
#                j <- svd(crossprod(target, amat))
#            }
            rot <- j$v %*% t(j$u)
#            if (magnify)
                beta <- sum(j$d)/sum(amat.c^2)
#            else beta <- 1

            B <- beta * amat.c %*% rot
            if (translate) B <- B + rep(as.vector(target.m), rep.int(p,dX[2]))   
            value <- list(rmat = B, tmat = rot, magnify = beta)
            if (translate) value$translate <- target.m - (rot %*% amat.m)
    
  }
#  else {
#            b <- solve(amat, target)
#            gamma <- sqrt(diag(solve(crossprod(b))))
#            rot <- b * rep(gamma, rep.int(dim(b)[1], length(gamma)))
#            B <- amat %*% rot
#            fcor <- solve(crossprod(rot))
#            value <- list(rmat = B, tmat = rot, correlation = fcor)
#        }
        return(value)
    }

  
  coord.sup<-function(G,Gsup,ncp,p,res.mca){
    #G tdc des indiv actifs
    #Gsup tdc des illustratifs
    #p nb de variables quali
    
    n <- nrow(G)
    Dc <- drop((rep(1, n)) %*% G)
    X <- t(t(G)/(sqrt(p * Dc)))
    Xsup <- t(t(Gsup)/(sqrt(p * Dc)))
    X.svd <- svd(X)
    sec <- 1 + (1L:ncp)
    rs <- sweep(Xsup %*% X.svd$v[, sec],1,STATS=X.svd$u[1:nrow(Gsup),1,drop=F],FUN="/")*(sweep(X %*% X.svd$v[, sec],1,STATS=X.svd$u[1:nrow(Gsup),1,drop=F],FUN="/")/res.mca$ind$coord[,1:ncp])
	return(rs)
  }
  
  plot.tmp<-function(res.pca,axes,col.ind.sup,label,col.quali,ellipse,title,invisible, new.plot,graph.type){
    #plot PCA pour modifier label axes
    if (graph.type=="classic"){
  	  par(col.lab="white")
      FactoMineR::plot.PCA(res.pca, axes = 1:2, col.ind.sup = col.ind.sup, label = label, col.quali = col.quali, ellipse = ellipse,
             title = title, invisible = invisible, new.plot = new.plot,graph.type=graph.type)
	  par(col.lab="black")
      mtext(paste0("Dim ", axes[1], " (", format(res.pca$eig[1,2], nsmall = 2, digits = 2), "%)"), side=1, line=3)
      mtext(paste0("Dim ", axes[2], " (", format(res.pca$eig[2,2], nsmall = 2, digits = 2), "%)"), side=2, line=3)
	  } else{
        p <- FactoMineR::plot.PCA(res.pca, axes = 1:2, col.ind.sup = col.ind.sup, label = label, col.quali = col.quali,  ellipse=ellipse,
               title = title, invisible = invisible, new.plot = new.plot,graph.type=graph.type)
	    return(p)
	}
  }
  
  res <- x
  if (!inherits(res, "MIMCA")) stop("non convenient data")
  choice <- tolower(choice)
  choice <- match.arg(choice,c("all","ind.supp","ind.proc","dim","mod.supp"))
#  graph.type <- "classic" ## graph too long with ggplot
  graph.type <- match.arg(graph.type[1],c("ggplot","classic"))
  graph <- list()
  ncp <- max(axes)
  reference <- FactoMineR::MCA(res$call$X, graph = FALSE, ncp = ncp,tab.disj=res$res.imputeMCA)
  res.dim  <- res$call$X
  res.procrustes <- reference$ind$coord[, 1:ncp]
  
  for (i in 1:length(res$res.MI)){
    rec.mca <- res$res.MI[[i]]
    acmfin <- FactoMineR::MCA(res$call$X, graph = FALSE, ncp = ncp,tab.disj = res$call$tab.disj[,,i])
    tourne <- procrustes(acmfin$ind$coord[, 1:ncp], reference$ind$coord[,1:ncp],
                         orthogonal = TRUE, translate = TRUE, magnify = TRUE)$rmat
    colnames(tourne) <- colnames(res.procrustes)
    res.procrustes <- rbind.data.frame(res.procrustes, tourne)
    res.dim <- cbind.data.frame(res.dim, acmfin$ind$coord[,1:ncp])
  }

  if (("all"%in%choice) | ("ind.supp"%in%choice) | ("mod.supp"%in%choice)) {
    coordinsup<-apply(res$call$tab.disj,3,FUN=function(tab.imp,don.na,axes,G,res.mca){
      p <-ncol(don.na)
      res.out<-coord.sup(G,Gsup=tab.imp,ncp=max(axes),p=p,res.mca=res.mca)
      return(list(res.out))
    },don.na=res$call$X,axes= axes,G=res$res.imputeMCA,res.mca=reference)
    coordinsup<-lapply(coordinsup,"[[",1)
  }

  if (!is.null(main)) title <- main
  if (("all"%in%choice) | ("ind.proc"%in%choice)) {
    if ((new.plot)&!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) dev.new()
    oo = FactoMineR::PCA(res.procrustes, ind.sup = c((nrow(res$call$X) + 1):nrow(res.procrustes)), scale.unit = FALSE, graph = FALSE)
    oo$eig = reference$eig
    if (is.null(main)){title = "Multiple imputation using Procrustes"}
 	el = coord.ellipse(cbind.data.frame(as.factor(rep(rownames(res$call$X), res$call$nboot)), oo$ind.sup$coord[, axes]), level.conf = level.conf)
    if (graph.type=="classic"){
	  plot(oo, axes = axes, col.ind.sup = rep(1:nrow(res$call$X), res$call$nboot), label = "ind", ellipse = el, col.quali = 1, title = title, invisible = "ind.sup", new.plot = FALSE,graph.type=graph.type)
    } else {
	  graph$PlotIndProc <- plot(oo, axes = axes, col.ind.sup = rep(1:nrow(res$call$X), res$call$nboot), label = "ind", ellipse = el, col.quali = 1, title = title, invisible = "ind.sup", new.plot = FALSE,graph.type=graph.type)
#	  graph$PlotIndProc <- plot(reference, axes = axes, label = "ind", invisible=c("quali","quali.sup"), title = title, new.plot = FALSE)
#	  graph$PlotIndProc <- graph$PlotIndProc + stat_ellipse(aes(x=oo$ind.sup$coord[, axes[1]],y=oo$ind.sup$coord[, axes[2]],color=as.factor(rep(rownames(res$call$X), res$call$nboot))))
	}
  }

  if (("all"%in%choice) | ("dim"%in%choice)) {
    if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) {dev.new()}
    colnames(res.dim) = paste("V", 1:ncol(res.dim))
    ooo = FactoMineR::MCA(res.dim, quanti.sup = (ncol(res$call$X) + 1):ncol(res.dim), graph = FALSE,tab.disj = cbind(res$res.imputeMCA,res.dim[,-(1:ncol(res$call$X))]),ncp = ncp)
    ooo$eig = reference$eig
    if (is.null(main)){title <- "Projection of the Principal Components"}
    graph$PlotDim <- plot(ooo, choi = "quanti.sup", axes = axes, title = title, label = "none", new.plot = FALSE, graph.type=graph.type)
  }
  
  if (("all"%in%choice) | ("ind.supp"%in%choice)) {
    if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) {dev.new()}
    if (is.null(main)){title = "Supplementary projection\n individuals"}
    oo<-FactoMineR::PCA(rbind(reference$ind$coord[, axes],do.call(rbind,coordinsup)[, axes]),
            ind.sup = c((nrow(res$call$X) + 1):(nrow(res$call$X)*(res$call$nboot+1))), 
            scale.unit = FALSE, graph = FALSE,ncp=2)
    oo$eig = reference$eig[axes,]
    oo$var<-lapply(reference$var,"[",,axes)
    
	el = coord.ellipse(cbind.data.frame(as.factor(rep(rownames(res$call$X), res$call$nboot)),oo$ind.sup$coord), level.conf = level.conf)
    if (graph.type=="classic"){
      plot.tmp(oo, axes = axes, col.ind.sup = rep(1:nrow(res$call$X), res$call$nboot),
        label = "ind", col.quali = 1, ellipse = el, title = title, invisible = c("ind.sup","var"), new.plot = FALSE, graph.type=graph.type)
     } else {
	   oo$ind.sup=NULL
	   graph$PlotIndSupp <- plot.tmp(oo, axes = axes, col.ind.sup = rep(1:nrow(res$call$X), res$call$nboot),
        label = "ind", col.quali = 1, ellipse = el, title = title, invisible = c("ind.sup","var"), new.plot = FALSE, graph.type=graph.type)
	}

  }
  
  if (("all"%in%choice) | ("mod.supp"%in%choice)) {
    if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) dev.new()

    coordmodsup<-sapply(1:length(coordinsup),FUN=function(indic,coord,tabdisj,don){
      res.mca<-FactoMineR::MCA(res$call$X,graph=F,tab.disj = tabdisj[,,indic])
      res.out<-sapply(1:ncol(coord[[1]]),FUN=function(ii,res.mca,tabdisj,coord){
        res.out<-sqrt(1/res.mca$eig[ii,1])*diag(1/colSums(tabdisj))%*%t(tabdisj)%*%coord[,ii]
        return(res.out)},res.mca=res.mca,tabdisj=tabdisj[,,indic],coord=coord[[indic]])
      return(res.out)
    },coord=coordinsup,tabdisj=res$call$tab.disj,don=res$call$X,simplify = FALSE)
    
    if (is.null(main)){title <- "Supplementary projection\ncategories"}
    
    oo<-FactoMineR::PCA(rbind(reference$var$coord[, axes],do.call(rbind,coordmodsup)[,axes]),
            ind.sup = c((nrow(reference$var$coord) + 1):(nrow(reference$var$coord)*(res$call$nboot+1))), 
            scale.unit = FALSE, graph = FALSE)
    oo$eig = reference$eig[axes,]
    oo$var<-lapply(reference$var,"[",,axes)
    
    el = coord.ellipse(cbind.data.frame(
      as.factor(rep(rownames(reference$var$coord), res$call$nboot)),
      oo$ind.sup$coord[, 1:2]), level.conf = level.conf)
    if (graph.type=="classic"){
	  plot.tmp(oo, axes = axes, col.ind.sup = rep(1:nrow(reference$var$coord), res$call$nboot),
          label = "ind", col.quali = 1, ellipse = el,title = title, invisible = c("ind.sup","var"), new.plot = FALSE, graph.type=graph.type)
    } else {
     graph$PlotModSupp <- plot.tmp(oo, axes = axes, col.ind.sup = rep(1:nrow(reference$var$coord), res$call$nboot),
          label = "ind", col.quali = 1, ellipse = el,title = title, invisible = c("ind.sup","var"), new.plot = FALSE,graph.type=graph.type)
    }
  }
   if (graph.type=="ggplot") return(graph)
}

