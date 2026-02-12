plot.MIFAMD<-function (x, choice = "all", axes = c(1, 2), new.plot = TRUE, 
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

  
  coord.sup<-function(G,Gsup,ncp,p,res.famd){
    #G tdc des indiv actifs
    #Gsup tdc des illustratifs
    G <- t(t(G)-res.famd$call$centre)
    ponderation <- c(res.famd$call$ecart.type[1:(length(res.famd$call$ecart.type)-length(res.famd$call$prop))]^2,res.famd$call$prop) ## ponderation = ecart-type OU sqrt(prop)
	X <- t(t(G)/sqrt(ponderation))
    X.svd <- svd(X)

    Gsup <- t(t(Gsup)-res.famd$call$centre)
	Xsup <- t(t(Gsup)/sqrt(ponderation))
#print(Xsup %*% X.svd$v[,1:2])
#print(res.famd$ind$coord)
# print(X %*% X.svd$v) # = print(res.famd$ind$coord)
#    sec <- 1 + (1L:ncp)
#print((X %*% X.svd$v[, 1:ncp])/res.famd$ind$coord[,1:ncp])
#    rs <- Xsup %*% X.svd$v[, sec]
    rs <- (Xsup %*% X.svd$v[, 1:ncp])*((X %*% X.svd$v[, 1:ncp])/res.famd$ind$coord[,1:ncp])
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
  if (!inherits(res, "MIFAMD")) stop("non convenient data")
  choice <- tolower(choice)
  choice <- match.arg(choice,c("all","ind.supp","ind.proc","dim","mod.supp","var"),several.ok =TRUE)
  graph.type <- match.arg(graph.type[1],c("ggplot","classic"))
  graph <- list()
  ncp <- max(axes)
  reference <- FactoMineR::FAMD(res$call$X, graph = FALSE, ncp = ncp,tab.disj=res$res.imputeFAMD$tab.disj)
  res.dim  <- res$call$X
  res.procrustes <- reference$ind$coord[, 1:ncp]
  
  for (i in 1:length(res$res.MI)){
    rec.famd <- res$res.MI[[i]]
    famdfin <- FactoMineR::FAMD(res$call$X, graph = FALSE, ncp = ncp,tab.disj = as.matrix(res$call$tab.disj[[i]]))
	tourne <- procrustes(famdfin$ind$coord[, 1:ncp], reference$ind$coord[,1:ncp],
                         orthogonal = TRUE, translate = TRUE, magnify = TRUE)$rmat
    colnames(tourne) <- colnames(res.procrustes)
    res.procrustes <- rbind.data.frame(res.procrustes, tourne)
    res.dim <- cbind.data.frame(res.dim, famdfin$ind$coord[,1:ncp])
  }

  if (("all"%in%choice) | ("ind.supp"%in%choice) | ("mod.supp"%in%choice)) {
    coordinsup<-lapply(res$call$tab.disj,FUN=function(tab.imp,don.na,axes,G,res.famd){
      p <-ncol(don.na)
      res.out<-coord.sup(G,Gsup=tab.imp,ncp=max(axes),p=p,res.famd=res.famd)
      return(list(res.out))
    },don.na=res$call$X,axes= axes,G=res$res.imputeFAMD$tab.disj,res.famd=reference)
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
	}
  }

  if (("all"%in%choice) | ("dim"%in%choice)) {
    if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) {dev.new()}
    colnames(res.dim) = paste0("V", 1:ncol(res.dim))
    colnames(res.dim)[1:ncol(res$call$X)] <- colnames(res$call$X)
    ooo = FactoMineR::FAMD(res.dim, sup.var = (ncol(res$call$X) + 1):ncol(res.dim), graph = FALSE,tab.disj = as.matrix(cbind(res$res.imputeFAMD$tab.disj,res.dim[,-(1:ncol(res$call$X))])),ncp = ncp)
    ooo$eig = reference$eig
    if (is.null(main)){title <- "Projection of the Principal Components"}
    graph$PlotDim <- plot(ooo, choi = "quanti", axes = axes, title = title, lab.var = FALSE, new.plot = FALSE, graph.type=graph.type)
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
	  res.famd<-FactoMineR::FAMD(res$call$X,graph=F,tab.disj = as.matrix(tabdisj[[indic]]))
      res.out<-sapply(1:ncol(coord[[1]]),FUN=function(ii,res.famd,tabdisj,coord){
        res.out <- diag(1/colSums(tabdisj[,-which(sapply(res$call$X,is.numeric))]))%*%t(tabdisj[,-which(sapply(res$call$X,is.numeric))])%*%coord[,ii]
        return(res.out)},res.famd=res.famd,tabdisj=tabdisj[[indic]],coord=coord[[indic]])
      return(res.out)
    },coord=coordinsup,tabdisj=res$call$tab.disj,don=res$call$X,simplify = FALSE)

    if (is.null(main)){title <- "Supplementary projection\ncategories"}
    oo<-FactoMineR::PCA(rbind(reference$quali.var$coord[, axes],do.call(rbind,coordmodsup)[,axes]),
            ind.sup = c((nrow(reference$quali.var$coord) + 1):(nrow(reference$quali.var$coord)*(res$call$nboot+1))), 
            scale.unit = FALSE, graph = FALSE)
    oo$eig = reference$eig[axes,]
    oo$quali.var<-lapply(reference$quali.var,"[",,axes)
    
    el = coord.ellipse(cbind.data.frame(
      as.factor(rep(rownames(reference$quali.var$coord), res$call$nboot)),
      oo$ind.sup$coord[, 1:2]), level.conf = level.conf)
    if (graph.type=="classic"){
	  plot.tmp(oo, axes = axes, col.ind.sup = rep(1:nrow(reference$quali.var$coord), res$call$nboot),
          label = "ind", col.quali = 1, ellipse = el,title = title, invisible = c("ind.sup","var"), new.plot = FALSE, graph.type=graph.type)
    } else {
     graph$PlotModSupp <- plot.tmp(oo, axes = axes, col.ind.sup = rep(1:nrow(reference$quali.var$coord), res$call$nboot),
          label = "ind", col.quali = 1, ellipse = el,title = title, invisible = c("ind.sup","var"), new.plot = FALSE,graph.type=graph.type)
    }
  }
  
  if (("all"%in%choice)|("quanti"%in%choice)){
  if ((new.plot)&!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) dev.new()
  color = c("black", "red", "green3", "blue", "cyan", "magenta", 
            "darkgray", "darkgoldenrod", "darkgreen", "violet", 
            "turquoise", "orange", "lightpink", "lavender", "yellow", 
            "lightgreen", "lightgrey", "lightblue", "darkkhaki", 
            "darkmagenta", "darkolivegreen", "lightcyan", "darkorange", 
            "darkorchid", "darkred", "darksalmon", "darkseagreen", 
            "darkslateblue", "darkslategray", "darkslategrey", 
            "darkturquoise", "darkviolet", "lightgray", "lightsalmon", 
            "lightyellow", "maroon")
  res.var <- res$res.imputeFAMD$completeObs
  tabdisj.var <- res$res.imputeFAMD$tab.disj
  for (i in 1:length(res$res.MI)){
    rec.pca <- res$res.MI[[i]][,sapply(res$call$X,is.numeric)]
    res.var <- cbind.data.frame(res.var,rec.pca)
    tabdisj.var <- cbind(tabdisj.var,rec.pca)
  }

  tabdisj.var <- as.matrix(tabdisj.var)
  colnames(res.var)=paste("V",1:ncol(res.var))
  colnames(res.var)[1:ncol(res$call$X)]=colnames(res$call$X)
  colnames(tabdisj.var)=paste("V",1:ncol(tabdisj.var))
  colnames(tabdisj.var)[1:sum(sapply(res$call$X,is.numeric))]=colnames(res$call$X)[sapply(res$call$X,is.numeric)]
  oo=FactoMineR::FAMD(res.var,sup.var=c((ncol(res$call$X)+1):ncol(res.var)),graph=FALSE,tab.disj=tabdisj.var)
  if (is.null(main)) title="Quantitative variables representation"    
  PlotVar <- plot(oo, axes=axes, choix = "quanti", title=title,invisible = "quanti.sup", col.hab = color[1:sum(sapply(res$call$X,is.numeric))],new.plot=FALSE,graph.type=graph.type)
  if (graph.type=="classic") {
    for (k in 1:res$call$nboot) points(oo$quanti.sup$coord[((k-1)*sum(sapply(res$call$X,is.numeric))+1):(k*sum(sapply(res$call$X,is.numeric))),axes[1]], oo$quanti.sup$coord[((k-1)*sum(sapply(res$call$X,is.numeric))+1):(k*sum(sapply(res$call$X,is.numeric))),axes[2]], col = color[1:sum(sapply(res$call$X,is.numeric))], pch = 15, cex = 0.3)
  } else {
    y <- variables <- NULL ## to avoid no visible binding for global variable
    dta <- cbind.data.frame(x=oo$quanti.sup$coord[,axes[1]],y=oo$quanti.sup$coord[,axes[2]],variables=rep(colnames(res$call$X[,sapply(res$call$X,is.numeric)]),res$call$nboot))
    PlotVar <- PlotVar + geom_point(data=dta,aes(x=x,y=y,color=variables),alpha=0.5,shape=20,size=2)+guides(color = guide_legend(title=NULL,override.aes = list(alpha=1, size=4)))
    print(PlotVar)
    graph$PlotVar <- PlotVar
  }
}

   if (graph.type=="ggplot") return(graph)
}

