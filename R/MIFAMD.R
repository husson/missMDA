MIFAMD <-
function(X,
                 ncp = 2,
                 method = c("Regularized","EM"),
                 coeff.ridge = 1,
                 threshold = 1e-06,
                 seed = NULL,
                 maxiter = 1000,
                 nboot=20,
                 verbose=T
){
  
  #intern functions
  
  estim.sigma2<-function(Xquanti,Xquali,M,Zhat,ncp,WW,D){
    # tab.disjonctif.NA <- function(tab) {
      # if(ncol(tab)==0){return(NULL)}
      # tab <- as.data.frame(tab)
      # modalite.disjonctif <- function(i) {
        # moda <- tab[, i]
        # nom <- names(tab)[i]
        # n <- length(moda)
        # moda <- as.factor(moda)
        # x <- matrix(0, n, length(levels(moda)))
        # ind <- (1:n) + n * (unclass(moda) - 1)
        # indNA <- which(is.na(ind))
        # x[(1:n) + n * (unclass(moda) - 1)] <- 1
        # x[indNA, ] <- NA
        # if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), 
                                                       # "n", "N", "y", "Y"))) 
          # dimnames(x) <- list(row.names(tab), paste(nom, 
                                                    # levels(moda), sep = "."))
        # else dimnames(x) <- list(row.names(tab), levels(moda))
        # return(x)
      # }
      # if (ncol(tab) == 1) 
        # res <- modalite.disjonctif(1)
      # else {
        # res <- lapply(1:ncol(tab), modalite.disjonctif)
        # res <- as.matrix(data.frame(res, check.names = FALSE))
      # }
      # return(res)
    # }
    
    reconst.FAMD<-function(xxquanti,xxquali,M=NULL,D=NULL,ncp,coeff.ridge=1,method="em"){
      zz<-cbind.data.frame(xxquanti,xxquali)
      if(is.null(M)){M<-c(1/apply(xxquanti,2,var),colMeans(zz)[-c(1:ncol(xxquanti))])}
      if(is.null(D)){D<-rep(1/nrow(zz),nrow(zz))}
      moy<-colMeans(zz)
      zzimp<-sweep(zz,MARGIN = 2,FUN = "-",STATS = moy)
      res.svd<-svd.triplet(zzimp,col.w = M,row.w = D,ncp=ncp)
      tmp<-seq(ncol(zz)-ncol(xxquali))
      if (nrow(zz) > length(tmp)){ 
        moyeig <- mean(res.svd$vs[tmp[-seq(ncp)]]^2)
      }else{
        moyeig <- mean(res.svd$vs[-c(1:ncp)]^2)
      }
      moyeig <- min(moyeig * coeff.ridge, res.svd$vs[ncp +1]^2)
      moyeigret<-moyeig
      if (method == "em"){
        moyeig <- 0
      }
      if(ncp>1){
        eig.shrunk <- (res.svd$vs[1:ncp]^2 - moyeig)/res.svd$vs[1:ncp]
      }else if(ncp==1){
        eig.shrunk <- matrix((res.svd$vs[1:ncp]^2 - moyeig)/res.svd$vs[1:ncp],1,1)
      }
      zzhat<-tcrossprod(res.svd$U%*%diag(eig.shrunk),res.svd$V[which(apply(is.finite(res.svd$V),1,any)),,drop=FALSE])
      zzhat<-sweep(zzhat,MARGIN = 2,FUN = "+",STATS = moy)
      return(list(zzhat=zzhat,moyeig=moyeigret,res.svd=res.svd,M=M))
    }
    
    
    nb.obs<-sum(WW[,-cumsum(sapply(Xquali,nlevels))])
    nb.obs.quanti<-sum(WW[,seq(ncol(Xquanti))])
    Zhat2<-reconst.FAMD(Zhat[,seq(ncol(Xquanti))],Zhat[,-seq(ncol(Xquanti))],ncp=ncp,D = D,M=M)
    Residu<-(Zhat-Zhat2$zzhat)%*%diag(M)^{1/2}
    Residu[is.na(cbind.data.frame(Xquanti,tab.disjonctif(Xquali)))]<-0
    sigma2<-sum(WW[,seq(ncol(Xquanti))]*((Residu[,seq(ncol(Xquanti))])^2))/(nb.obs.quanti-(ncol(Xquanti)+ncp*(sum(D)-1)+ncol(Xquanti)-ncp))
    return(sigma2)
  }
  
  
  imputeFAMD.stoch<-function(don, 
                             ncp = 4, 
                             method = c("Regularized", "EM"), 
                             row.w = NULL, 
                             coeff.ridge = 1, 
                             threshold = 1e-06, 
                             seed = NULL,
                             maxiter = 1000,
                             nboot,
                             verbose=FALSE){
    normtdc <- function(tab.disj, data.na) {
      tdc <- tab.disj
      tdc[tdc < 0] <- 0
      tdc[tdc > 1] <- 1
      col.suppr <- cumsum(sapply(data.na, function(x) {nlevels(x)}))
      tdc <- t(apply(tdc, 1, FUN = function(x, col.suppr) {
        if (sum(x[1:col.suppr[1]]) != 1) {
          x[1:col.suppr[1]] <- x[1:col.suppr[1]]/sum(x[1:col.suppr[1]])
        }
        if(length(col.suppr)>1){
          for (i in 2:length(col.suppr)) {
            x[(col.suppr[i - 1] + 1):(col.suppr[i])] <- x[(col.suppr[i - 
                                                                       1] + 1):(col.suppr[i])]/sum(x[(col.suppr[i - 
                                                                                                                  1] + 1):col.suppr[i]])
          }}
        return(x)
      }, col.suppr = col.suppr))
      return(tdc)
    }
    draw <- function(tabdisj, Don) {
      nbdummy <- rep(1, ncol(Don))
      is.quali <- which(!sapply(Don, is.numeric))
      nbdummy[is.quali] <- sapply(Don[, is.quali, drop = FALSE], nlevels)
      vec = c(0, cumsum(nbdummy))
      Donres <- Don
      for (i in is.quali) {
        Donres[, i] <- as.factor(levels(Don[, i])[apply(tabdisj[, 
                                                                (vec[i] + 1):vec[i + 1]], 1, function(x) {
                                                                  sample(1:length(x), size = 1, prob = x)
                                                                })])
        Donres[, i] <- factor(Donres[, i], levels(Don[, is.quali,drop=FALSE][, 
                                                                             i]))
      }
      return(don.imp = Donres)
    }
    
    quanti<-which(sapply(don,is.numeric))
    Ncol<-length(quanti)+sum(sapply(don[,-quanti],nlevels))
    W<-matrix(0,nrow(don),Ncol);W[!is.na(don)]<-1
    
    if(is.null(row.w)){
      Row.w<-as.integer(rep(1,nrow(don)))
    }else{
      Row.w<-row.w
    }
    
    if(is.integer(Row.w)){
      D<-Row.w/sum(Row.w)
      D[which(Row.w==0)]<-1/(1000*nrow(don))#D without 0
      WW<-diag(Row.w)%*%W
    }else{
      stop(paste0("row.w of class ",class(Row.w),", it needs to be an integer"))
    }
    
    
    res.imp<-imputeFAMD(X=don, ncp = ncp, method = method, row.w = D,
                        coeff.ridge = coeff.ridge, threshold = threshold, seed = seed, maxiter = maxiter)
    
    
    var_homo<-estim.sigma2(Xquanti=don[,quanti,drop=FALSE],
                           Xquali=don[,-quanti,drop=FALSE],
                           M=1/apply(res.imp$tab.disj,2,var),
                           Zhat=res.imp$tab.disj,
                           ncp=ncp,
                           D=D,WW=WW)
    sigma2<-var_homo/(1/apply(res.imp$tab.disj,2,var)[quanti])
    res.imp$fittedX<-res.imp$tab.disj
    res.imp$quanti.act<-which(sapply(don,is.numeric))
    
    
    
    classvar<-unlist(lapply(lapply(don,class),"[",1))#quand le type est "ordered" il y a 2 classes pour la variable
    if("integer"%in%classvar){classvar[classvar=="integer"]<-"numeric"}
    if("ordered"%in%classvar){classvar[classvar=="ordered"]<-"factor"}
    
    donimp<-don
    donimp[,which(classvar=="numeric")]<-res.imp$tab.disj[,seq(length(res.imp$quanti.act))]#les quanti active sont les premieres variables du tdc
    missing.quanti <-is.na(don[,res.imp$quanti.act])
    res.MI<-vector("list",length=nboot);names(res.MI)<-paste("nboot=",1:nboot,sep="")
    for(i in seq(nboot)){ 
      if(verbose){cat(paste(i, "...", sep = ""))}
      donimp.tmp<-donimp
      if(any("factor"%in% classvar)){
        tdc.imp <- res.imp$tab.disj[,(length(c(res.imp$quanti.act,res.imp$quanti.sup))+1):(ncol(res.imp$tab.disj)-length(res.imp$quali.sup)),drop=FALSE]
        tdc.norm <- normtdc(tab.disj = tdc.imp, data.na = don[,which(classvar=="factor"),drop=FALSE])
        donimp.quali<-draw(tdc.norm,don[,which(classvar=="factor"),drop=FALSE])
        donimp.tmp[,which(classvar=="factor")]<-donimp.quali[,names(which(classvar=="factor"))]
      }
      donimp.tmp[,res.imp$quanti.act][missing.quanti]<-res.imp$fittedX[,res.imp$quanti.act][missing.quanti]+rmvnorm(nrow(don),sigma=diag(sigma2))[missing.quanti]
      res.MI[[paste("nboot=",i,sep="")]]<-donimp.tmp
    }
    if (verbose) {
      cat("\ndone!\n")
    }
    res.MI<-list(res.MI=res.MI,sigma2=sigma2)
    class(res.MI)<-"MIFAMD"
    return(res.MI)
  }
  
  imputeFAMD.stoch.print <- function(don, ncp, method = c("Regularized", 
                                                          "EM"), row.w = NULL, coeff.ridge = 1, threshold = 1e-06, 
                                     seed = NULL, maxiter = 1000, verbose, printm) {
    if (verbose) {
      cat(paste(printm, "...", sep = ""))
    }
    res <- imputeFAMD.stoch(don = don, ncp = ncp, method = method, 
                            row.w = row.w, coeff.ridge = coeff.ridge, threshold = threshold, 
                            seed = seed, maxiter = maxiter,nboot=1)$res.MI[[1]]
    return(res)
  }
  
  #check if data are mixed
  if(sum(sapply(X,is.numeric))==ncol(X)){stop("All variables are numeric, use MIPCA")
  }else if(sum(sapply(X,is.numeric))==0){stop("No variable is numeric, use MIMCA")}
  
  #variables are ordered
  don<-X[,c(which(sapply(X,is.numeric)),which(!sapply(X,is.numeric)))]
  
  #print 
  temp <- if (coeff.ridge == 1) {
    "regularized"
  }
  else if ((coeff.ridge == 0) |(method=="EM")) {
    "EM"
  }else {
    paste("coeff.ridge=", coeff.ridge)
  }
  
  if (verbose) {
    cat("Multiple Imputation using", temp, "FAMD using", nboot, 
        "imputed arrays", "\n")
  }
  
  #multiple imputation
  n <- nrow(don)
  Boot <- matrix(sample(1:n, size = nboot * n, replace = T), n, nboot)
  Boot<-lapply(as.data.frame(Boot),FUN=function(xx){
    yy<-as.factor(xx)
    levels(yy)<-c(levels(yy),xx[which(!xx%in%as.numeric(as.character(levels(yy))))])
    return(yy)})
  Weight <- as.data.frame(matrix(0, n, nboot, dimnames = list(1:n,paste("nboot=", 1:nboot, sep = ""))))
  Boot.table <- lapply(Boot, table)
  for (i in 1:nboot) {
    Weight[names(Boot.table[[i]]), i] <- Boot.table[[i]]
  }
  Weight <-do.call(cbind.data.frame,lapply(Weight,as.integer))
  res.imp <- mapply(Weight,
                    FUN = imputeFAMD.stoch.print,
                    MoreArgs = list(don = don,ncp = ncp,
                                    coeff.ridge = coeff.ridge,
                                    method =method,
                                    threshold = threshold,
                                    maxiter = maxiter,
                                    verbose=verbose,
                                    seed=NULL),printm = as.character(1:nboot),
                    SIMPLIFY = FALSE)
  res <- list(res.MI = res.imp)
  res <- list(res.MI = lapply(res$res.MI,function(xx,nom){return(xx[,nom])},nom=colnames(X)),
              res.imputeFAMD = imputeFAMD(X,ncp = ncp, coeff.ridge = coeff.ridge, method =method,  threshold = threshold,  maxiter = maxiter,seed=seed),
              call=list(X = X,nboot = nboot, ncp = ncp, coeff.ridge = coeff.ridge, threshold = threshold, seed = seed, maxiter = maxiter))
  class(res) <- c("MIFAMD", "list")
  if (verbose) {
    cat("\ndone!\n")
  }
  return(res)
}
