imputeMultilevel <- function(X, ifac=1, ncpB = 2, ncpW = 2, method  = c("Regularized","EM"), 
      scale=TRUE, row.w=NULL, threshold=1e-4, maxiter = 1000, ...){
					 
  find.category <- function (X,tabdisj){
  # X matrix of categoriccal variables
  # tabdisj fuzzy disjunctive table of the categorical variables
  nbdummy <- NULL
  if (ncol(X)>1){
     for (i in 1:ncol(X)) nbdummy <- c(nbdummy,nlevels(X[,i]))
  } else {
    nbdummy <- nlevels(X[,1])
  }
#  if (ncol(X)>1) nbdummy <- unlist(lapply(X,nlevels))
#  else nbdummy <- nlevels(X[,1,drop=FALSE])
  vec = c(0,cumsum(nbdummy))
  for (i in 1:ncol(X)) {
    temp <- as.factor(levels(X[, i])[apply(tabdisj[,(vec[i] + 1):vec[i + 1]], 1, which.max)])
    X[,i]<-factor(temp,levels(X[,i]))
#    temp <- apply(tabdisj[,(vec[i] + 1):vec[i + 1]], 1, which.max)
#    X[,i]<-factor(as.character(temp),labels=levels(X[,i]))
  }  
  return(X)
}

# tab.disjonctif.NA <- function(tab) {
  # tab <- as.data.frame(tab)
  # modalite.disjonctif <- function(i) {
    # moda <- tab[, i]
    # nom <- names(tab)[i]
    # n <- length(moda)
    # moda <- as.factor(moda)
    # x <- matrix(0L, n, length(levels(moda)))
    # ind <- (1:n) + n * (unclass(moda) - 1L)
    # indNA <- which(is.na(ind))
    # x[(1:n) + n * (unclass(moda) - 1)] <- 1L
    # x[indNA, ] <- NA
    # dimnames(x) <- list(row.names(tab), levels(moda))
    # return(x)
  # }
  # if (ncol(tab) == 1) res <- modalite.disjonctif(1)
  # else {
    # res <- lapply(1:ncol(tab), modalite.disjonctif)
    # res <- as.matrix(data.frame(res, check.names = FALSE))
  # }
  # return(res)
# }

moy.p <- function(V, poids) {
  res <- sum(V * poids,na.rm=TRUE)/sum(poids[!is.na(V)])
}
ec <- function(V, poids) {
  res <- sqrt(sum(V^2 * poids,na.rm=TRUE)/sum(poids[!is.na(V)]))
}
					 
X <- as.data.frame(X)
method <- match.arg(method,c("Regularized","regularized","EM","em"),several.ok=T)[1]
  method=tolower(method)	
  if (is.null(row.w)) {row.w = rep(1, nrow(X))/nrow(X)} # Weights useful for bootstrap. If no bootstrap they are set to 1/n.
  if (any(is.na(X[,ifac]))){
    warning("Rows with missing valued on the group variables are dropped")
	row.w = row.w[!is.na(X[,ifac])] # Remove rows where group membership is missing
	X = X[!is.na(X[,ifac]),] # Remove rows where group membership is missing
  }
  if(is.null(rownames(X)))  rownames(X)=1:nrow(X)
  indexrow = rownames(X)
  if(is.null(colnames(X)))  colnames(X)=1:ncol(X)
  indexcol = colnames(X)
  ## reorganize the variables
  is.quali <- which(!unlist(lapply(X,is.numeric)))
  nb.quali <- length(is.quali)
  is.quanti <- which(unlist(lapply(X,is.numeric)))
  nb.quanti <- length(is.quanti)
  X <- X[,c(ifac,is.quanti,setdiff(is.quali,ifac))]
  if (nb.quali>1){
    X[,(2+nb.quanti):ncol(X)] <- do.call(cbind.data.frame, lapply((2+nb.quanti):ncol(X), function(i) as.factor(X[,i])))
    levels.number <- unlist(lapply((2+nb.quanti):ncol(X), function(i) nlevels(X[,i])))
    levels.number.dropped <- unlist(lapply((2+nb.quanti):ncol(X), function(i) nlevels(droplevels(X)[,i])))
    if(! all(levels.number == levels.number.dropped)){
      X = droplevels(X)
      warning("Empty levels dropped")
    }
    liste.levels <- lapply(X,levels)
    for (j in (2+nb.quanti):ncol(X)) levels(X[,j])<- sort(paste("v",j,".",1:nlevels(X[,j]),sep=""))  # In case two cat variables have the same levels
  }
  if(ncpB > nlevels(X[, 1])-1){
    ncpB= nlevels(X[, 1])-1
    warning("Number of between components larger than number of groups. By default it was set to the number of groups minus 1.  ")
  }
  
  X = X[ order(X[,1]),] # Order data frame by group
  fac = X[, 1]
QualiAct <- QuantiAct <- NULL
  if (nb.quanti>0) QuantiAct <- as.matrix(X[,2:(1+nb.quanti),drop=FALSE]) # Data table containing all active numerical variables
  if (nb.quali>1){
    QualiAct <- as.matrix(X[,(2+nb.quanti):ncol(X),drop=FALSE]) # Table containing categorical variables
    Z <- Zold <- FactoMineR::tab.disjonctif.prop(QualiAct) # Fuzzy disjunctive table where missing values are replaced by proportions of categories
    colnames(Z) <- unlist(sapply(X[,(2+nb.quanti):ncol(X),drop=FALSE], function(cc) levels(droplevels(cc))))
  }
  
  if (nb.quali>1) {indNA = is.na(cbind(QuantiAct, tab.disjonctif(QualiAct)))} else {
    indNA = is.na(QuantiAct)} # Indicator matrix for missingness: missing->TRUE, observed->FALSE.
  
  Xhat <- NULL
  # Center, scale, initialize for numerical variables
  if (nb.quanti>0){
    mean.p <- apply(QuantiAct, 2, moy.p,row.w) # Compute weighted mean
    QuantiAct <- t(t(QuantiAct)-mean.p) # Center
    if (scale){
      et <- apply(QuantiAct, 2, ec,row.w)
      QuantiAct <- t(t(QuantiAct)/et) # Scale
    }
    QuantiAct[is.na(QuantiAct)]=0 # Initialize missing values by mean (0 because columns are centered)
    Xhat<-QuantiAct
  }
  
  # Center, scale, initialize for categorical variables
  if (nb.quali>1){
    Mglobal = apply(Z,2,mean) # nc/n
    Zcentered <- t(t(Z)-Mglobal)# Z - between  Z - nck/nk
    A = t(t(Zcentered)/sqrt(Mglobal)) # Zcentererd* (nck)^{-1/2}
    colnames(A) <- colnames(Z)
    Xhat<-cbind(Xhat,A)
  }
  
  Xhatold <- Xhat # Store previous value of Xhat
  continue = TRUE
  nb.iter <- 1
  critere = NULL
  while (continue){
    nb.iter = nb.iter+1
    if (nb.quanti>0){
      # Quanti: recompute standard deviation and mean
      if (scale) Xhat[,1:ncol(QuantiAct)]=t(t(Xhat[,1:ncol(QuantiAct)])*et)
      Xhat[,1:ncol(QuantiAct)] <- t(t(Xhat[,1:ncol(QuantiAct)])+mean.p)
      mean.p <- apply(Xhat[,1:ncol(QuantiAct)], 2, moy.p,row.w)
      Xhat[,1:ncol(QuantiAct)] <- t(t(Xhat[,1:ncol(QuantiAct)])-mean.p)
      if (scale){
        et <- apply(Xhat[,1:ncol(QuantiAct)], 2, ec,row.w)
        Xhat[,1:ncol(QuantiAct)] <- t(t(Xhat[,1:ncol(QuantiAct)])/et)
      }
    }
    
    if (nb.quali>1){
      # Quali: recompute margins
      Z <-  t(t(Xhat[, (nb.quanti+1):ncol(Xhat)])*sqrt(Mglobal) + Mglobal)
      Mglobal = apply(Z,2,mean) # nc/n
      Zcentered <- t(t(Z)-Mglobal)# Z - between  Z - nck/nk
#      A = Zcentered%*%diag(1/sqrt(Mglobal)) # Zcentered* (nck)^{-1/2}
      A = t(t(Zcentered)/sqrt(Mglobal)) # Zcentererd* (nck)^{-1/2}
      Xhat[,(nb.quanti+1):ncol(Xhat)] <-A
    }
    BETWEEN = aggregate(Xhat, list(fac), mean, na.rm=TRUE)
	rownames(BETWEEN) <- levels(fac)
    BETWEEN = as.matrix(BETWEEN[, -1])
    
    WITHIN_global = NULL
    for (k in 1:nlevels(fac)){
      #print(k)
      moy = c(as.matrix(BETWEEN[k,])) # Group means
      inter= t(t(Xhat[ fac==levels(fac)[k], , drop = FALSE]) - moy) # Center by group
      WITHIN_global = rbind(WITHIN_global,inter) 
    }
    
    svd.WITHIN_global <- FactoMineR::svd.triplet(WITHIN_global,ncp=ncpW) # SVD of the Within part
    moyeigG <- mean(svd.WITHIN_global$vs[-c(1:ncpW)]^2) # Check when n<p
    if (method=="em") moyeigG <- 0
    eig.shrunkG <- ((svd.WITHIN_global$vs[1:ncpW]^2 - moyeigG)/svd.WITHIN_global$vs[1:ncpW])
    
    ReconWITHIN_global = t(t(svd.WITHIN_global$U[,1:ncpW, drop = FALSE])* eig.shrunkG)%*% t(svd.WITHIN_global$V[,1:ncpW, drop = FALSE])
    
    W = sqrt(table(fac))
    mat_diag_bet= sweep(BETWEEN,1,W,FUN="*") # Weight groups by square of effective
    ### Svd of BETEWEEN
    svd_b=svd(mat_diag_bet)#,nu=qb,nv=qb)
    Fb= sweep(svd_b$u[,1:ncpB,drop=F],1,1/W,FUN="*")
    if(length(svd_b$d)>ncpB) sigma2_b <- mean(svd_b$d[setdiff(1:length(svd_b$d),(1:ncpB))]^2) else sigma_b = 0  
    if (method=="em")  sigma2_b <- 0          
    lambda.shrinked_b = (svd_b$d[1:ncpB]^2 - sigma2_b)/(svd_b$d[1:ncpB])            
    Bb=t(t(svd_b$v[,1:ncpB,drop=FALSE])*lambda.shrinked_b)           
    Fsupb=NULL
    for (i in 1: nlevels(fac)) {
      Fsupb=rbind.data.frame(Fsupb,rep(1,length= table(fac)[i]) %*%t(Fb[i,]))
    }
    Fsupb=as.matrix(Fsupb)
    ReconBETWEEN  = Fsupb%*%t(Bb)  
    
    RECONGLOB =  ReconWITHIN_global + ReconBETWEEN
    
    Xhat[indNA] <- RECONGLOB[indNA]
    
    #print(sum((Xhatold[indNA]-Xhat[indNA])^2))

    if ((sum((Xhatold[indNA]-Xhat[indNA])^2)/length(indNA))<threshold || nb.iter>=maxiter) continue=FALSE
    Xhatold=Xhat
  }
  
  # Quanti recompute sds and means
  if (nb.quanti>0){
    if (scale) Xhat[,1:ncol(QuantiAct)]=t(t(Xhat[,1:ncol(QuantiAct)])*et)
    Xhat[,1:ncol(QuantiAct)] <- t(t(Xhat[,1:ncol(QuantiAct)])+mean.p)
  }
  # Quali recompute margins
  if (nb.quali>1){
    Xhat[, (nb.quanti+1):ncol(Xhat)] <- t(t(Xhat[, (nb.quanti+1):ncol(Xhat)])*sqrt(Mglobal) + Mglobal)
  }

  completeObs = X[, 1,drop=FALSE]
  if (nb.quanti>0) completeObs <- cbind.data.frame(completeObs,Xhat[,1:ncol(QuantiAct)])  
  if (nb.quali>1){
    completeObs <- cbind.data.frame(completeObs,find.category(X[,(2+nb.quanti):ncol(X),drop=FALSE], Xhat[, (nb.quanti+1):ncol(Xhat)]))
    colnames(completeObs)[1]=colnames(X)[1]
    for (j in (2+nb.quanti):ncol(X)) levels(completeObs[,j]) <-  liste.levels[[j]]
  }
  completeObs=completeObs[indexrow,indexcol]
  if(nb.iter>=maxiter){
    warning(paste("Stopped after", maxiter, "iterations"))
  }
  return(list(completeObs=completeObs,Xhat=Xhat))
}

