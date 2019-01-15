estim_ncpMCA <- function(don,ncp.min=0,ncp.max=5,method=c("Regularized","EM"),method.cv=c("Kfold","loo"),nbsim=100,pNA=0.05,threshold=1e-4,verbose=TRUE){

#### Debut tab.disjonctif.NA
tab.disjonctif.NA<-function (tab) {
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) {
        moda <- tab[, i]
        nom <- names(tab)[i]
        n <- length(moda)
        moda <- as.factor(moda)
        x <- matrix(0, n, length(levels(moda)))
	      ind<-(1:n) + n * (unclass(moda) - 1)
	      indNA<-which(is.na(ind))
		
        x[(1:n) + n * (unclass(moda) - 1)] <- 1
        x[indNA,]<-NA 
	  if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), "n", "N", "y", "Y"))) dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), sep = "."))
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
#### Fin tab.disjonctif.NA

  prodna<-function (x, noNA){
    n <- nrow(x)
    p <- ncol(x)
    NAloc <- rep(FALSE, n * p)
    NAloc[sample(n * p, floor(n * p * noNA))] <- TRUE
    x[matrix(NAloc, nrow = n, ncol = p)] <- NA
    return(x)
  }
  

########## Debut programme principal
don <- as.data.frame(don)
  method <- match.arg(method,c("Regularized","regularized","EM","em"),several.ok=T)[1]
method.cv <- match.arg(method.cv,c("loo","Kfold","kfold","LOO"),several.ok=T)[1]
method <- tolower(method)
method.cv <- tolower(method.cv)
auxi = NULL
don <- droplevels(don)
for (j in 1:ncol(don)) if (is.numeric(don[,j])) auxi = c(auxi,colnames(don)[j])
if (!is.null(auxi)) stop(paste("\nAll variables are not categorical, the following ones are numeric: ", auxi))
vrai.tab=tab.disjonctif.NA(don)

if (method.cv=="kfold"){
res = matrix(NA,ncp.max-ncp.min+1,nbsim)
if(verbose) pb <- txtProgressBar(min=1/nbsim*100, max=100,style=3)

for (sim in 1:nbsim){
 compteur<-1
 while(compteur<50){
    donNA <- prodna(don, pNA)
    for (i in 1:ncol(don)) donNA[,i]=as.factor(as.character(donNA[,i]))
    compteur <- 1+100*(sum(unlist(sapply(as.data.frame(donNA),nlevels)))==sum(unlist(sapply(don,nlevels))))
  }
  if (compteur<100) stop('It is too difficult to suppress some cells.\nMaybe several categories are taken by only 1 individual. You should uppress these variables or try with method.cv="loo".')
 for (nbaxes in ncp.min:ncp.max){
  tab.disj.comp <- imputeMCA(as.data.frame(donNA),ncp=nbaxes,method=method,threshold=threshold)$tab.disj
  if (sum(is.na(donNA))!=sum(is.na(don))) res[nbaxes-ncp.min+1,sim] <- sum((tab.disj.comp-vrai.tab)^2,na.rm=TRUE)/(sum(is.na(tab.disjonctif.NA(donNA)))-sum(is.na(tab.disjonctif.NA(don))))
 }
 if(verbose) setTxtProgressBar(pb, sim/nbsim*100)
}
if(verbose) close(pb)
crit=apply(res,1,mean,na.rm=TRUE)
 names(crit) <- c(ncp.min:ncp.max)
result = list(ncp = as.integer(which.min(crit)+ncp.min-1),criterion=crit)
return(result)
}

if (method.cv=="loo"){
  if(verbose) pb <- txtProgressBar(min = 0, max = 100, style = 3)
crit <- NULL
tab.disj.hat <- vrai.tab
col.in.indicator <- c(0,sapply(don,nlevels))
 for (nbaxes in ncp.min:ncp.max){
   for (i in 1:nrow(don)){
    for (j in 1:ncol(don)){
     if (!is.na(don[i,j])){
       donNA <- as.matrix(don)
       donNA[i,j] <- NA
       if(!any(unlist(sapply(as.data.frame(donNA),summary))==0)){
	     for (k in 1:ncol(donNA)) donNA[,k]=as.factor(as.character(donNA[,k]))
         tab.disj.hat[i,(cumsum(col.in.indicator)[j]+1):(cumsum(col.in.indicator)[j+1])] <- imputeMCA(as.data.frame(donNA),ncp=nbaxes,method=method,threshold=threshold)$tab.disj[i,(cumsum(col.in.indicator)[j]+1):(cumsum(col.in.indicator)[j+1])]
	   }
 }
}
    if(verbose) setTxtProgressBar(pb, round((((1:length(ncp.min:ncp.max))[which(nbaxes==(ncp.min:ncp.max))]-1)*nrow(don)+i)/(length(ncp.min:ncp.max)*nrow(don))*100))  
   }
crit <- c(crit,mean((tab.disj.hat-vrai.tab)^2,na.rm=TRUE))
}
  if(verbose) close(pb)
names(crit) <- c(ncp.min:ncp.max)
return(list(ncp = as.integer(which.min(crit)+ncp.min-1),criterion=crit))
}
}
