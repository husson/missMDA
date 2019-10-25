imputeFAMD <- function (X, ncp = 2, method=c("Regularized","EM"),row.w=NULL,coeff.ridge=1,threshold = 1e-6,
ind.sup=NULL, sup.var=NULL, seed = NULL,maxiter=1000,...){

  X <- as.data.frame(X)
  method <- match.arg(method,c("Regularized","regularized","EM","em"),several.ok=T)[1]
  method <- tolower(method)
  type = rep("s",ncol(X))
  type[!sapply(X,is.numeric)]="n"
  res <- imputeMFA(X=X,group=rep(1,ncol(X)),type=type, ncp = ncp, method=method,row.w=row.w,coeff.ridge=coeff.ridge,ind.sup=ind.sup, num.group.sup=sup.var,threshold = threshold,seed = seed,maxiter=maxiter)
  return(res)
}

