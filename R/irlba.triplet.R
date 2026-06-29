irlba.triplet <- function (X, row.w = NULL, col.w = NULL,ncp=5) {

   if (is.null(row.w)) row.w <- rep(1/nrow(X), nrow(X))
   if (is.null(col.w)) col.w <- rep(1, ncol(X))
   ncp <- min(ncp,nrow(X)-1,ncol(X))
   row.w <- row.w / sum(row.w)
   X <- t(t(X)*sqrt(col.w))*sqrt(row.w)
   if (ncp>=0.5*min(length(row.w),length(col.w))) svd.usuelle <- svd(X,nu=ncp,nv=ncp)
   else svd.usuelle <- irlba(X,nu=ncp,nv=ncp)
   U <- svd.usuelle$u/sqrt(row.w)
   V <- svd.usuelle$v/sqrt(col.w)
   vs <- svd.usuelle$d[1:min(ncol(X),nrow(X)-1,ncp)]
   num <- which(vs[1:ncp]<1e-15)
   if (length(num)==1){
	 U[,num] <- U[,num,drop=FALSE]*vs[num]
     V[,num] <- V[,num,drop=FALSE]*vs[num]
   } 
   if (length(num)>1){
	 U[,num] <- t(t(U[,num])*vs[num])
     V[,num] <- t(t(V[,num])*vs[num])
   }
   res <- list(vs = vs, U = U, V = V,sumvp=sum(X^2))
   return(res)
}