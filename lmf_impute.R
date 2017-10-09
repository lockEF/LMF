## SVD Imputation Approaches
# Initialize using column or row averages
# Compute SVD
# Regress missing values on SVD components

lmf.impute <- function (x,y,z,r) {
  # 1) Initialize using row and column averages
  mindex <- is.na(x)
  mcol <- colMeans(x, na.rm=T)
  mrow <- rowMeans(x, na.rm=T)
  mtot <- mean(x, na.rm=T)
  last <- 0
  xtemp <- x
  for (i in 1:nrow(x)) {
   for (j in 1:ncol(x)) {
    # For single missing values, average col and row means
    xtemp[i,j] <- (mrow[i] + mcol[j]) / 2
   }
  }
  for (i in 1:nrow(x)) {
    if (is.na(mrow[i])) {
       xtemp[i,] <- mcol
    }
  }
  for (j in 1:ncol(x)) {
    if (is.na(mcol[j])) {
       xtemp[,j] <- mrow
    }
  }
  for (i in 1:nrow(x)) {
   for (j in 1:ncol(x)) {
    if (is.na(mcol[j]) & is.na(mrow[i])) {
       xtemp[i,j] <- mtot
    }
   }
  }
  xtemp[!mindex] <- x[!mindex]
  ## Replace missing values with average of row and column mean
  #rcmeans <- (matrix(rep(mcol,nrow(x)), nrow(x), ncol(x)) + 
  #            matrix(rep(mcol,ncol(x)), nrow(x), ncol(x))) / 2
  #xtemp[mindex] <- rcmeans[mindex]
  #colmiss <- is.na(mcol)
  #rowmiss <- is.na{mrow}
  

 conv <- F
 while(!conv) {
  # 2) Compute rank r Joint LMF
  stemp <- lmfJ(xtemp, y, z, r, scale=F)
  
  # 3) Impute missing values using LMF in (2) (EM algorithm)
  last <- xtemp
  xfull <- stemp$Jx
  xtemp[mindex] <- xfull[mindex]

  # Check for convergence
  if (norm(xtemp - last, type='f')^2 / norm(xtemp, type='f')^2 < .0001) { conv <- T }
 }
  return(xtemp)
}

# Try first coding these up for a single matrix (SVD)
svd.impute <- function (x, r) {
  # 1) Initialize using column averages
  mindex <- is.na(x)
  mcol <- colMeans(x, na.rm=T)
  mrow <- rowMeans(x, na.rm=T)
  mtot <- mean(x, na.rm=T)
  last <- 0
  xtemp <- x
  for (i in 1:nrow(x)) {
   for (j in 1:ncol(x)) {
    # For single missing values, average col and row means
    xtemp[i,j] <- (mrow[i] + mcol[j]) / 2
   }
  }
  for (i in 1:nrow(x)) {
    if (is.na(mrow[i])) {
       xtemp[i,] <- mcol
    }
  }
  for (j in 1:ncol(x)) {
    if (is.na(mcol[j])) {
       xtemp[,j] <- mrow
    }
  }
  for (i in 1:nrow(x)) {
   for (j in 1:ncol(x)) {
    if (is.na(mcol[j]) & is.na(mrow[i])) {
       xtemp[i,j] <- mtot
    }
   }
  }
  xtemp[!mindex] <- x[!mindex]

### LOOP 2 and 3 until convergence
 conv <- F
 while(!conv) {
  # 2) Compute rank r SVD
  stemp <- svd(xtemp, nu=r, nv=r)
  
  # 3) Impute missing values using SVD in (2) (EM algorithm)
  if (r>0) {
    xfull <- stemp$u[,1:r] %*% diag(x=stemp$d[1:r], nrow=r) %*% t(stemp$v[,1:r])
  } else if (r==0) {
    xfull <- matrix(0,nrow(x),ncol(x))
  }
  last <- xtemp
  xtemp[mindex] <- xfull[mindex]

  # Check for convergence
  if (norm(xtemp - last, type='f')^2  / norm(xtemp, type='f')^2 < .0001) { conv <- T }
 }
  return(xtemp)
}

# SVD only (no regression iteration)
svd.only.impute <- function (x, r) {
  # 1) Initialize using column averages
  mindex <- is.na(x)
  mcol <- colMeans(x, na.rm=T)
  last <- 0
  xtemp <- x
  for (i in 1:ncol(x)) {
    xtemp[which(mindex[,i]),i] <- mcol[i]
  }

  # 2) Compute rank r SVD
  stemp <- svd(xtemp, nu=r, nv=r)
  
  # 3) Impute missing values using SVD in (2) (EM algorithm)
  if (r>0) {
    xfull <- stemp$u[,1:r] %*% diag(x=stemp$d[1:r], nrow=r) %*% t(stemp$v[,1:r])
  } else if (r==0) {
    xfull <- matrix(0,nrow(x),ncol(x))
  }
  xtemp[mindex] <- xfull[mindex]

  return(xtemp)
}











