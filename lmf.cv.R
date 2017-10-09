lmf.cv <- function (x,y,z,maxiter=10) {
  # 1) Generate missing values (~18% of x set to missing)
  nmissR <- 20 # ~5% of missing rows
  nmissC <- 5  # ~5% of missing cols
  nmissI <- 4555  # total missing vals from R/C
  # Set a few rows to missing
  xmiss <- x
  rmiss <- sample(1:nrow(x), nmissR, replace=F)
  cmiss <- sample(1:ncol(x), nmissC, replace=F)
  xmiss[rmiss,] <- NA
  xmiss[,cmiss] <- NA
  # Set some individual entires to missing
  rmissI <- sample(1:nrow(x), nmissI, replace=T)
  cmissI <- sample(1:ncol(x), nmissI, replace=T)
  for (j in 1:nmissI) {
    xmiss[rmissI[j], cmissI[j]] <- NA
  }
  mindex <- is.na(xmiss)

  # 2) Loop for multiple rank choices
  r <- 1; rx <- 2; ry <- 0; rz <- 1
  # 2a) LMF Imputation Algorithm for Rank 0
  stemp <- lmf.impute.ji(xmiss, y, z, r, rx, ry, rz)
  save(stemp, file=paste("lmf",r,rx,ry,rz,".R",sep="_"))
  cveC <- vnorm(x[mindex] - stemp$x[mindex])
  for (i in 1:maxiter) {
    # 2b) LMF Imputation Algorithm at Each r+1
    stemp <- lmf.impute.ji(xmiss, y, z, r+1, rx, ry, rz)
    save(stemp, file=paste("lmf",r+1,rx,ry,rz,".R",sep="_"))
    cveJ <- vnorm(x[mindex] - stemp$x[mindex])
    stemp <- lmf.impute.ji(xmiss, y, z, r, rx+1, ry, rz)
    save(stemp, file=paste("lmf",r,rx+1,ry,rz,".R",sep="_"))
    cveX <- vnorm(x[mindex] - stemp$x[mindex])
    stemp <- lmf.impute.ji(xmiss, y, z, r, rx, ry+1, rz)
    save(stemp, file=paste("lmf",r,rx,ry+1,rz,".R",sep="_"))
    cveY <- vnorm(x[mindex] - stemp$x[mindex])
    stemp <- lmf.impute.ji(xmiss, y, z, r, rx, ry, rz+1)
    save(stemp, file=paste("lmf",r,rx,ry,rz+1,".R",sep="_"))
    cveZ <- vnorm(x[mindex] - stemp$x[mindex])

    cv.vec <- c(cveC,cveJ,cveX,cveY,cveZ)
    if(which.min(cv.vec)==1) { break }
    if(which.min(cv.vec)==2) { r <- r+1; cveC <- cveJ }
    if(which.min(cv.vec)==3) { rx <- rx+1; cveC <- cveX }
    if(which.min(cv.vec)==4) { ry <- ry+1; cveC <- cveY }
    if(which.min(cv.vec)==5) { rz <- rz+1; cveC <- cveZ }
  }

  return(list(c(r,rx,ry,rz),xmiss,rmiss,cmiss))
}

# Function to return a vector norm
vnorm <- function(x) {
  return(sqrt(sum(x^2)))
}






