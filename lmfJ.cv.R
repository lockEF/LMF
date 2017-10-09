lmfJ.cv <- function (x,y,z,mindex,maxiter=10) {
  xmiss <- x
  xmiss[mindex] <- NA

  # 2) Loop for multiple rank choices
  r <- 0
  # 2a) LMF Imputation Algorithm for Rank 0
  stemp <- lmf.impute(xmiss, y, z, r)
  save(stemp, file=paste("lmfJ",r,".R",sep="_"))
  cveC <- vnorm(x[mindex] - stemp[mindex])
  for (i in 1:maxiter) {
    # 2b) LMF Imputation Algorithm at Each r+1
    stemp <- lmf.impute(xmiss, y, z, r+1)
    save(stemp, file=paste("lmfJ",r+1,".R",sep="_"))
    cveJ <- vnorm(x[mindex] - stemp[mindex])

    cv.vec <- c(cveC,cveJ)
    if(which.min(cv.vec)==1) { break }
    if(which.min(cv.vec)==2) { r <- r+1; cveC <- cveJ }
  }

  return(r)
}

svd.cv <- function (x,mindex,maxiter=10) {
  xmiss <- x
  xmiss[mindex] <- NA

  # 2) Loop for multiple rank choices
  r <- 0
  # 2a) LMF Imputation Algorithm for Rank 0
  stemp <- svd.impute(xmiss, r)
  save(stemp, file=paste("svd",r,".R",sep="_"))
  cveC <- vnorm(x[mindex] - stemp[mindex])
  for (i in 1:maxiter) {
    # 2b) LMF Imputation Algorithm at Each r+1
    stemp <- svd.impute(xmiss, r+1)
    save(stemp, file=paste("svd",r+1,".R",sep="_"))
    cveJ <- vnorm(x[mindex] - stemp[mindex])

    cv.vec <- c(cveC,cveJ)
    if(which.min(cv.vec)==1) { break }
    if(which.min(cv.vec)==2) { r <- r+1; cveC <- cveJ }
  }

  return(r)
}


