lmf.impute.ji <- function (x,y,z,r,rx,ry,rz,maxiter=500) {
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
  niter <- 0
  convlist <- c()
  while(!conv & niter < maxiter) {
    # 2) Compute rank r Joint LMF
    stemp <- lmf(xtemp, y, z, r, rx, ry, rz, scale=F)
    
    # 3) Impute missing values using LMF in (2) (EM algorithm)
    last <- xtemp
    xfull <- stemp$Jx + stemp$Ax
    xtemp[mindex] <- xfull[mindex]
    
    # Check for convergence
    convlist[niter+1] <- norm(xtemp - last, type='f')^2 / norm(xtemp, type='f')^2
    if (norm(xtemp - last, type='f')^2 / norm(xtemp, type='f')^2 < .0001) { conv <- T }
    print(paste("Iteration ", niter,": ",convlist[niter+1],sep=""))
    niter <- niter + 1
  }

  return(list(x=xtemp,nconv=niter,convlist=convlist))
}


# LMF with Joint and Individual Structure
# Scale convergence criterion by either data or previous estimate (leave tol for now)
lmf <- function(x, y, z, r, rx, ry, rz, tol=.00001, maxiter=500, scale=T, orth=0, anneal=F, nanneal=500, Tmax=0.1) {
  if (scale) {
    # Center and Scale data
    x <- x - matrix(rep(apply(x,1,mean,na.rm=T),ncol(x)),nrow=nrow(x))
    x <- x / norm(x, type="f")
    y <- y - matrix(rep(apply(y,1,mean,na.rm=T),ncol(y)),nrow=nrow(y))
    y <- y / norm(y, type="f")
    z <- z - matrix(rep(apply(z,1,mean,na.rm=T),ncol(z)),nrow=nrow(z))
    z <- z / norm(z, type="f")
  }
  
  # Initialize individual structure
  Ax <- matrix(0,nrow(x),ncol(x))
  Ay <- matrix(0,nrow(y),ncol(y))
  Az <- matrix(0,nrow(z),ncol(z))
  
  # Dimensions:
  m <- nrow(x) + nrow(y); n <- ncol(x) + ncol(z)
  mx <- nrow(x); nx <- ncol(x);
  if (r > 0) {
    # vtilde: starting value for vtilde
    ztilde <- cbind(x,z)
    vtilde <- svd(ztilde,nu=r,nv=r)$v
    v <- vtilde[1:nx,]
    sx <- diag(1,r)
    uy <- t(solve(t(v) %*% v)%*%t(v)%*%t(y))
  } else {
    sx <- NULL
  }
  
  current <- matrix(0,m,n)
  sse <- c()
  convs <- c()
  
  for (i in 1:maxiter) {
    if (r > 0) {
      ### Joint Structure
      x_j <- x - Ax
      y_j <- y - Ay
      z_j <- z - Az
      # Combined data sets
      xtilde <- rbind(x_j,y_j)
      ztilde <- cbind(x_j,z_j)
      xvec <- as.vector(x_j)
      
      u <- t(solve(t(vtilde) %*% vtilde)%*%t(vtilde)%*%t(ztilde))
      # Scale u
      for (j in 1:r) {
        u[,j] <- u[,j] / norm(matrix(u[,j],ncol=1), type="f")
      }
      utilde <- rbind(u %*% sx, uy)
      v <- t(solve(t(utilde) %*% utilde)%*%t(utilde)%*%xtilde)
      vz <- t(solve(t(u) %*% u)%*%t(u)%*%z_j)
      # Scale v
      for (j in 1:r) {
        v[,j] <- v[,j] / norm(matrix(v[,j],ncol=1), type="f")
      }
      #vtilde <- rbind(v %*% sx, vz)
      uy <- t(solve(t(v) %*% v)%*%t(v)%*%t(y_j))
      
      # Calculate scaling matrix Sx
      uv <- matrix(0,mx*nx,r)
      cnt <- 1
      for (k in 1:nx) {
        for (j in 1:mx) {
          uv[cnt,] <- u[j,] * v[k,]
          cnt <- cnt + 1
        }
      }
      sx_vec <- as.vector(solve(t(uv) %*% uv) %*% t(uv) %*% xvec)
      sx <- diag(sx_vec, r)
      # Re-estimate (to incorporate scaling)
      utilde <- rbind(u %*% sx, uy)
      vtilde <- rbind(v %*% sx, vz)
      
      # Set joint estimates
      Jx <- u %*% sx %*% t(v)
      Jy <- uy %*% t(v)
      Jz <- u %*% t(vz)
    } else {
      # Set joint estimates to 0 when r=0
      Jx <- matrix(0,nrow(x),ncol(x))
      Jy <- matrix(0,nrow(y),ncol(y))
      Jz <- matrix(0,nrow(z),ncol(z))
    }
    
    sse_current <- norm(x-Jx-Ax,type='f')^2 + norm(y-Jy-Ay,type='f')^2 + norm(z-Jz-Az,type='f')^2
    
    if (anneal) {
      ### Simulated annealing
      #Tmax <- .1
      if (i <= nanneal) {
        rx <- matrix(rnorm(nrow(Jx)*ncol(Jx),0,Tmax/i),nrow(Jx),ncol(Jx))
        ry <- matrix(rnorm(nrow(Jx)*ncol(Jx),0,Tmax/i),nrow(Jy),ncol(Jy))
        rz <- matrix(rnorm(nrow(Jx)*ncol(Jx),0,Tmax/i),nrow(Jz),ncol(Jz))
        
        sseA <- norm(x-Jx-rx-Ax,type='f')^2 + norm(y-Jy-ry-Ay,type='f')^2 + norm(z-Jz-rz-Az,type='f')^2
        p <- exp(-(sseA-sse_current)/(Tmax/i))
        #   if (sseA < sse_current | runif(1) < p) {
        temp <- svd(Jx + rx, nu=r, nv=r)
        Jx <- temp$u[,1:r] %*% t(temp$v[,1:r]) 
        temp <- svd(Jy + ry, nu=r, nv=r)
        Jy <- temp$u[,1:r] %*% t(temp$v[,1:r]) 
        temp <- svd(Jz + rz, nu=r, nv=r)
        Jz <- temp$u[,1:r] %*% t(temp$v[,1:r])
        
        #      print("accepted") 
        #   }
      }
    }
    
    sse <- c(sse,(norm(x-Jx-Ax,type='f')^2 + norm(y-Jy-Ay,type='f')^2 + norm(z-Jz-Az,type='f')^2))
    
    ### Individual Structure
    x_i <- (x - Jx)
    y_i <- (y - Jy)
    z_i <- (z - Jz)
    
    if (rx > 0) {
      s <- svd(x_i, nu=rx, nv=rx)
      Ax <- s$u[,1:rx] %*% diag(x=s$d[1:rx], nrow=rx) %*% t(s$v[,1:rx])
    }
    if (ry > 0) {
      s <- svd(y_i, nu=ry, nv=ry)
      Ay <- s$u[,1:ry] %*% diag(x=s$d[1:ry], nrow=ry) %*% t(s$v[,1:ry])
    }
    if (rz > 0) {
      s <- svd(z_i, nu=rz, nv=rz)
      Az <- s$u[,1:rz] %*% diag(x=s$d[1:rz], nrow=rz) %*% t(s$v[,1:rz])
    }
    
    ### Check convergence
    if (i > 1) { 
      conv <- (norm(Jx-last[[1]],type='f')^2 + norm(Jy-last[[2]],type='f')^2 + norm(Jz-last[[3]],type='f')^2 +
        norm(Ax-last[[4]],type='f')^2 + norm(Ay-last[[5]],type='f')^2 + norm(Az-last[[6]],type='f')^2) / 
        (norm(x,type='f')^2 + norm(y,type='f')^2 + norm(z,type='f')^2)
    } else {
      conv <- (norm(Jx,type='f')^2 + norm(Jy,type='f')^2 + norm(Jz,type='f')^2 +
        norm(Ax,type='f')^2 + norm(Ay,type='f')^2 + norm(Az,type='f')^2) / 
        (norm(x,type='f')^2 + norm(y,type='f')^2 + norm(z,type='f')^2)
    }
    convs <- c(convs,conv)
    last <- list(Jx,Jy,Jz,Ax,Ay,Az)
    sse <- c(sse,(norm(x-Jx-Ax,type='f')^2 + norm(y-Jy-Ay,type='f')^2 + norm(z-Jz-Az,type='f')^2))
    print(paste(conv,(norm(x-Jx-Ax,type='f')^2 + norm(y-Jy-Ay,type='f')^2 + norm(z-Jz-Az,type='f')^2),sep=" "))
    if (conv < tol) { break }
  }
  
  temp <- svd(Jx)
  u <- temp$u[,1:r]
  v <- temp$v[,1:r]
  
  Jy_perp = Jy + (Ay %*% v %*% t(v))
  Jz_perp = Jz + (u %*% t(u) %*% Az)
  Ay_perp = Ay %*% (diag(ncol(x)) - v %*% t(v))
  Az_perp = (diag(nrow(x)) - u %*% t(u)) %*% Az
  
  Jy=Jy_perp; Jz=Jz_perp;
  Ay=Ay_perp; Az=Az_perp;
  
  return(list(Jx=Jx,Jy=Jy,Jz=Jz,Ax=Ax,Ay=Ay,Az=Az,sx=sx,sse=sse,convs=convs))
}


# LMF with Individual Structure Estimated Before Joint
# Scale convergence criterion by either data or previous estimate (leave tol for now)
lmfR <- function(x, y, z, r, rx, ry, rz, tol=.00001, maxiter=500, scale=T, orth=0, anneal=F, nanneal=500, Tmax=0.1) {
  if (scale) {
    # Center and Scale data
    x <- x - matrix(rep(apply(x,1,mean,na.rm=T),ncol(x)),nrow=nrow(x))
    x <- x / norm(x, type="f")
    y <- y - matrix(rep(apply(y,1,mean,na.rm=T),ncol(y)),nrow=nrow(y))
    y <- y / norm(y, type="f")
    z <- z - matrix(rep(apply(z,1,mean,na.rm=T),ncol(z)),nrow=nrow(z))
    z <- z / norm(z, type="f")
  }
  
  # Initialize joint structure
  Jx <- matrix(0,nrow(x),ncol(x))
  Jy <- matrix(0,nrow(y),ncol(y))
  Jz <- matrix(0,nrow(z),ncol(z))
  
  # Dimensions:
  m <- nrow(x) + nrow(y); n <- ncol(x) + ncol(z)
  mx <- nrow(x); nx <- ncol(x);
  if (r > 0) {
    # vtilde: starting value for vtilde
    ztilde <- cbind(x,z)
    vtilde <- svd(ztilde,nu=r,nv=r)$v
    v <- vtilde[1:nx,]
    sx <- diag(1,r)
    uy <- t(solve(t(v) %*% v)%*%t(v)%*%t(y))
  } else {
    sx <- NULL
  }
  
  current <- matrix(0,m,n)
  sse <- c()
  convs <- c()
  
  for (i in 1:maxiter) {
    ### Individual Structure
    x_i <- (x - Jx)
    y_i <- (y - Jy)
    z_i <- (z - Jz)
    
    if (rx > 0) {
      s <- svd(x_i, nu=rx, nv=rx)
      Ax <- s$u[,1:rx] %*% diag(x=s$d[1:rx], nrow=rx) %*% t(s$v[,1:rx])
    }
    if (ry > 0) {
      s <- svd(y_i, nu=ry, nv=ry)
      Ay <- s$u[,1:ry] %*% diag(x=s$d[1:ry], nrow=ry) %*% t(s$v[,1:ry])
    }
    if (rz > 0) {
      s <- svd(z_i, nu=rz, nv=rz)
      Az <- s$u[,1:rz] %*% diag(x=s$d[1:rz], nrow=rz) %*% t(s$v[,1:rz])
    }
    
    sse <- c(sse,(norm(x-Jx-Ax,type='f')^2 + norm(y-Jy-Ay,type='f')^2 + norm(z-Jz-Az,type='f')^2))

    if (r > 0) {
      ### Joint Structure
      x_j <- x - Ax
      y_j <- y - Ay
      z_j <- z - Az
      # Combined data sets
      xtilde <- rbind(x_j,y_j)
      ztilde <- cbind(x_j,z_j)
      xvec <- as.vector(x_j)
      
      u <- t(solve(t(vtilde) %*% vtilde)%*%t(vtilde)%*%t(ztilde))
      # Scale u
      for (j in 1:r) {
        u[,j] <- u[,j] / norm(matrix(u[,j],ncol=1), type="f")
      }
      utilde <- rbind(u %*% sx, uy)
      v <- t(solve(t(utilde) %*% utilde)%*%t(utilde)%*%xtilde)
      vz <- t(solve(t(u) %*% u)%*%t(u)%*%z_j)
      # Scale v
      for (j in 1:r) {
        v[,j] <- v[,j] / norm(matrix(v[,j],ncol=1), type="f")
      }
      #vtilde <- rbind(v %*% sx, vz)
      uy <- t(solve(t(v) %*% v)%*%t(v)%*%t(y_j))
      
      # Calculate scaling matrix Sx
      uv <- matrix(0,mx*nx,r)
      cnt <- 1
      for (k in 1:nx) {
        for (j in 1:mx) {
          uv[cnt,] <- u[j,] * v[k,]
          cnt <- cnt + 1
        }
      }
      sx_vec <- as.vector(solve(t(uv) %*% uv) %*% t(uv) %*% xvec)
      sx <- diag(sx_vec, r)
      # Re-estimate (to incorporate scaling)
      utilde <- rbind(u %*% sx, uy)
      vtilde <- rbind(v %*% sx, vz)
      
      # Set joint estimates
      Jx <- u %*% sx %*% t(v)
      Jy <- uy %*% t(v)
      Jz <- u %*% t(vz)
    } else {
      # Set joint estimates to 0 when r=0
      Jx <- matrix(0,nrow(x),ncol(x))
      Jy <- matrix(0,nrow(y),ncol(y))
      Jz <- matrix(0,nrow(z),ncol(z))
    }
    
    sse_current <- norm(x-Jx-Ax,type='f')^2 + norm(y-Jy-Ay,type='f')^2 + norm(z-Jz-Az,type='f')^2
    
    if (anneal) {
      ### Simulated annealing
      #Tmax <- .1
      if (i <= nanneal) {
        rx <- matrix(rnorm(nrow(Jx)*ncol(Jx),0,Tmax/i),nrow(Jx),ncol(Jx))
        ry <- matrix(rnorm(nrow(Jx)*ncol(Jx),0,Tmax/i),nrow(Jy),ncol(Jy))
        rz <- matrix(rnorm(nrow(Jx)*ncol(Jx),0,Tmax/i),nrow(Jz),ncol(Jz))
        
        sseA <- norm(x-Jx-rx-Ax,type='f')^2 + norm(y-Jy-ry-Ay,type='f')^2 + norm(z-Jz-rz-Az,type='f')^2
        p <- exp(-(sseA-sse_current)/(Tmax/i))
        #   if (sseA < sse_current | runif(1) < p) {
        temp <- svd(Jx + rx, nu=r, nv=r)
        Jx <- temp$u[,1:r] %*% t(temp$v[,1:r]) 
        temp <- svd(Jy + ry, nu=r, nv=r)
        Jy <- temp$u[,1:r] %*% t(temp$v[,1:r]) 
        temp <- svd(Jz + rz, nu=r, nv=r)
        Jz <- temp$u[,1:r] %*% t(temp$v[,1:r])
        
        #      print("accepted") 
        #   }
      }
    }
    
    ### Check convergence
    if (i > 1) { 
      conv <- (norm(Jx-last[[1]],type='f')^2 + norm(Jy-last[[2]],type='f')^2 + norm(Jz-last[[3]],type='f')^2 +
        norm(Ax-last[[4]],type='f')^2 + norm(Ay-last[[5]],type='f')^2 + norm(Az-last[[6]],type='f')^2) / 
        (norm(x,type='f')^2 + norm(y,type='f')^2 + norm(z,type='f')^2)
    } else {
      conv <- (norm(Jx,type='f')^2 + norm(Jy,type='f')^2 + norm(Jz,type='f')^2 +
        norm(Ax,type='f')^2 + norm(Ay,type='f')^2 + norm(Az,type='f')^2) / 
        (norm(x,type='f')^2 + norm(y,type='f')^2 + norm(z,type='f')^2)
    }
    convs <- c(convs,conv)
    last <- list(Jx,Jy,Jz,Ax,Ay,Az)
    sse <- c(sse,(norm(x-Jx-Ax,type='f')^2 + norm(y-Jy-Ay,type='f')^2 + norm(z-Jz-Az,type='f')^2))
    print(paste(conv,(norm(x-Jx-Ax,type='f')^2 + norm(y-Jy-Ay,type='f')^2 + norm(z-Jz-Az,type='f')^2),sep=" "))
    if (conv < tol) { break }
  }
  
  temp <- svd(Jx)
  u <- temp$u[,1:r]
  v <- temp$v[,1:r]
  
  Jy_perp = Jy + (Ay %*% v %*% t(v))
  Jz_perp = Jz + (u %*% t(u) %*% Az)
  Ay_perp = Ay %*% (diag(ncol(x)) - v %*% t(v))
  Az_perp = (diag(nrow(x)) - u %*% t(u)) %*% Az
  
  Jy=Jy_perp; Jz=Jz_perp;
  Ay=Ay_perp; Az=Az_perp;
  
  return(list(Jx=Jx,Jy=Jy,Jz=Jz,Ax=Ax,Ay=Ay,Az=Az,sx=sx,sse=sse,convs=convs))
}


