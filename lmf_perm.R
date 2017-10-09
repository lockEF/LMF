# LMF with Only Joint Structure
lmfJ <- function(x, y, z, r, tol=.00001, maxiter=5000, scale=TRUE) {
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
  
  # Special Case r=0
  if (r == 0) {
    Jx <- Ax; Jy <- Ay; Jz <- Az
    u <- NULL; sx <- NULL; v <- NULL
    sse <- norm(x,type='f')^2 + norm(y,type='f')^2 + norm(z,type='f')^2; convs = NA
    return(list(Jx=Jx,Jy=Jy,Jz=Jz,Ax=Ax,Ay=Ay,Az=Az,u=u,sx=sx,v=v,sse=sse,convs=convs))
  }
  
  # Dimensions:
  m <- nrow(x) + nrow(y); n <- ncol(x) + ncol(z)
  mx <- nrow(x); nx <- ncol(x);
  # vtilde: starting value for vtilde
  ztilde <- cbind(x,z)
  vtilde <- svd(ztilde,nu=r,nv=r)$v
  v <- vtilde[1:nx,]
  sx <- diag(1,r)
  uy <- t(solve(t(v) %*% v)%*%t(v)%*%t(y))
  
  current <- matrix(0,m,n)
  sse <- c()
  convs <- c()
  
  for (i in 1:maxiter) {
    ### Joint Structure
    x_j <- x
    y_j <- y
    z_j <- z
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
    
    ### Check convergence
    if (i > 1) { 
      conv <- norm(Jx-last[[1]],type='f')^2 + norm(Jy-last[[2]],type='f')^2 + norm(Jz-last[[3]],type='f')^2
    } else {
      conv <- norm(Jx,type='f')^2 + norm(Jy,type='f')^2 + norm(Jz,type='f')^2
    }
    convs <- c(convs,conv)
    last <- list(Jx,Jy,Jz)
    sse <- c(sse,(norm(x-Jx,type='f')^2 + norm(y-Jy,type='f')^2 + norm(z-Jz,type='f')^2))
    if (conv < tol) { break }
  }
  
  return(list(Jx=Jx,Jy=Jy,Jz=Jz,Ax=Ax,Ay=Ay,Az=Az,u=u,sx=sx,v=v,sse=sse,convs=convs))
}


# LMF with Joint and Individual Structure
lmf <- function(x, y, z, r, rx, ry, rz, tol=.00001, maxiter=5000, scale=T, orth=0, anneal=F, nanneal=500, Tmax=0.1) {
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
      conv <- norm(Jx-last[[1]],type='f')^2 + norm(Jy-last[[2]],type='f')^2 + norm(Jz-last[[3]],type='f')^2 +
        norm(Ax-last[[4]],type='f')^2 + norm(Ay-last[[5]],type='f')^2 + norm(Az-last[[6]],type='f')^2
    } else {
      conv <- norm(Jx,type='f')^2 + norm(Jy,type='f')^2 + norm(Jz,type='f')^2 +
        norm(Ax,type='f')^2 + norm(Ay,type='f')^2 + norm(Az,type='f')^2
    }
    convs <- c(convs,conv)
    last <- list(Jx,Jy,Jz,Ax,Ay,Az)
    sse <- c(sse,(norm(x-Jx-Ax,type='f')^2 + norm(y-Jy-Ay,type='f')^2 + norm(z-Jz-Az,type='f')^2))
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


# Function to return a vector norm
vnorm <- function(x) {
  return(sqrt(sum(x^2)))
}


# Permutation test for rank selection (ad hoc method)
lmf.perm <- function(x, y, z, nperms=100, alpha=.05) {
  # Center and Scale data
  x <- x - matrix(rep(apply(x,1,mean,na.rm=T),ncol(x)),nrow=nrow(x))
  x <- x / norm(x, type="f")
  y <- y - matrix(rep(apply(y,1,mean,na.rm=T),ncol(y)),nrow=nrow(y))
  y <- y / norm(y, type="f")
  z <- z - matrix(rep(apply(z,1,mean,na.rm=T),ncol(z)),nrow=nrow(z))
  z <- z / norm(z, type="f")
  
  # Set maxr (include as function param later)
  #maxr=min(mx,nx)
  maxr <- 5
  
  # Initializations for permutation test
  nrun <- 0
  last <- rep(-2, 4)
  current <- rep(-1, 4)
  Jx <- matrix(0, nrow=nrow(x), ncol=ncol(x))
  Jy <- matrix(0, nrow=nrow(y), ncol=ncol(y))
  Jz <- matrix(0, nrow=nrow(z), ncol=ncol(z))
  Ax <- matrix(0, nrow=nrow(x), ncol=ncol(x))
  Ay <- matrix(0, nrow=nrow(y), ncol=ncol(y))
  Az <- matrix(0, nrow=nrow(z), ncol=ncol(z))
  
  while (!isTRUE(all.equal(last, current)) & nrun < 10) {
    # STEP THROUGH TO SEE IF FIRST ITERATION ALSO UNDERESTIMATES JOINT STRUCTURE
    last <- current
    
    # Call joint factorization on data with ind structure removed
    xj <- x - Ax
    yj <- y - Ay
    zj <- z - Az
    actual <- Jperm2(xj,yj,zj)
    
    # Call on data with permuted cols Y, rows Z
    # for each permutation
    perm <- matrix(NA,maxr,nperms)
    mx <- nrow(x); nx <- ncol(x);
    for(i in 1:nperms) {
      # Doing these simultaneously makes sense
      Yperm <- yj[, sample(1:nx, nx, replace=F)]
      Zperm <- zj[sample(1:mx, mx, replace=F), ]
      perm[,i] <- Jperm2(xj,Yperm,Zperm)
    }
    
    # Get ranks
    r <- 0
    for (i in 1:maxr) {
      if (actual[i] > quantile(perm[i,], 1-alpha)) {
        r <- r + 1
      } else {
        break
      }
    }
    # Could try (in initial step only): 1) estimate rank(J)   2) estimate rank(Xi)   3) rank(Ai) = rank(Xi) - rank(J)
    # Update Eric by Sunday night on perm simulation
    
    # Individual rank permutation test
    xi <- x - Jx
    yi <- y - Jy
    zi <- z - Jz
    
    # X (permute rows and columns)
    actual <- svd(xi,nu=0,nv=0)$d[1:maxr]
    perm <- matrix(NA,maxr,nperms)
    for (i in 1:nperms) {
      temp <- matrix(sample(xi, nrow(xi) * ncol(xi)), nrow(xi), ncol(xi))
      perm[,i] <- svd(temp,nu=0,nv=0)$d[1:maxr]
    }
    rx <- 0
    for (i in 1:maxr) {
      if (actual[i] > quantile(perm[i,], 1-alpha)) { 
        rx <- rx + 1
      } else {
        break
      }
    }
    
    # Y (permute rows)
    actual <- svd(yi,nu=0,nv=0)$d
    perm <- matrix(NA,maxr,nperms)
    for (i in 1:nperms) {
      temp <- t(yi)
      # Permute Rows
      pind <- order(c(col(temp)), runif(length(temp))) 
      temp <- matrix(temp[pind], nrow=nrow(yi), ncol=ncol(yi), byrow=TRUE)
      perm[,i] <- svd(temp,nu=0,nv=0)$d[1:maxr]
    }
    ry <- 0
    for (i in 1:maxr) {
      if (actual[i] > quantile(perm[i,], 1-alpha)) { 
        ry <- ry + 1
      } else {
        break
      }
    }
    
    # Z (permute columns)
    actual <- svd(zi,nu=0,nv=0)$d
    perm <- matrix(NA,maxr,nperms)
    for (i in 1:nperms) {
      temp <- zi
      # Permute Cols
      pind <- order(c(col(temp)), runif(length(temp))) 
      temp <- matrix(temp[pind], nrow=nrow(zi), ncol=ncol(zi), byrow=FALSE)
      perm[,i] <- svd(temp,nu=0,nv=0)$d[1:maxr]
    }
    rz <- 0
    for (i in 1:maxr) {
      if (actual[i] > quantile(perm[i,], 1-alpha)) { 
        rz <- rz + 1
      } else {
        break
      }
    }
    
    fit <- lmf(x,y,z,r=r,rx=rx,ry=ry,rz=rz)
    Jx=fit$Jx; Jy=fit$Jy; Jz=fit$Jz;
    Ax=fit$Ax; Ay=fit$Ay; Az=fit$Az;
    current <- c(r,rx,ry,rz)
    if (nrun == 0) {
      rx <- rx - r
      ry <- ry - r
      rz <- rz - r
    }
    nrun <- nrun + 1
  } # (end of while loop)
  
  #return(c(fit,r))
  return(c(r,rx,ry,rz))
}

### TRY CROSS-VALIDATION APPROACH
# Pseudocode
# 1) Partition Rows and Columns of Data
# 2) Estimate models of rank 1-(5) for each joint and individual rank (20 models for example case) using all but 1 partition
# 3) Use partitioned sample to compare models
# 4) Repeat for other partitions of sample
JpermCV <- function(x,y,z) {
   # Choose partitions
   # Estimate models
}


# Joint Permutation Test Statistic
## Look at gain in SSE for each model
Jperm <- function(x,y,z) {
   temp <- lmfJ(x,y,z,1)
   ss <- norm(temp$Jx,type='f')^2 + norm(temp$Jy,type='f')^2 + norm(temp$Jz,type='f')^2
   for (i in 2:5) {
      temp <- lmfJ(x,y,z,i)
      ss <- c(ss, norm(temp$Jx,type='f')^2 + norm(temp$Jy,type='f')^2 + norm(temp$Jz,type='f')^2 - sum(ss[1:(i-1)]))
   }
   return(ss)
}

# Joint Permutation Test Statistic (method 2)
## Look at gain in SSE for each model
Jperm2 <- function(x,y,z) {
  # Center and Scale data
  x <- x - matrix(rep(apply(x,1,mean,na.rm=T),ncol(x)),nrow=nrow(x))
  x <- x / norm(x, type="f")
  y <- y - matrix(rep(apply(y,1,mean,na.rm=T),ncol(y)),nrow=nrow(y))
  y <- y / norm(y, type="f")
  z <- z - matrix(rep(apply(z,1,mean,na.rm=T),ncol(z)),nrow=nrow(z))
  z <- z / norm(z, type="f")
  
  temp <- lmfJ(x,y,z,1,scale=F)
  ss <- norm(x,type='f')^2 + norm(y,type='f')^2 + norm(z,type='f')^2 -
    (norm(x - temp$Jx,type='f')^2 + norm(y - temp$Jy,type='f')^2 + norm(z - temp$Jz,type='f')^2)
  x_cum <- x - temp$Jx
  y_cum <- y - temp$Jy
  z_cum <- z - temp$Jz
  for (i in 2:5) {
    temp <- lmfJ(x_cum,y_cum,z_cum,1,scale=F)
    ss <- c(ss, norm(x_cum,type='f')^2 + norm(y_cum,type='f')^2 + norm(z_cum,type='f')^2 -
              (norm(x_cum - temp$Jx,type='f')^2 + norm(y_cum - temp$Jy,type='f')^2 + norm(z_cum - temp$Jz,type='f')^2))
    x_cum <- x_cum - temp$Jx
    y_cum <- y_cum - temp$Jy
    z_cum <- z_cum - temp$Jz
  }
  return(ss)
}


