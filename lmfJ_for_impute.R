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
