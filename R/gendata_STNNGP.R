rstnngp <- function(coords, 
                    params = list(tau.sq = tau.sq,
                                  sigma.sq = sigma.sq,
                                  phi = phi, 
                                  rho = rho,
                                  beta = beta, 
                                  A = A, 
                                  n.neighbors = n.neighbors), 
                    seed = seed) {
  # coords: nx2
  # beta: pxq matrix of coefficients
  # tau.sq: q vector, nugget
  # sigma.sq: scalar, spatial variability
  # phi: scalar, decay
  # rho: q-1 vector, correlation coefficients
  # A: qxq adjacency matrix of a spanning tree
  
  # set seed
  set.seed(seed)
  
  # parameters
  tau.sq <- params$tau.sq
  sigma.sq <- params$sigma.sq
  phi <- params$phi
  rho <- params$rho
  beta <- params$beta
  A <- params$A
  m <- params$n.neighbors
  
  # ordering
  ord <- order(coords[,1])
  coords <- coords[ord,]
  
  # dimensions
  n <- nrow(coords)
  D <- as.matrix(dist(coords))
  q <- length(tau.sq)
  
  # exponential covariance
  M <- sigma.sq*exp(-phi*D)
  
  # covariates
  p <- nrow(beta)
  X <- cbind(1, matrix(rnorm(n*(p-1)), n, p-1))
  
  # neighbors
  dir_nn_idx = undir_nn_idx <- list()
  dir_nn_idx[[1]] <- NULL
  for (i in 1:n) {
    if (i > 1) {
      dir_nn <- RANN::nn2(data = matrix(coords[1:(i-1),], ncol = 2),              # where
                          query = matrix(coords[i,], ncol = 2),                   # whose
                          k = min(m, i-1))                                        # how many
      dir_nn_idx[[i]] <- as.numeric(dir_nn$nn.idx)
    }
    undir_nn <- RANN::nn2(data = coords,                                          # where
                          query = matrix(coords[i,], ncol = 2),                   # whose
                          k = m+1)                                                # how many, including itself
    undir_nn_idx[[i]] <- as.numeric(undir_nn$nn.idx)
  }
  
  # generate Y
  cov_list <- list()
  for (i in 1:n) {
    nn_idx_i <- c(i,undir_nn_idx[[i]],dir_nn_idx[[i]])
    cov_list[[i]] <- M[nn_idx_i, nn_idx_i]
  }
  
  Y <- matrix(0, n, q)
  Y[1,1] <- sqrt(sigma.sq + tau.sq[1])*rnorm(1)
  for (i in 2:n) {
    cov_mat <- cov_list[[i]][c(1,(m+2)+1:min(m,i-1)), c(1,(m+2)+1:min(m,i-1))] +
      diag(tau.sq[1], 1+min(m,i-1))
    
    coef <- cov_mat[1,-1]%*%solve(cov_mat[-1,-1])
    Y[i,1] <- coef%*%Y[dir_nn_idx[[i]],1] + 
      sqrt(cov_mat[1,1]-coef%*%cov_mat[-1,1])*rnorm(1)
  }
  
  for (k in 2:q) {
    j <- which(A[,k] == 1)
    for (i in 1:n) {
      cov_mat <- cov_list[[i]]
      size <- nrow(cov_mat)
      rho_mat <- matrix(1, size, size)
      rho_mat[1, 2:(m+2)] <- rho[k-1]
      rho_mat[2:(m+2),1] <- rho[k-1]
      if (i > 1) {
        rho_mat[2:(m+2), (m+3):size] <- rho[k-1] 
        rho_mat[(m+3):size, 2:(m+2)] <- rho[k-1]
      }
      cov_mat <- rho_mat*cov_mat + 
        diag(c(tau.sq[k], rep(tau.sq[j], m+1), rep(tau.sq[k], size-m-2)))
      
      coef <- cov_mat[1,-1]%*%solve(cov_mat[-1,-1])
      Y[i,k] <- coef%*%c(Y[undir_nn_idx[[i]],j],Y[dir_nn_idx[[i]],k]) + 
        sqrt(cov_mat[1,1]-coef%*%cov_mat[-1,1])*rnorm(1)
    }
  }
  
  Y <- X%*%beta + Y
  
  out <- list()
  out$X <- X[order(ord),,drop=FALSE]
  out$Y <- Y[order(ord),,drop=FALSE]
  
  return(out)
}

rstnngp_ns <- function(coords, 
                    params = list(tau.sq = tau.sq,
                                  sigma.sq = sigma.sq,
                                  phi = phi, 
                                  cross_phi = cross_phi,
                                  rho = rho,
                                  beta = beta, 
                                  A = A, 
                                  n.neighbors = n.neighbors), 
                    seed = seed) {
  # coords: nx2
  # beta: pxq matrix of coefficients
  # tau.sq: q vector, nugget
  # sigma.sq: q vector, spatial variability
  # phi: q vector, decay
  # rho: q-1 vector, correlation coefficients
  # A: qxq adjacency matrix of a spanning tree
  
  # set seed
  set.seed(seed)
  
  # parameters
  tau.sq <- params$tau.sq
  sigma.sq <- params$sigma.sq
  phi <- params$phi
  cross_phi <- params$cross_phi
  rho <- params$rho
  beta <- params$beta
  A <- params$A
  m <- params$n.neighbors
  
  # dimensions
  n <- nrow(coords)
  q <- length(tau.sq)
  
  # generate cross-covariance parameters
  cross_sigma <- rep(0, q-1)
  for (me in 2:q) {
    parent <- which(A[,me] == 1)
    cross_sigma[me-1] <- sqrt(sigma.sq[parent]*sigma.sq[me])
  }
  
  # ordering
  ord <- order(coords[,1])
  coords <- coords[ord,]
  D <- as.matrix(dist(coords))
  
  # exponential covariance
  marg_cov = cross_cov <- list()
  for (j in 1:q) {
    marg_cov[[j]] <- sigma.sq[j]*exp(-phi[j]*D) + diag(tau.sq[j],n)
  }
  for (j in 2:q) {
    cross_cov[[j-1]] <- rho[j-1]*cross_sigma[j-1]*exp(-cross_phi[j-1]*D)
  }
  
  # covariates
  p <- nrow(beta)
  X <- cbind(1, matrix(rnorm(n*(p-1)), n, p-1))
  
  # neighbors
  dir_nn_idx = undir_nn_idx <- list()
  dir_nn_idx[[1]] <- NULL
  for (i in 1:n) {
    if (i > 1) {
      dir_nn <- RANN::nn2(data = matrix(coords[1:(i-1),], ncol = 2),              # where
                          query = matrix(coords[i,], ncol = 2),                   # whose
                          k = min(m, i-1))                                        # how many
      dir_nn_idx[[i]] <- as.numeric(dir_nn$nn.idx)
    }
    undir_nn <- RANN::nn2(data = coords,                                          # where
                          query = matrix(coords[i,], ncol = 2),                   # whose
                          k = m+1)                                                # how many, including itself
    undir_nn_idx[[i]] <- as.numeric(undir_nn$nn.idx)
  }
  
  # generate Y
  cov_list <- list()
  for (j in 1:q) {
    covj_list <- list()
    if (j == 1) {
      for (i in 1:n) {
        covj_list[[i]] <- marg_cov[[j]][c(i,dir_nn_idx[[i]]), c(i,dir_nn_idx[[i]])]
      }
    } else {
      parent <- which(A[,j]==1)
      for (i in 1:n) {
        size <- 1+(m+1)+min(m,i-1)
        covj_tmp <- matrix(0, size, size)
        covj_tmp[1,1+1:(m+1)] = covj_tmp[1+1:(m+1),1] <- cross_cov[[j-1]][i,undir_nn_idx[[i]]]
        covj_tmp[1+1:(m+1),1+1:(m+1)] <- marg_cov[[parent]][undir_nn_idx[[i]],undir_nn_idx[[i]]]
        if (i == 1) {
          covj_tmp[1,1] <- marg_cov[[j]][i,i]
        } else {
          covj_tmp[c(1,(m+2)+1:min(m, i-1)),c(1,(m+2)+1:min(m, i-1))] = 
            marg_cov[[j]][c(i,dir_nn_idx[[i]]),c(i,dir_nn_idx[[i]])]
          covj_tmp[1+1:(m+1),(m+2)+1:min(m, i-1)] = cross_cov[[j-1]][undir_nn_idx[[i]], dir_nn_idx[[i]]]
          covj_tmp[(m+2)+1:min(m, i-1),1+1:(m+1)] = cross_cov[[j-1]][dir_nn_idx[[i]], undir_nn_idx[[i]]]
        }
        covj_list[[i]] <- covj_tmp
      }
    }
    cov_list[[j]] <- covj_list
  }
  
  Y <- matrix(0, n, q)
  Y[1,1] <- sqrt(sigma.sq[1] + tau.sq[1])*rnorm(1)
  for (i in 2:n) {
    cov_mat <- cov_list[[1]][[i]]
    
    coef <- cov_mat[1,-1]%*%solve(cov_mat[-1,-1])
    Y[i,1] <- coef%*%Y[dir_nn_idx[[i]],1] + 
      sqrt(cov_mat[1,1]-coef%*%cov_mat[-1,1])*rnorm(1)
  }
  
  for (k in 2:q) {
    j <- which(A[,k] == 1)
    for (i in 1:n) {
      cov_mat <- cov_list[[k]][[i]] 
      
      coef <- cov_mat[1,-1]%*%solve(cov_mat[-1,-1])
      Y[i,k] <- coef%*%c(Y[undir_nn_idx[[i]],j],Y[dir_nn_idx[[i]],k]) + 
        sqrt(cov_mat[1,1]-coef%*%cov_mat[-1,1])*rnorm(1)
    }
  }
  
  Y <- X%*%beta + Y
  
  out <- list()
  out$X <- X[order(ord),,drop=FALSE]
  out$Y <- Y[order(ord),,drop=FALSE]
  
  return(out)
}