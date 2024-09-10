find_circle <- function(A) {
  # A: adjacency matrix
  # output: the first circle if found, otherwise return NULL
  
  p <- nrow(A)
  for (start in 1:p) {
    visited <- NULL
    stack <- list(start)
    while (length(stack) > 0) {
      v <- stack[[length(stack)]]
      stack <- stack[-length(stack)]  
      
      if (v %in% visited) {
        C <- NULL
        while (!(v %in% C)) {
          C <- c(C, v)
          v <- which(A[,v]==1)  
        }
        return(C)
      }
      
      visited <- c(visited, v)
      if (v %in% 1:p) {
        stack <- c(stack, which(A[v,]==1))
      }
    }
  }
  return(NULL)
}

directedMST <- function(W, root) {
  # Edmonds' algorithm 
  # W: pairwise weight matrix (we assume a fully connected graph)
  # root: should be pre-specified
  # output: adjacency matrix 
  
  # prepare
  last_nd = p <- nrow(W)
  nd_set <- 1:p
  non_r <- nd_set[-root]                                                        # index without root
  W_tmp <- W
  A_tmp <- matrix(0, p, p)                                                      # no in-coming edges to the root
  A = C <- list()
  history_in = history_out <- list()
  
  ## step 1: find the parent with minimum weight 
  pa <- apply(W_tmp[,non_r], 2, which.min)
  A_tmp[cbind(pa,non_r)] <- 1
  
  ## check if A contains circle(s) 
  C_tmp <- find_circle(A_tmp)
  if (is.null(C_tmp)) { return(A_tmp) } 
  
  # recursive starts
  while (!is.null(C_tmp)) {
    C[[last_nd-(p-1)]] <- C_tmp
    A_old <- A_tmp
    A[[last_nd-(p-1)]] <- A_old
    
    ## step 2: create new A and W with vc 
    v_c <- last_nd + 1
    nd_set <- nd_set[!(nd_set %in% C_tmp)]
    non_r <- nd_set[-root]
    nd_count <- length(nd_set)
    
    A_tmp <- matrix(0, v_c, v_c)
    W_tmp <- cbind(rbind(W_tmp, Inf), Inf)
    
    ### edges going away from the cycle 
    w_out <- matrix(0, length(C_tmp), nd_count-1, 
                     dimnames = list(C_tmp, non_r))
    for (m in non_r) {
      for (l in C_tmp) {
        w_out[as.character(l),as.character(m)] <- W_tmp[l,m] 
      }
    }
    W_tmp[v_c,non_r] <- apply(w_out, 2, min)
    
    ### edges coming in to the cycle 
    w_in <- matrix(0, nd_count, length(C_tmp), 
                   dimnames = list(nd_set, C_tmp))
    for (m in nd_set) {
      for (l in C_tmp) {
        w_in[as.character(m),as.character(l)] <- W_tmp[m,l] - W_tmp[which(A_old[,l]==1),l]
      }
    }
    W_tmp[nd_set,v_c] <- apply(w_in, 1, min)
    
    ### history save
    history_in[[last_nd-(p-1)]] <- cbind(nd_set, 
                                      C_tmp[apply(w_in, 1, which.min)])
    history_out[[last_nd-(p-1)]] <- cbind(C_tmp[apply(w_out, 2, which.min)], 
                                           non_r)
    
    ## back to step 1
    pa <- apply(as.matrix(W_tmp[c(nd_set,v_c),c(non_r,v_c)]), 2, which.min)
    A_tmp[cbind(c(nd_set,v_c)[pa],c(non_r,v_c))] <- 1
    
    ## check if A contains circle(s) 
    C_tmp <- find_circle(A_tmp)
    last_nd <- last_nd + 1                                                      # possible to add one more node if C_tmp is not null
    nd_set <- c(nd_set, v_c)
  }
  
  ## step 3: marking
  k <- length(history_in)
  A_out <- A_tmp                                                                # naturally include original node -> original node

  while (k > 0) {
    A_old <- A[[k]]
    
    ### edges to the cycle
    edges_in <- history_in[[k]]
    u <- which(A_out[1:v_c,v_c] == 1)                         
    v <- edges_in[which(edges_in[,1] == u),2]
    A_out[u,v] <- 1
    
    ### mark edges in C
    A_out[C[[k]],C[[k]]] <- A_old[C[[k]],C[[k]]]
    
    ### remove the edge towards v in C 
    A_out[-u,v] <- rep(0, last_nd-1)
    
    ### edges away from C
    edges_out <- history_out[[k]]
    if (nrow(edges_out) > 0) {
      for (i in 1:nrow(edges_out)) {
        if (A_out[v_c,edges_out[i,2]] == 1) {
          A_out[edges_out[i,1],edges_out[i,2]] <- 1
        }
      }
    }
    
    k <- k-1
    v_c <- v_c-1
  }
  
  return(A_out[1:p,1:p])
}

