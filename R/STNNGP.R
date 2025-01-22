STNNGP <- function(formula, data = parent.frame(), coords, method = "response", 
                   mv.model = "separable",                                      # BJ
                   family = "gaussian", weights, n.neighbors = 15, 
                   starting, tuning, priors, cov.model = "exponential",
                   n.samples, n.omp.threads = 1, search.type = "cb", ord, 
                   return.neighbor.info = FALSE, neighbor.info,
                   verbose = TRUE, n.report = 100, 
                   debug = list(),                                              # BJ
                   ...){
  
  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- c(names(formals(sys.function(sys.parent()))), "u.search.type")
  
  elip.args <- list(...)
  for(i in names(elip.args)){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
  
  # BJ: debug
  if (is.null(debug$update.sigma.sq)) {
    update.sigma.sq <- TRUE
  } else {
    update.sigma.sq <- debug$update.sigma.sq
  }
  if (is.null(debug$update.w)) {
    update.w <- TRUE
  } else {
    update.w <- debug$update.w
  }
  if (is.null(debug$update.beta)) {
    update.beta <- TRUE
  } else {
    update.beta <- debug$update.beta
  }
  storage.mode(update.sigma.sq) <- "integer"
  storage.mode(update.w) <- "integer"
  storage.mode(update.beta) <- "integer"
  
  ####################################################
  ##Formula
  ####################################################
  if(missing(formula)){stop("error: formula must be specified")}
  
  if(is(formula, "formula")){
    if (formula[[2]] != "Y") {stop("error: response should be 'Y'.")}           # BJ
    
    # BJ: misaligned case 
    # We assume X is common for all variables and has no missingness
    if(missing(data)) data <- sys.frame(sys.parent())
    Y = Y_org <- data$Y
    misalign <- FALSE
    if (sum(is.na(Y)) > 0) { 
      misalign <- TRUE 
      if(method == "response") {stop("error: only the latent method can be used for misaligned data.")}
      Y[is.na(Y)] <- -9999                                                      # BJ: not accessed
      data$Y <- Y
      print("Misaligned multivariate data is provided.")
    }
    
    holder <- parseFormula(formula, data)
    Y <- as.matrix(holder[[1]])                                                 # BJ
    X <- as.matrix(holder[[2]])
    x.names <- holder[[3]]
    
    if (misalign) {
      data$Y <- Y_org
    }
  }else{
    stop("error: formula is misspecified")
  }
  
  q <- ncol(Y)                                                                  # BJ
  p <- ncol(X)
  n <- nrow(X)
  
  ##Coords
  if(missing(coords)){stop("error: coords must be specified")}
  
  if(is.vector(coords)){
    if(is.character(coords)){
      if(all(coords %in% colnames(data))){
        coords <- as.matrix(data[,coords])
      }else{
        stop(paste0("error: coords name ", paste(coords[!(coords %in% colnames(data))], collapse=" and "), " not in data"))
      }
    }else if(all(coords %in% (1:ncol(data)))){
      coords <- as.matrix(data[,coords])
    }else{
      stop(paste0("error: coords column index ", paste(coords[!(coords %in% (1:ncol(data)))], collapse=" and "), " not in data"))
    }
  }else{
    if(!any(is.matrix(coords), is.data.frame(coords))){
      stop("error: coords must n-by-2 matrix or dataframe of coordinates or vector indicating the column names or integer indexes in data") # BJ
    }
    coords <- as.matrix(coords)
  }
  
  if(nrow(coords) != n || ncol(coords) != 2){
    stop("error: either coords has more than two columns or the number of rows is different than data used in the model formula") # BJ
  }
  
  ## This can be slow if n is large, the onus is on the user to check
  ## if(any(duplicated(coords))){
  ##     stop("error: duplicated coordinates found. Remove duplicates.")
  ## }
  
  ####################################################
  ##Family
  #################################################### 
  family.names <- c("gaussian","binomial")
  family <- tolower(family)
  
  if(!family%in%family.names){stop("error: specified family '",family,"' is not a valid option; choose from ", paste(family.names, collapse=", ", sep="") ,".")}    
  if(method == "response" & family != "gaussian"){stop("error: only the latent method can be used with non-Gaussain family.")}
  
  if(family == "binomial"){
    if(missing(weights)){
      weights <- rep(1, n*q)
    }else if(length(weights) == 1){
      weights <- rep(weights, n*q)
    }else{
      weights <- c(weights)
      if(length(weights) != n*q){
        stop("error: for non-Gaussian models weights must be a n-by-q matrix or a scalar (in which case the specified value is repeted nq times)")
      }
    }
  }else{
    weights <- NA
  }
  
  ####################################################
  ## BJ: Multivariate model 
  #################################################### 
  mv.model.names <- c("separable", "nonseparable")
  mv.model <- sub("-", "", x = tolower(mv.model))
  
  if(!mv.model %in% mv.model.names){stop("error: specified multivariate model '", mv.model, "' is not a valid option; choose from ", paste(mv.model.names, collapse=", ", sep="") ,".")}    
  
  ####################################################
  ##Neighbors and ordering
  ####################################################
  neighbor.info.provided <- FALSE
  
  u.search.type <- 2 ##2 is very fast, 1 is slightly faster than 0, and 0 is the original slow one (2 and 0 should match, 1 is also corrected just different opposite sorting among the children)
  if("u.search.type" %in% names(elip.args)){
    u.search.type <- elip.args[["u.search.type"]]
  }
  
  if(!missing(neighbor.info)){
    
    warning("Using user defined neighbor.info. No checks are done on the supplied neighbor information.")
    
    if(!all(c("n.neighbors","nn.indx","nn.indx.lu",
              "nn.dist", "nn.dist.parent",                                      # BJ
              "nn.indx.parent", "nn.indx.lu.parent", "nn.indx.lu.all",          # BJ
              "ord") %in% names(neighbor.info))){stop("The supplied neighbor.info is malformed.")}

    nn.indx <- neighbor.info$nn.indx
    nn.indx.lu <- neighbor.info$nn.indx.lu
    nn.dist <- neighbor.info$nn.dist                                            # BJ
    nn.dist.parent <- neighbor.info$nn.dist.parent                              # BJ
    nn.indx.parent <- neighbor.info$nn.indx.parent                              # BJ
    nn.indx.lu.parent <- neighbor.info$nn.indx.lu.parent                        # BJ
    nn.indx.lu.all <- neighbor.info$nn.indx.lu.all                              # BJ
    ord <- neighbor.info$ord
    n.neighbors <- neighbor.info$n.neighbors
    nn.indx.run.time <- neighbor.info$nn.indx.run.time
    neighbor.info.provided <- TRUE
    
    if(method == "latent"){
      
      if(!all(c("u.indx", "u.indx.lu", "ui.indx") %in% names(neighbor.info))){
        ##must be neighbor.info from a response or conjugate model
        neighbor.info <- c(neighbor.info, mkUIndx(n, n.neighbors, nn.indx, nn.indx.lu, u.search.type))
        print("Computing additional index needed for method = 'latent' using user defined neighbor.info. This might take a while if n is large.")
      }
      
      u.indx <- neighbor.info$u.indx
      u.indx.lu <- neighbor.info$u.indx.lu
      ui.indx <- neighbor.info$ui.indx
      u.indx.run.time <- neighbor.info$u.indx.run.time
    }
    
  }else{
    
    if(missing(ord)){   
      ord <- order(coords[,1])##default order by x column
    }else{
      if(length(ord) != n){stop("error: supplied order vector ord must be of length n")}
      ## if(search.type == "cb"){
      ##     warning("switching to search.type='brute' given user defined ordering ord, this could be slow if n is large")
      ##     search.type <- "brute"
      ## }
    }
  }
  
  coords <- coords[ord,]
  X <- X[ord,,drop=FALSE]
  Y <- Y[ord,,drop=FALSE]                                                       # BJ
  Xbd <- diag(q)%x%X                                                            # BJ
  
  if(family == "binomial"){
    weights <- weights[rep(0:(q-1)*n, each = n)+ord]                            # BJ                                          
  }
  
  # BJ: misaligned case 
  if (misalign) { 
    Y_missing <- matrix(c(1*(Y == -9999)), nrow = n, ncol = q)
    nj <- n - colSums(Y_missing)
    XtX <- matrix(0, nrow = p^2, ncol = q)
    for (j in 1:q) {
      Xtmp <- X*(1-Y_missing[,j])
      XtX[,j] <- crossprod(Xtmp)
    }
  } else {
    Y_missing <- matrix(0, nrow = n, ncol = q)
  }
  
  storage.mode(Y) <- "double"                                                   # BJ
  storage.mode(X) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"
  storage.mode(coords) <- "double"
  storage.mode(Y_missing) <- "integer"                                          # BJ
  if (misalign) {
    storage.mode(nj) <- "integer"                                               # BJ
    storage.mode(XtX) <- "double"                                               # BJ
  }
  storage.mode(weights) <- "integer"
  
  ####################################################
  ##NNGP method
  ####################################################
  method.names <- c("response","latent")
  method <- tolower(method)
  
  if(!method%in%method.names){stop("error: specified method '",method,"' is not a valid option; choose from ", paste(method.names, collapse=", ", sep="") ,".")}
  
  ####################################################
  ##Covariance model
  ####################################################
  if(missing(cov.model)){stop("error: cov.model must be specified")}
  
  cov.model.names <- c("exponential","spherical","matern","gaussian")##order much match util.cpp spCor
  
  if(!cov.model%in%cov.model.names)
  {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose from ", paste(cov.model.names, collapse=", ", sep="") ,".")}
  
  cov.model.indx <- which(cov.model == cov.model.names)-1##obo for cov model lookup on c side
  storage.mode(cov.model.indx) <- "integer"
  
  ####################################################
  ##Priors
  ####################################################
  if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
  
  names(priors) <- tolower(names(priors))
  
  if(!"sigma.sq.ig" %in% names(priors)){stop("error: sigma.sq.IG must be specified")}
  if (mv.model == "separable") {                                                # BJ
    sigma.sq.IG <- priors[["sigma.sq.ig"]]
    if(!is.vector(sigma.sq.IG) || length(sigma.sq.IG) != 2){
      stop("error: sigma.sq.IG must be a vector of length 2")
    }
    if(any(sigma.sq.IG <= 0)){stop("error: sigma.sq.IG must be positive")}      # BJ
    storage.mode(sigma.sq.IG) <- "double"
  } else {
    sigma.sq.IG <- as.matrix(priors[["sigma.sq.ig"]])                           # BJ
    if(ncol(sigma.sq.IG) != 2 || nrow(sigma.sq.IG) != q){                       # BJ
      stop("error: sigma.sq.IG must be a q-by-2 matrix")                        # BJ
    }                                                                           # BJ
    if(any(sigma.sq.IG <= 0)){stop("error: sigma.sq.IG must be positive")}      # BJ
    sigma.sq.IGa <- sigma.sq.IG[,1]                                             # BJ
    sigma.sq.IGb <- sigma.sq.IG[,2]                                             # BJ
    storage.mode(sigma.sq.IGa) = storage.mode(sigma.sq.IGb) <- "double"         # BJ
  }
  
  if(family == "gaussian"){
    if(!"tau.sq.ig" %in% names(priors)){stop("error: tau.sq.IG must be specified")}
    tau.sq.IG <- as.matrix(priors[["tau.sq.ig"]])                               # BJ
    if(ncol(tau.sq.IG) != 2 || nrow(tau.sq.IG) != q){                           # BJ
      stop("error: tau.sq.IG must be q-by-2 matrix")                            # BJ
    }                                                                           # BJ
    if(any(tau.sq.IG <= 0)){stop("error: tau.sq.IG must be positive")}          # BJ
    tau.sq.IGa <- tau.sq.IG[,1]                                                 # BJ
    tau.sq.IGb <- tau.sq.IG[,2]                                                 # BJ
    storage.mode(tau.sq.IGa) = storage.mode(tau.sq.IGb) <- "double"             # BJ
  }
  
  if(!"phi.unif" %in% names(priors)){stop("error: phi.Unif must be specified")}
  if (mv.model == "separable") {                                                # BJ
    phi.Unif <- priors[["phi.unif"]]
    if(!is.vector(phi.Unif) || length(phi.Unif) != 2){
      stop("error: phi.Unif must be a vector of length 2")
      }
    if(any(phi.Unif <= 0, phi.Unif[1] >= phi.Unif[2])){
      stop("error: phi.Unif must be a positive vector of length 2 with element 1 < element 2")
    }
    storage.mode(phi.Unif) <- "double"
  } else {
    phi.Unif <- as.matrix(priors[["phi.unif"]])                                 # BJ
    if(ncol(phi.Unif) != 2 || nrow(phi.Unif) != q){                             # BJ
      stop("error: phi.Unif must be a q-by-2 matrix")                           # BJ
    }                                                                           # BJ
    phi.Unifa <- phi.Unif[,1]                                                   # BJ
    phi.Unifb <- phi.Unif[,2]                                                   # BJ
    if(any(phi.Unif <= 0, phi.Unifa - phi.Unifb >= 0)){                         # BJ
      stop("error: phi.Unif must be positive with element 1 < element 2")       # BJ
    }                                                                           # BJ
    storage.mode(phi.Unifa) = storage.mode(phi.Unifb) <- "double"               # BJ
  }

  nu.Unif <- 0
  nu.Unifa = nu.Unifb <- rep(0, q)                                              # BJ
  if(cov.model == "matern"){
    if(!"nu.unif" %in% names(priors)){stop("error: nu.Unif must be specified")}
    if (mv.model == "separable") {                                              # BJ
      nu.Unif <- priors[["nu.unif"]]
      if(!is.vector(nu.Unif) || length(nu.Unif) != 2){
        stop("error: nu.Unif must be a vector of length 2")
        }
      if(any(nu.Unif <= 0, nu.Unif[1] >= nu.Unif[2])){
        stop("error: nu.Unif must be a positive vector of length 2 with element 1 < element 2")
        }
    } else {
      nu.Unif <- as.matrix(priors[["nu.unif"]])                                 # BJ
      if(ncol(nu.Unif) != 2 || nrow(nu.Unif) != q){                             # BJ
        stop("error: nu.Unif must be a q-by-2 matrix")                          # BJ
      }                                                                         # BJ
      nu.Unifa <- nu.Unif[,1]                                                   # BJ
      nu.Unifb <- nu.Unif[,2]                                                   # BJ
      if(any(nu.Unif <= 0, nu.Unifa - nu.Unifb >= 0)){                          # BJ
        stop("error: nu.Unif must be positive with element 1 < element 2")      # BJ
      }                                                                         # BJ
    }
  }
  storage.mode(nu.Unif) <- "double"
  storage.mode(nu.Unifa) = storage.mode(nu.Unifb) <- "double"                   # BJ

  if (mv.model == "separable") {                                                # BJ
    if ("rho.unif" %in% names(priors)) {                                        # BJ
      rho.Unif <- priors[["rho.unif"]]                                          # BJ
    } else {                                                                    # BJ
      rho.Unif <- c(-1, 1)                                                      # BJ
    }                                                                           # BJ
    if(!is.vector(rho.Unif) || length(rho.Unif) != 2){                          # BJ
      stop("error: rho.Unif must be a vector of length 2")                      # BJ
    }                                                                           # BJ
    storage.mode(rho.Unif) <- "double"                                          # BJ
  }                                                                             # BJ
  
  ####################################################
  ##Starting values
  ####################################################
  beta.starting <- 0
  sigma.sq.starting <- 0
  tau.sq.starting <- 0
  phi.starting <- 0
  cross.phi.starting <- rep(0, q-1)                                             # BJ
  nu.starting <- rep(0, q)                                                      # BJ
  rho.starting <- 0                                                             # BJ
  w.starting <- c(starting$w)                                                   # BJ: debug
  
  if(missing(starting)){stop("error: starting value list for the parameters must be specified")}
  
  names(starting) <- tolower(names(starting))   
  
  if(!"beta" %in% names(starting)){                                             # BJ
    beta.starting <- as.vector(coefficients(lm(as.vector(Y)~Xbd-1)))            # BJ
  }else{                                                                        # BJ
    beta.starting <- starting[["beta"]]                                         # BJ
    if (nrow(beta.starting) != p) {                                             # BJ
      stop("error: supplied beta.starting has the wrong length")                # BJ
    }                                                                           # BJ
    beta.starting <- as.vector(beta.starting)                                   # BJ
  }                                                                             # BJ 
  
  if(!"sigma.sq" %in% names(starting)){stop("error: sigma.sq must be specified in starting value list")}
  if (mv.model == "separable") {                                                # BJ 
    sigma.sq.starting <- starting[["sigma.sq"]][1]
  } else {
    sigma.sq.starting <- starting[["sigma.sq"]]                                 # BJ
    if (!is.vector(sigma.sq.starting) || length(sigma.sq.starting) != q) {      # BJ
      stop("error: sigma.sq.starting must be a vector of length q")             # BJ
    }                                                                           # BJ
  }
  
  if(family == "gaussian"){
    if(!"tau.sq" %in% names(starting)){stop("error: tau.sq must be specified in starting value list")}
    tau.sq.starting <- starting[["tau.sq"]]                                     # BJ: q vector
    if (!is.vector(tau.sq.starting) || length(tau.sq.starting) != q) {          # BJ
      stop("error: tau.sq.starting must be a vector of length q")               # BJ
    }  
  }
  
  if(!"phi" %in% names(starting)){stop("error: phi must be specified in starting value list")}
  if (mv.model == "separable") {                                                # BJ
    phi.starting <- starting[["phi"]][1]
  } else {
    phi.starting <- starting[["phi"]]                                           # BJ
    if (!is.vector(phi.starting) || length(phi.starting) != q) {                # BJ
      stop("error: phi.starting must be a vector of length q")                  # BJ
    }                                                                           # BJ
  }
  
  if(cov.model == "matern"){
    if(!"nu" %in% names(starting)){stop("error: nu must be specified in starting value list")}
    if (mv.model == "separable") {                                              # BJ
      nu.starting <- starting[["nu"]][1]                
    } else {
      nu.starting <- starting[["nu"]]                                           # BJ
      if (!is.vector(nu.starting) || length(nu.starting) != q) {                # BJ
        stop("error: nu.starting must be a vector of length q")                 # BJ
      }                                                                         # BJ
    }
  }
  
  if(!"adjmat" %in% names(starting)) {                                          # BJ
    stop("error: adjmat (graph adjacency matrix) must be specified in starting value list")} # BJ                     
  adjmat.starting <- starting[["adjmat"]]                                       # BJ: a qxq matrix
  emp_rho <- cor(Y)                                                             # BJ
  if (isSymmetric(adjmat.starting)) {                                           # BJ
    print("An undirected graph is provided. Its directed minimum spanning tree with root 1 will be used whose edge weights are proportional to the negative of absolute empirical correlation coefficients.") # BJ
    adjmat.starting <- directedMST(-abs(emp_rho)*adjmat.starting, root = 1)     # BJ
  }                                                                             # BJ
  
  if(!"rho" %in% names(starting)) {                                             # BJ
    rho.starting <- colSums(emp_rho*adjmat.starting)[-1]                        # BJ
  } else {
    rho.starting <- starting[["rho"]]                                           # BJ
  }
  if (!is.vector(rho.starting) || length(rho.starting) != q-1) {                # BJ
    stop("error: rho.starting must be a vector of length q-1")                  # BJ: without root
  }                                                                             # BJ
  
  if (mv.model == "nonseparable") {
    if(!"cross.phi" %in% names(starting)) {                                     # BJ
      for (j in 2:q) {                                                          # BJ
        parent <- which(adjmat.starting[,j] == 1)                               # BJ
        cross.phi.starting[j-1] <- sqrt(.5*phi[parent]^2+.5*phi[j]^2)           # BJ
      }                                                                         # BJ
    } else {                                                                    # BJ
      cross.phi.starting <- starting[["cross.phi"]]                             # BJ
    }
    storage.mode(cross.phi.starting) <- "double"                                # BJ
  }
  
  storage.mode(w.starting) <- "double"                                          # BJ: debug
  storage.mode(beta.starting) <- "double"
  storage.mode(sigma.sq.starting) <- "double"
  storage.mode(tau.sq.starting) <- "double"
  storage.mode(phi.starting) <- "double"
  storage.mode(nu.starting) <- "double"
  storage.mode(rho.starting) <- "double"                                        # BJ
  storage.mode(adjmat.starting) <- "integer"                                    # BJ
  
  ####################################################
  ##Tuning values
  ####################################################
  sigma.sq.tuning <- 0
  tau.sq.tuning <- 0
  phi.tuning <- 0
  nu.tuning <- rep(0, q)
  rho.tuning <- 0                                                               # BJ
  
  if(missing(tuning)){stop("error: tuning value vector for the spatial parameters must be specified")}
  
  names(tuning) <- tolower(names(tuning))
  
  if(!"sigma.sq" %in% names(tuning) & method == "response"){stop("error: sigma.sq must be specified in tuning value list")}
  if(!"sigma.sq" %in% names(tuning) & method == "latent" & mv.model == "nonseparable"){ # BJ
    stop("error: sigma.sq must be specified in tuning value list")
  }
  if (mv.model == "nonseparable") {
    sigma.sq.tuning <- tuning[["sigma.sq"]]                                     # BJ
    if (!is.vector(sigma.sq.tuning) || length(sigma.sq.tuning) != q) {          # BJ
      stop("error: sigma.sq.tuning must be a vector of length q")               # BJ
    }
  }
  if (method == "response") {                                                   # BJ
    if (mv.model == "separable") {                                              # BJ
      sigma.sq.tuning <- tuning[["sigma.sq"]][1]
    } 
  } 
  
  if(family == "gaussian"){
    if(!"tau.sq" %in% names(tuning) & method == "response"){stop("error: tau.sq must be specified in tuning value list")}
    if (method == "response") {                                                 # BJ
      tau.sq.tuning <- tuning[["tau.sq"]]                                       # BJ: q vector
      if (!is.vector(tau.sq.tuning) || length(tau.sq.tuning) != q) {            # BJ
        stop("error: tau.sq.tuning must be a vector of length q")               # BJ
      }                                                                         # BJ
    }                                                                           # BJ
  }
  
  if(!"phi" %in% names(tuning)){stop("error: phi must be specified in tuning value list")}
  if (mv.model == "separable") {                                                # BJ
    phi.tuning <- tuning[["phi"]][1]
  } else {
    phi.tuning <- tuning[["phi"]]                                               # BJ
    if (!is.vector(phi.tuning) || length(phi.tuning) != q) {                    # BJ
      stop("error: phi.tuning must be a vector of length q")                    # BJ
    }                                                                           # BJ
  }
  
  if (mv.model == "nonseparable") {
    if (!"cross.phi" %in% names(tuning)) {                                      # BJ
      stop("error: cross.phi (decay parameter for cross-covariance) must be specified in tuning value list")  # BJ
    }                                                                           # BJ
    cross.phi.tuning <- tuning[["cross.phi"]]                                   # BJ
    if (!is.vector(cross.phi.tuning) || length(cross.phi.tuning) != q-1) {      # BJ: without root
      stop("error: cross.phi.tuning must be a vector of length q-1")            # BJ
    }                                                                           # BJ
    storage.mode(cross.phi.tuning) <- "double"                                    # BJ
  }

  if(cov.model == "matern"){
    if(!"nu" %in% names(tuning)){stop("error: nu must be specified in tuning value list")}
    if (mv.model == "separable") {                                              # BJ
      nu.tuning <- tuning[["nu"]][1]
    } else {
      nu.tuning <- tuning[["nu"]]                                               # BJ
      if (!is.vector(nu.tuning) || length(nu.tuning) != q) {                    # BJ
        stop("error: nu.tuning must be a vector of length q")                   # BJ
      }                                                                         # BJ
    }
  }    
  
  if (!"rho" %in% names(tuning)) {                                              # BJ
    stop("error: rho must be specified in tuning value list")                   # BJ
  }                                                                             # BJ
  rho.tuning <- tuning[["rho"]]                                                 # BJ
  if (!is.vector(rho.tuning) || length(rho.tuning) != q-1) {                    # BJ: without root
    stop("error: rho.tuning must be a vector of length q-1")                    # BJ
  }                                                                             # BJ
  
  storage.mode(sigma.sq.tuning) <- "double"
  storage.mode(tau.sq.tuning) <- "double"
  storage.mode(phi.tuning) <- "double"
  storage.mode(nu.tuning) <- "double"
  storage.mode(rho.tuning) <- "double"                                          # BJ
  
  ####################################################
  ##nn search 
  ####################################################
  if(!neighbor.info.provided){
    
    if(verbose){
      cat("----------------------------------------\n");
      cat("\tBuilding the neighbor list\n");
      cat("----------------------------------------\n");
    }
    
    search.type.names <- c("brute", "cb")
    
    if(!search.type %in% search.type.names){
      stop("error: specified search.type '",search.type,"' is not a valid option; choose from ", paste(search.type.names, collapse=", ", sep="") ,".")
    }
    
    ##indexes
    if(search.type == "brute"){
      indx <- mkNNIndx(coords, n.neighbors, n.omp.threads)
    }else{
      indx <- mkNNIndxCB(coords, n.neighbors, n.omp.threads)
    }
    
    nn.indx <- indx$nnIndx
    nn.indx.lu <- indx$nnIndxLU
    nn.dist <- indx$nnDist                                                      # BJ
    nn.indx.run.time <- indx$run.time
    
    nn.indx.parent.run.time <- system.time({
      nn.parent.res <- nn2(coords, k = n.neighbors+1)                           # BJ
      nn.indx.parent <- c(t(nn.parent.res$nn.idx)-1)                            # BJ: m neighbors and itself, -1 for cpp index
      nn.dist.parent <- c(t(nn.parent.res$nn.dists))                            # BJ
    })                                                                          # BJ
    nn.indx.lu.parent <- c(cumsum(c(0,rep(n.neighbors+1, n-1))),   
                           rep(n.neighbors+1, n))                               # BJ
    nn.indx.lu.all <- nn.indx.lu + nn.indx.lu.parent                            # BJ
    
    storage.mode(nn.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"
    storage.mode(nn.dist) <- "double"                                           # BJ
    
    storage.mode(nn.indx.parent) <- "integer"                                   # BJ
    storage.mode(nn.indx.lu.parent) <- "integer"                                # BJ
    storage.mode(nn.dist.parent) <- "double"                                    # BJ
    storage.mode(nn.indx.lu.all) <- "integer"                                   # BJ
    
    if(method == "latent"){
      
      if(verbose){
        cat("----------------------------------------\n");
        cat("Building the neighbors of neighbors list\n");
        cat("----------------------------------------\n");
      }
      
      indx <- mkUIndx(n, n.neighbors, nn.indx, nn.indx.lu, u.search.type)
      
      u.indx <- indx$u.indx
      u.indx.lu <- indx$u.indx.lu
      ui.indx <- indx$ui.indx
      u.indx.run.time <- indx$run.time
      
      storage.mode(u.indx) <- "integer"
      storage.mode(u.indx.lu) <- "integer"
      storage.mode(ui.indx) <- "integer"
      
      uindx <- mkuUIndx(n, n.neighbors, nn.indx.parent, nn.indx.lu.parent) # BJ: test

      uu.indx <- uindx$uu.indx
      uu.indx.lu <- uindx$uu.indx.lu
      uui.indx <- uindx$uui.indx
      uu.indx.run.time <- uindx$run.time

      storage.mode(uu.indx) <- "integer"
      storage.mode(uu.indx.lu) <- "integer"
      storage.mode(uui.indx) <- "integer"
      
      c.indx <- unlist(apply(adjmat.starting, 1, function(x) which(x == 1)-1))  # BJ
      c.indx.len <- unlist(apply(adjmat.starting, 1, function(x) sum(x == 1)))  # BJ
      c.indx.lu <- c(cumsum(c(0,c.indx.len[-q])), c.indx.len)                   # BJ: 2q 
      
      storage.mode(c.indx) <- "integer"                                         # BJ
      storage.mode(c.indx.lu) <- "integer"                                      # BJ
    }
  }  
  
  ####################################################
  ##Other stuff
  ####################################################
  storage.mode(n.samples) <- "integer"
  storage.mode(n.omp.threads) <- "integer"
  storage.mode(n.neighbors) <- "integer"
  storage.mode(n.report) <- "integer"
  storage.mode(verbose) <- "integer"
  
  ####################################################
  ##Pack it up and off it goes
  ####################################################

  ptm <- proc.time()

  if(family == "gaussian"){

    if(method == "response"){
      if (mv.model == "separable") {
        out <- .Call("rSTNNGP", Y, X, q, p, n, n.neighbors, coords,
                     cov.model.indx, nn.indx, nn.indx.lu,
                     nn.indx.parent, nn.indx.lu.parent, nn.indx.lu.all,
                     sigma.sq.IG,
                     tau.sq.IGa, tau.sq.IGb,
                     phi.Unif, nu.Unif,
                     rho.Unif,
                     beta.starting, sigma.sq.starting, tau.sq.starting,
                     phi.starting, nu.starting,
                     rho.starting,
                     adjmat.starting,
                     sigma.sq.tuning, tau.sq.tuning, phi.tuning, nu.tuning,
                     rho.tuning,
                     n.samples, n.omp.threads, verbose, n.report)
      } else {
        out <- .Call("rSTNNGP_NS", Y, X, q, p, n, n.neighbors, coords,
                     cov.model.indx, nn.indx, nn.indx.lu,
                     nn.indx.parent, nn.indx.lu.parent, nn.indx.lu.all,
                     sigma.sq.IGa, sigma.sq.IGb,
                     tau.sq.IGa, tau.sq.IGb,
                     phi.Unifa, phi.Unifb,
                     nu.Unifa, nu.Unifb,
                     beta.starting, sigma.sq.starting, tau.sq.starting,
                     phi.starting,
                     cross.phi.starting,
                     nu.starting,
                     rho.starting,
                     adjmat.starting,
                     sigma.sq.tuning, tau.sq.tuning, phi.tuning, cross.phi.tuning,
                     nu.tuning, rho.tuning,
                     n.samples, n.omp.threads, verbose, n.report)
      }
    }else{##sequential
      if (mv.model == "separable") {
        if (misalign) {
          out <- .Call("sSTNNGP_misalign", Y, X, 
                       XtX,                                                     # BJ
                       q, p, n, n.neighbors, coords,
                       nj, Y_missing,                                           # BJ
                       cov.model.indx, nn.indx, nn.indx.lu,
                       nn.indx.parent, nn.indx.lu.parent, nn.indx.lu.all,
                       u.indx, u.indx.lu, ui.indx,                              # BJ: added for sSTNNGP compared to rSTNNGP
                       uu.indx, uu.indx.lu, uui.indx,                           # BJ: for case 3
                       c.indx, c.indx.lu,                                       # BJ
                       sigma.sq.IG,
                       tau.sq.IGa, tau.sq.IGb,
                       phi.Unif, nu.Unif,
                       rho.Unif,
                       beta.starting, sigma.sq.starting, tau.sq.starting,
                       phi.starting, nu.starting,
                       rho.starting,
                       adjmat.starting,
                       phi.tuning, nu.tuning,
                       rho.tuning,
                       n.samples, n.omp.threads, verbose, n.report)
        } else {
          out <- .Call("sSTNNGP", Y, X, q, p, n, n.neighbors, coords,
                       cov.model.indx, nn.indx, nn.indx.lu,
                       nn.indx.parent, nn.indx.lu.parent, nn.indx.lu.all,
                       u.indx, u.indx.lu, ui.indx,                              # BJ: added for sSTNNGP compared to rSTNNGP
                       uu.indx, uu.indx.lu, uui.indx,                           # BJ: for case 3
                       c.indx, c.indx.lu,                                       # BJ
                       sigma.sq.IG,
                       tau.sq.IGa, tau.sq.IGb,
                       phi.Unif, nu.Unif,
                       rho.Unif,
                       beta.starting, sigma.sq.starting, tau.sq.starting,
                       phi.starting, nu.starting,
                       rho.starting,
                       adjmat.starting,
                       phi.tuning, nu.tuning,
                       rho.tuning,
                       n.samples, n.omp.threads, verbose, n.report)
        }
      } else {
        if (misalign) {
          out <- .Call("sSTNNGP_NS_misalign", Y, X, 
                       XtX,                                                     # BJ
                       q, p, n, n.neighbors, coords,
                       nj, Y_missing,                                           # BJ
                       cov.model.indx, nn.indx, nn.indx.lu,
                       nn.indx.parent, nn.indx.lu.parent, nn.indx.lu.all,
                       u.indx, u.indx.lu, ui.indx,                              # BJ: added for sSTNNGP compared to rSTNNGP
                       uu.indx, uu.indx.lu, uui.indx,                           # BJ: for case 3
                       c.indx, c.indx.lu,                                       # BJ
                       sigma.sq.IGa, sigma.sq.IGb,                              # BJ: added compared to sSTNNGP
                       tau.sq.IGa, tau.sq.IGb,
                       phi.Unifa, phi.Unifb,                                    # BJ: added compared to sSTNNGP
                       nu.Unifa, nu.Unifb,                                      # BJ: added compared to sSTNNGP
                       beta.starting, sigma.sq.starting, tau.sq.starting,
                       phi.starting,
                       cross.phi.starting,                                      # BJ: added compared to sSTNNGP
                       nu.starting,
                       rho.starting,
                       adjmat.starting,
                       sigma.sq.tuning,
                       phi.tuning,
                       cross.phi.tuning,                                        # BJ: added compared to sSTNNGP
                       nu.tuning,
                       rho.tuning,
                       n.samples, n.omp.threads, verbose, n.report)
        } else {
          out <- .Call("sSTNNGP_NS", Y, X, q, p, n, n.neighbors, coords,
                       cov.model.indx, nn.indx, nn.indx.lu,
                       nn.indx.parent, nn.indx.lu.parent, nn.indx.lu.all,
                       u.indx, u.indx.lu, ui.indx,                              # BJ: added for sSTNNGP compared to rSTNNGP
                       uu.indx, uu.indx.lu, uui.indx,                           # BJ: for case 3
                       c.indx, c.indx.lu,                                       # BJ
                       sigma.sq.IGa, sigma.sq.IGb,                              # BJ: added compared to sSTNNGP
                       tau.sq.IGa, tau.sq.IGb,
                       phi.Unifa, phi.Unifb,                                    # BJ: added compared to sSTNNGP
                       nu.Unifa, nu.Unifb,                                      # BJ: added compared to sSTNNGP
                       beta.starting, sigma.sq.starting, tau.sq.starting,
                       phi.starting,
                       cross.phi.starting,                                      # BJ: added compared to sSTNNGP
                       nu.starting,
                       rho.starting,
                       adjmat.starting,
                       sigma.sq.tuning,
                       phi.tuning,
                       cross.phi.tuning,                                        # BJ: added compared to sSTNNGP
                       nu.tuning,
                       rho.tuning,
                       n.samples, n.omp.threads, verbose, n.report)
        }
      }
    }
    
    col.names <- c("sigma.sq", "phi")
    
    if(cov.model == "matern"){
      col.names <- c("sigma.sq", "phi", "nu")
    }
    
  }else{
    if (mv.model == "separable") {
      out <- .Call("sSTNNGPLogit", Y, X, q, p, n, n.neighbors, coords, 
                   Y_missing, weights, cov.model.indx, nn.indx, nn.indx.lu, 
                   nn.indx.parent, nn.indx.lu.parent, nn.indx.lu.all,
                   u.indx, u.indx.lu, ui.indx,
                   uu.indx, uu.indx.lu, uui.indx,                               # BJ: for case 3
                   c.indx, c.indx.lu,   
                   sigma.sq.IG, phi.Unif, nu.Unif, rho.Unif,
                   w.starting,                                                  # BJ: debug
                   beta.starting, sigma.sq.starting, phi.starting, nu.starting,
                   rho.starting, adjmat.starting,
                   phi.tuning, nu.tuning, rho.tuning,
                   n.samples, n.omp.threads, verbose, n.report, 
                   update.sigma.sq, update.beta, update.w)
      
      col.names <- c("sigma.sq", "phi")
      
      if(cov.model == "matern"){
        col.names <- c("sigma.sq", "phi", "nu")
      }
    } else {
      out <- .Call("sSTNNGPLogit_NS", Y, X, 
                   q, p, n, n.neighbors, coords,
                   Y_missing,                                                   # BJ
                   weights, 
                   cov.model.indx, nn.indx, nn.indx.lu,
                   nn.indx.parent, nn.indx.lu.parent, nn.indx.lu.all,
                   u.indx, u.indx.lu, ui.indx,                                  # BJ: added for sSTNNGP compared to rSTNNGP
                   uu.indx, uu.indx.lu, uui.indx,                               # BJ: for case 3
                   c.indx, c.indx.lu,                                           # BJ
                   sigma.sq.IGa, sigma.sq.IGb,                                  # BJ: added compared to sSTNNGP
                   phi.Unifa, phi.Unifb,                                        # BJ: added compared to sSTNNGP
                   nu.Unifa, nu.Unifb,                                          # BJ: added compared to sSTNNGP
                   w.starting,                                                  # BJ: debug
                   beta.starting, sigma.sq.starting, 
                   phi.starting,
                   cross.phi.starting,                                          # BJ: added compared to sSTNNGP
                   nu.starting,
                   rho.starting,
                   adjmat.starting,
                   sigma.sq.tuning,
                   phi.tuning,
                   cross.phi.tuning,                                            # BJ: added compared to sSTNNGP
                   nu.tuning,
                   rho.tuning,
                   n.samples, n.omp.threads, verbose, n.report)
    }
  }

  out$run.time <- proc.time() - ptm

  out$p.beta.samples <- mcmc(t(out$p.beta.samples))
  if (family == "gaussian") {
    out$p.tausq.samples <- mcmc(t(out$p.tausq.samples))
  }
  out$p.rho.samples <- mcmc(t(out$p.rho.samples))
  if (mv.model == "separable") {                                                # BJ
    out$p.theta.samples <- mcmc(t(out$p.theta.samples))
    colnames(out$p.theta.samples) <- col.names
  } else {
    out$p.sigmasq.samples <- mcmc(t(out$p.sigmasq.samples))                     # BJ
    out$p.phi.samples <- mcmc(t(out$p.phi.samples))                             # BJ
    out$p.crossphi.samples <- mcmc(t(out$p.crossphi.samples))                   # BJ
    if (cov.model == "matern") {                                                # BJ
      out$p.nu.samples <- mcmc(t(out$p.nu.samples))                             # BJ
    }                                                                           # BJ
  }

  if(return.neighbor.info){
    if(method == "response"){
      out$neighbor.info <- list(n.neighbors = n.neighbors,
                                n.indx = mk.n.indx.list(nn.indx, n, n.neighbors),
                                nn.indx = nn.indx, nn.indx.lu = nn.indx.lu,
                                nn.dist = nn.dist,                              # BJ
                                nn.indx.parent = nn.indx.parent,                # BJ
                                nn.indx.lu.parent = nn.indx.lu.parent,          # BJ
                                nn.dist.parent = nn.dist.parent,                # BJ
                                nn.indx.lu.all = nn.indx.lu.all,                # BJ
                                ord = ord,
                                nn.indx.run.time = nn.indx.run.time,
                                nn.indx.parent.run.time = nn.indx.parent.run.time) # BJ
    }else{
      out$neighbor.info <- list(n.neighbors = n.neighbors,
                                n.indx = mk.n.indx.list(nn.indx, n, n.neighbors),
                                nn.indx = nn.indx, nn.indx.lu = nn.indx.lu,
                                nn.dist = nn.dist,                              # BJ
                                nn.indx.parent = nn.indx.parent,                # BJ
                                nn.indx.lu.parent = nn.indx.lu.parent,          # BJ
                                nn.dist.parent = nn.dist.parent,                # BJ
                                nn.indx.lu.all = nn.indx.lu.all,                # BJ
                                u.indx = u.indx,                                # BJ: added compared to rSTNNGP
                                u.indx.lu = u.indx.lu,                          # BJ: added compared to rSTNNGP
                                ui.indx = ui.indx,                              # BJ: added compared to rSTNNGP
                                uu.indx = uu.indx,                              # BJ: added compared to rSTNNGP
                                uu.indx.lu = uu.indx.lu,                        # BJ: added compared to rSTNNGP
                                uui.indx = uui.indx,                            # BJ: added compared to rSTNNGP
                                c.indx = c.indx,                                # BJ: added compared to rSTNNGP
                                c.indx.lu = c.indx.lu,                          # BJ: added compared to rSTNNGP
                                ord = ord,
                                nn.indx.run.time = nn.indx.run.time,
                                nn.indx.parent.run.time = nn.indx.parent.run.time, # BJ
                                u.indx.run.time = u.indx.run.time,              # BJ
                                uu.indx.run.time = uu.indx.run.time)            # BJ  
    }
  }

  ##put everthing back in the original order
  out$coords <- coords[order(ord),]
  if (misalign) {
    Y[Y == -9999] <- NA
  }
  out$Y <- Y[order(ord),,drop=FALSE]
  out$X <- X[order(ord),,drop=FALSE]
  out$weights <- matrix(weights, nrow = n, ncol = q)[order(ord),,drop=FALSE]    # BJ

  if(method == "latent"){
    out$p.w.samples <- out$p.w.samples[rep(0:(q-1)*n, each = n)+order(ord),,drop=FALSE] # BJ
    if (misalign) {
      p.ypred.samples <- (diag(q)%x%out$X)%*%t(out$p.beta.samples) + 
        out$p.w.samples + 
        bind_rows(replicate(n, data.frame(sqrt(t(out$p.tausq.samples))), simplify = FALSE))*
        matrix(rnorm(n*q), nrow = n*q, ncol = n.samples)                        # BJ: posterior predictive sample
      out$p.ypred.samples <- as.matrix(p.ypred.samples)
    }
  }

  out$n.neighbors <- n.neighbors
  out$cov.model <- cov.model
  out$cov.model.indx <- cov.model.indx
  out$mv.model <- mv.model                                                      # BJ
  out$adjvec <- c(adjmat.starting)                                              # BJ
  out$starting <- starting
  out$priors <- priors
  out$tuning <- tuning
  out$type <- c(method, family)
  class(out) <- "STNNGP"

  out

}
