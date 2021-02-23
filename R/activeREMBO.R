#' REMBO for unconstrained problems, with changing embedding
#' @title Random Embedding Bayesian Optimization
#' @param par vector whose length is used to define the initial DoE or just the low dimension
#' @param fn function to be minimized over
#' @param ... additional parameters of fn
#' @param lower,upper bounds for optimization
#' @param budget total number of calls to the objective function
#' @param highDimGP should the GP use true design values?
# #' @param NashSearch should the search for the optimum be a Nash equilibrium for EI vs AS criterion?
# #' @param NashOptions list with parameters: \code{ns} the vector giving the number of strategies for EI and AS, respectively. 
# #' Note: the maximum of EI is added by default.
#' @param useAScrit should an optimization on a sequential active subspace identification criterion be performed in the orthogonal? 
#' @param kmcontrol an optional list of control parameters to be passed to the \code{\link[DiceKriging]{km}} model:
#' \code{iso}, \code{covtype}, \code{formula}. In addition, boolean \code{codereestim} is passed to \code{\link[DiceKriging]{update.km}}
#' @param control an optional list of control parameters. See "Details"
#' @param init optional list with elements \code{Amat} to provide a random matrix, \code{low_dim_design} for an initial design in the low-dimensional space and \code{fvalues} for the corresponding response.
#' When passing initial response values, care should be taken that the mapping with \code{Amat} of the design actually correspond to high-dimensional designs giving \code{fvalues}.
#' @details Options available from \code{control} are:
#' \itemize{
#' \item \code{Atype} see \code{\link[RRembo]{selectA}};
#' \item \code{reverse} if \code{TRUE}, use the new mapping from the zonotope,
#'  otherwise the original mapping with convex projection;
#' \item \code{bxsize} scalar controling the box size in the low-dimensional space;
#' \item \code{testU} with the regular mapping, set to \code{TRUE} to check that points are in U (to avoid non-injectivity);
#' \item \code{standard} for using settings of the original REMBO method;
#' \item \code{maxitOptA} if \code{Atype} is \code{optimized}, number of optimization iterations;
#' \item \code{lightreturn} only returns \code{par} and \code{value};
#' \item \code{warping} either \code{"standard"} for kY, \code{"kX"} or \code{"Psi"};
#' \item \code{designtype} one of "\code{LHS}", "\code{maximin}" and "\code{unif}",
#'  see \code{\link[RRembo]{designZ}} or \code{\link[RRembo]{designU}};
#' \item \code{tcheckP} minimal distance to an existing solution, see \code{\link[GPareto]{checkPredict}}
#' \item \code{roll} to alternate between optimization methods;
#' \item \code{inneroptim} optimization method for EI
#' \item \code{popsize, gen} population size and number of optimization generations of EI
#' }
#' @importFrom rgenoud genoud
#' @importFrom DiceOptim EI
#' @importFrom GPareto checkPredict
#' @import DiceKriging
#' @importFrom pso psoptim
#' @importFrom DEoptim DEoptim
#' @importFrom rgenoud genoud
#' @importFrom graphics axis filled.contour points title
#' @import OOR
#' @author Mickael Binois
#' @export
#' @references
#' O. Roustant, D. Ginsbourger & Y. Deville, DiceKriging, DiceOptim: Two R Packages for the Analysis of Computer Experiments by Kriging-Based Metamodeling and Optimization Journal of Statistical Software, 2012, 51, 1-55\cr \cr
#' Z. Wang, F. Hutter, M. Zoghi, D. Matheson, N. de Freitas (2016), Bayesian Optimization in a Billion Dimensions via Random Embeddings, JAIR. \cr \cr
#' M. Binois, D. Ginsbourger, O. Roustant (2015), A Warped Kernel Improving Robustness in Bayesian Optimization Via Random Embeddings, Learning and Intelligent Optimization, Springer \cr \cr
#' M. Binois, D. Ginsbourger, O. Roustant (2018), On the choice of the low-dimensional domain for global optimization via random embeddings, arXiv:1704.05318 \cr \cr
#' M. Binois (2015), Uncertainty quantification on Pareto fronts and high-dimensional strategies in Bayesian optimization, with applications in multi-objective automotive design, PhD thesis, Mines Saint-Etienne.
#' @examples
#' \dontrun{
#' set.seed(42)
#' library(rgl)
#' library(DiceKriging)
#'
#' lowd <- 2
#' highD <- 25
#'
#' ntest <- 50
#' maxEval <- 100
#'
#' res <- rep(0, ntest)
#' for(i in 1:ntest){
#'   #ii <- sample(1:highD, 2)
#'   ii <- c(1,2)
#'   branin_mod2 <- function(X){
#'     if(is.null(nrow(X))) X <- matrix(X, nrow = 1)
#'     X <- X[, c(ii[1], ii[2]), drop = FALSE]
#'     return(apply(X, 1, branin))
#'   }
#'   sol <- activeREMBO(par = rep(NA, lowd), branin_mod2, lower = rep(0, highD),
#'                    upper = rep(1, highD), budget = maxEval, highDimGP = TRUE)
#'   res[i] <- sol$value
#'   cat(sol$value, " ", i, "\n")
#' }
#'
#' plot(res - 0.397887, type = "b")
#' boxplot(res - 0.397887)
#' }
activeREMBO <- function(par, fn, lower, upper, budget, ..., highDimGP, useAScrit = FALSE, #NashSearch = FALSE, NashOptions = list(ns = c(50, 2)), 
                        homcontrol = list(beta0 = 0, covtype = "Matern5_2"),
                        kmcontrol = list(covtype = "matern5_2", iso = TRUE, covreestim = TRUE, formula =~1),
                        control = list(Atype = 'isotropic', reverse = TRUE, bxsize = NULL, testU = TRUE, standard = FALSE,
                                       maxitOptA = 100, lightreturn = FALSE, warping = 'Psi', designtype = 'unif',
                                       tcheckP = 1e-3, roll = F,
                                       inneroptim = "pso", popsize = 80, gen = 40),
                        init = NULL){
  # Initialisation
  d <- length(par)
  D <- length(lower)
  
  # if(is.null(control$Atype)) control$Atype <- 'isotropic'
  # if(is.null(control$Ysetting)) control$Ysetting <- 'max'
  # if(is.null(control$testU)) control$testU <- TRUE
  # if(is.null(control$standard)) control$standard <- FALSE
  # if(is.null(control$maxitOptA)) control$maxitOptA <- 100
  if(is.null(control$lightreturn)) control$lightreturn <- FALSE
  if(is.null(control$warping)) control$warping <- 'Psi'
  if(is.null(control$inneroptim)) control$inneroptim <- 'pso'
  if(is.null(control$popsize)) control$popsize <- 80
  if(is.null(control$gen)) control$gen <- 40
  # if(is.null(control$designtype)) control$designtype <- 'unif'
  if(is.null(control$reverse)) control$reverse <- TRUE
  if(is.null(control$maxf)) control$maxf <- control$popsize * control$gen
  if(is.null(control$tcheckP)) control$tcheckP <- 1e-3
  if(is.null(control$roll)) control$roll <- FALSE
  
  if(is.null(kmcontrol$covtype)) kmcontrol$covtype <- 'matern5_2'
  if(is.null(kmcontrol$iso)) kmcontrol$iso <- TRUE
  if(is.null(kmcontrol$covreestim)) kmcontrol$covreestim <- TRUE
  if(is.null(kmcontrol$formula)) kmcontrol$formula <- ~1
  
  # Add for activeREMBO
  if(is.null(homcontrol$known)) homcontrol$known <- list(beta0 = 0, g = 1e-6)
  if(is.null(homcontrol$covtype)) homcontrol$covtype <- "Matern5_2"
  if(is.null(homcontrol$lower)) homcontrol$lower <- rep(0.05, D)
  if(is.null(homcontrol$upper)) homcontrol$upper <- rep(2, D)
  if(is.null(homcontrol$init)) homcontrol$init <- list(theta = rep(0.5, D))
  if(highDimGP) nuggetestim <- FALSE else nuggetestim <- TRUE
  
  
  # Number of points of the DoE
  if(is.null(init$n)){
    if(!is.null(init$low_dim_design)){
      n.init <- 0 # For update error consistency
    }else{
      n.init <- max(4 * d, round(budget/3))
    }
  }else{
    n.init <- init$n  
  }
  
  DoE <- maximinLHS(n = n.init, k = D)
  
  fvalues <- apply(DoE, 1, fn, ...)
  
  library(hetGP)
  library(activegp)
  
  
  modelD <- mleHomGP(X = DoE, Z = fvalues, lower = homcontrol$lower, upper = homcontrol$upper,
                     known = homcontrol$known, covtype = homcontrol$covtype,
                     init = homcontrol$init)
  
  C_hat <- C_GP(modelD)
  A_hat <- eigen(C_hat$mat)$vectors[, 1:d, drop = F]
  tA <- t(A_hat)
  Amat <- cbind(A_hat, matrix(rep(c(1, 0), times = c(1, D-1)), D, D), matrix(rep(c(-1, 0), times = c(1, D-1)), D, D))
  Aind <- cbind(matrix(c(D, 1:D), D + 1, d), rbind(rep(1, D * 2), c(1:D, 1:D), matrix(0, D-1, D*2)))
  
  # Initial mapping
  if(control$warping == 'kX'){
    map <- mapZX
    formals(map)$Amat <- Amat
    formals(map)$Aind <- Aind
  }
  if(control$warping == 'Psi'){
    map <- Psi_Z
    formals(map)$Amat <- Amat
    formals(map)$Aind <- Aind
  }
  
  # now build design for Rembo
  design <- (DoE * 2 - 1) # rescale to [-1,1]^D
  design <- ortProj(design, tA) # now design is in Z
  
  # best value so far, slightly perturbed, to increase exploitation in EI optimization
  boundsEIopt <- rowSums(abs(tA)) ## box surrounding the zonotope Z
  spartan <- pmin(boundsEIopt, pmax(-boundsEIopt, 
                                    design[which.min(fvalues),] + rnorm(d, sd = 0.05)))
  
  
  # change design for GP fitting
  if(!highDimGP) design <- map(design, A_hat) else design <- (DoE * 2 - 1)
  
  model <- km(kmcontrol$formula, design = design, response = fvalues, covtype = kmcontrol$covtype, 
              iso = kmcontrol$iso, nugget.estim = nuggetestim,
              control = list(trace = FALSE))
  
  
  while(model@n < budget){
    
    # Expected Improvement function in the REMBO case
    # Not evaluated and penalized if not in U/Z or
    # mval: max penalty value (negative value)
    
    plugin <- min(predict(model, newdata = model@X, type = "UK")$mean)
    
    # TODO: issue at the border with map(x) not equal to x  
    EI_Rembo <- function(x, model, mval = -10){
      if(is.null(nrow(x)))
        x <- matrix(x, nrow = 1)
      
      inDomain <- rep(TRUE, nrow(x))
      
      inDomain <- testZ(x, tA)
      
      res <- rep(NA, nrow(x))
      if(any(inDomain)){
        
        if(highDimGP) xtmp <- mapZX(x[inDomain,], A_hat, Amat = Amat, Aind = Aind) else xtmp <- map(x[inDomain,], A_hat)
        
  
        # identify too close points
        tmp <- checkPredict(xtmp, list(model), control$tcheckP, distance = 'euclidean')
        
        res[inDomain[tmp]] <- 0
        if(any(!tmp)) res[inDomain[!tmp]] <- EI(x = xtmp[!tmp,], model = model, plugin = plugin)
      }
      
      if(any(!inDomain)) res[!inDomain] <- mval * apply(x[!inDomain,,drop = FALSE], 1, distance, x2 = 0)
      return(res)
    }
    
    boundsEIopt <- rowSums(abs(tA)) ## box surrounding the zonotope Z
    
    
    # Change optimization method every time
    if(control$roll){
      if(model@n %% 4 == 0) control$inneroptim <- 'genoud'
      if(model@n %% 4 == 1) control$inneroptim <- 'DEoptim'
      if(model@n %% 4 == 2) control$inneroptim <- 'genoud'
      if(model@n %% 4 == 3) control$inneroptim <- 'DEoptim'
    }
    
    if(control$inneroptim == 'DEoptim'){
      EI_2_DE <- function(x, model, mval = -10){
        return(-EI_Rembo(x, model, mval))
      }
      parinit <- 2*matrix(runif(d*control$popsize), control$popsize) - 1
      parinit[1,] <- spartan
      opt <- DEoptim(fn = EI_2_DE, lower = -boundsEIopt, upper = boundsEIopt, model = model,
                     control = list(initialpop = parinit, NP = control$popsize, itermax = control$gen, trace = F))
      opt$value <- -opt$optim$bestval
      opt$par <- opt$optim$bestmem
    }
    
    if(control$inneroptim == 'pso'){
      # parinit <- 2*matrix(runif(d), 1) - 1
      parinit <- spartan
      opt <- psoptim(parinit, fn = EI_Rembo, lower = -boundsEIopt, upper = boundsEIopt,
                     control = list(fnscale = -1, maxit = control$gen, s = control$popsize,
                                    vectorize = T), model = model)
      opt$value <- -opt$value
    }
    
    if(control$inneroptim == "genoud"){
      parinit <- 2*matrix(runif(d*control$popsize), control$popsize) - 1
      parinit[1,] <- spartan
      opt <- genoud(EI_Rembo, nvars = d, max = TRUE, starting.values = parinit,
                    Domains = cbind(-boundsEIopt, boundsEIopt), print.level = 0,
                    model = model,
                    boundary.enforcement = 2, pop.size = control$popsize,
                    max.generations = control$gen)
    }
    
    if(control$inneroptim == 'StoSOO'){
      opt <- StoSOO(par = rep(NA, d), fn = EI_Rembo, lower = -boundsEIopt, upper = boundsEIopt,
                    settings = list(max = TRUE, type = 'det'), nb_iter = control$maxf, model = model)
    }
    
    if(opt$value <= 0){ # no improvement found or not even in domain
      cat("No improvement or not even in Z \n")
    }
    
    # if(NashSearch){
    #   
    #   if(max(abs((range(t(A) %*% A - diag(ncol(A)))))) > 1e-10) print("A is expected to have orthonormal columns but it does not seem to be the case here")
    #   
    #   ## Define set of strategies with territory splitting: EI is given variable on the AS, AS is given the orthogonal
    #   X_EI <- matrix(runif(NashOptions$ns[1] * D), ncol = D) * 2 - 1
    #   
    #   
    #   
    #   # Basis for the orthogonal complement of A (i.e., W)
    #   W <- orthonormalization(A, basis = TRUE,norm = TRUE)[,-c(1:d),drop = F]
    #   # X_AS <- designZ(p = NashOptions$ns[2], pA = t(W), bxsize = sqrt(D), type = "unif")
    #   X_AS <- matrix(runif(NashOptions$ns[2] * D), ncol = D) * 2 - 1
    #   
    #   
    #   ## Add minimum of EI in alternatives
    #   X_EI <- rbind(X_EI, ((mapZX(opt$par, A_hat, Amat = Amat, Aind = Aind) + 1)/2) %*% diag(upper - lower) + lower)
    #   X_AS <- rbind(X_AS, ((mapZX(opt$par, A_hat, Amat = Amat, Aind = Aind) + 1)/2) %*% diag(upper - lower) + lower)
    #   
    #   ## Get coordinates in original domain [-1,1]^D
    #   expanded.indices <- expand.grid(seq(1, NashOptions$ns[1]), seq(1, NashOptions$ns[2]))
    #   Xn <- cbind((X_EI %*% A)[expanded.indices[,1],, drop = F], (X_AS %*% W)[expanded.indices[,2],, drop = F])
    #   Xn <- Xn %*% t(cbind(A, W))
    #   
    #   # plot3d(Xn)
    #   # points3d((matrix(runif(2*1000), 1000) * 2 - 1) %*% t(A), col = "red")
    #   
    #   ## Identify points not outside of the domain
    #   ids_in <- which(rowSums(apply(Xn, c(1,2), function(x) max(abs(x) > 1))) > 0)
    #   
    # }
    
    newX <- ((mapZX(opt$par, A_hat, Amat = Amat, Aind = Aind) + 1)/2) %*% diag(upper - lower) + lower
    
    if(useAScrit){
      
      W_hat <- orthonormalization(A_hat, basis = TRUE,norm = TRUE)[,-c(1:d),drop = F]
      yEI <- (newX * 2 - 1) %*% A_hat # solution in R^d of EI optimization
      
      # plot3d(((matrix((runif(1000*2) * 4 - 2),1000) %*% t(A_hat))))
      # XXX <- NULL
      
      #@param w component in the orthogonal space of A (dimension D - d)
      #@param yEI solution of EI in [-1,1]^D project on ran(A) (dimension d)
      #@param C AS matrix (see activegp)
      #@param mval penalty to go back
      af_ort <- function(w, yEI, C, A, W, mval = -10){
        if(is.null(nrow(w)))
          w <- matrix(w, nrow = 1)
        
        # recreate full x vector
        x <- cbind(matrix(yEI, nrow = nrow(w), ncol = length(yEI), byrow = TRUE), w) %*% t(cbind(A, W))
        # XXX <<- rbind(XXX, x)
        x <- ((x + 1)/2)  %*% diag(upper - lower) + lower
        
        inDomain <- (rowSums(apply(x, c(1,2), function(x) max(abs(x) > 1))) > 0)
        
        res <- rep(NA, nrow(x))
        if(any(inDomain)){
          res[inDomain] <- C_var(C, x[inDomain,,drop=F], grad = FALSE)
        }
        
        if(any(!inDomain)) res[!inDomain] <- mval * apply(x[!inDomain,,drop = FALSE], 1, distance, x2 = 0)
        
        return(res)
      }
      
      opt_af <- psoptim(rep(NA, D-d), fn = af_ort, lower = rep(-sqrt(D), D - d), upper = rep(sqrt(D), D - d),
                        control = list(fnscale = -1, maxit = control$gen, s = control$popsize,
                                       vectorize = T), C = C_hat, A = A_hat, W = W_hat, yEI = yEI)
      
      newX2 <- cbind(yEI, opt_af$par) %*% t(cbind(A_hat, W_hat))
      newX2 <- pmin(1, pmax(-1, newX2))
      newX2 <- ((newX2 + 1)/2) %*% diag(upper - lower) + lower
      
      # Verif: the projection of newX2 should match the one of newX
      # (newX * 2 - 1) %*% A_hat
      # (newX2 * 2 - 1) %*% A_hat
      
      newX <- newX2
    }
    
    
    newY <- fn(newX, ...)
    
    fvalues <- c(fvalues, newY)
    DoE <- rbind(DoE, newX)
    
    modelD <- mleHomGP(X = DoE, Z = fvalues, lower = homcontrol$lower, upper = homcontrol$upper,
                       known = homcontrol$known, covtype = homcontrol$covtype, init = homcontrol$init)
    
    C_hat <- C_GP(modelD)
    A_hat <- eigen(C_hat$mat)$vectors[, 1:d, drop = F]
    
    # Adapt to change of A matrix
    ## Reverse mode only
    ## precomputations
    tA <- t(A_hat)
    Amat <- cbind(A_hat, matrix(rep(c(1, 0), times = c(1, D-1)), D, D), matrix(rep(c(-1, 0), times = c(1, D-1)), D, D))
    Aind <- cbind(matrix(c(D, 1:D), D + 1, d), rbind(rep(1, D * 2), c(1:D, 1:D), matrix(0, D-1, D*2)))
    
    if(control$warping == 'kX'){
      map <- mapZX
      formals(map)$Amat <- Amat
      formals(map)$Aind <- Aind
    }
    if(control$warping == 'Psi'){
      map <- Psi_Z
      formals(map)$Amat <- Amat
      formals(map)$Aind <- Aind
    }
    ##
    # now build design for Rembo
    design <- (DoE * 2 - 1) # rescale to [-1,1]^D
    design <- ortProj(design, tA)
    
    
    ind <- which.min(predict(model, newdata = model@X, type = "UK")$mean)
    spartan <- pmin(boundsEIopt, pmax(-boundsEIopt, design[ind,] + rnorm(d, sd = 0.05)))
    
    if(!highDimGP) design <- map(design, A_hat) else design <- DoE * 2 - 1
    
    minpar <- DoE[which.min(fvalues),]
    
    # newmodel <- try(update(model, newDesign, newY, cov.reestim = kmcontrol$covreestim,
    #                        kmcontrol = list(control = list(trace = FALSE))))
    
    newmodel <- try(km(kmcontrol$formula, design = design, response = fvalues, covtype = kmcontrol$covtype, 
                       iso = kmcontrol$iso, nugget.estim = nuggetestim,
                       control = list(trace = FALSE)))
    
    # if(typeof(newmodel) == "character"){
    #   cat("Error in hyperparameter estimation - old hyperparameter values used instead for model at iteration", model@n - n.init, "\n")
    #   newmodel <- try(update(object = model, newX = newDesign, newy = newY, newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
    # }
    
    
    if(model@n %% 10 == 0){
      print(model@n)
      plot(modelD)
    }
    if(model@n %% 10 == 2 && d == 2){
      X_grid <- as.matrix(expand.grid(seq(-boundsEIopt[1], boundsEIopt[1], length.out = 101),
                                      seq(-boundsEIopt[2], boundsEIopt[2], length.out = 101)))
      EI_grid <- apply(X_grid, 1, EI_Rembo, model = model)
      filled.contour(seq(-boundsEIopt[1], boundsEIopt[1], length.out = 101),
                     seq(-boundsEIopt[2], boundsEIopt[2], length.out = 101),
                     matrix(EI_grid, 101),
                     plot.axes = { axis(1); axis(2); points(opt$par[1], opt$par[2], col = "blue", pch = 20) })
    }
    if(model@n %% 10 == 5){
      print(model@n)
      plot(model)
    }
    if(model@n %% 10 == 7){
      plot(DoE[,1:2])
      points(DoE[which.min(fvalues), 1:2, drop = F], col = "red", pch = 20)
      points(DoE[which.min(predict(model, newdata = model@X, type = "UK")$mean), 1:2, drop = F], col = "orange", pch = 1)
    }
    
    if(model@n %% 10 == 9 && d == 2){
      X_grid <- as.matrix(expand.grid(seq(-boundsEIopt[1], boundsEIopt[1], length.out = 101),
                                      seq(-boundsEIopt[2], boundsEIopt[2], length.out = 101)))
      if(!highDimGP){
        indin2 <- apply(X_grid, 1, testZ, pA = tA)
        p_grid <- rep(NA, nrow(X_grid))
        p_grid[indin2] <- predict(model, map(X_grid[indin2,], A_hat), type = "UK", checkNames = F)$mean
        
        filled.contour(seq(-boundsEIopt[1], boundsEIopt[1], length.out = 101),
                       seq(-boundsEIopt[2], boundsEIopt[2], length.out = 101),
                       matrix(p_grid, 101),
                       plot.axes = { axis(1); axis(2); points(opt$par, col = "blue", pch = 20) })
      } 
    }
    
    if (typeof(newmodel) == "character") {
      cat("Unable to udpate kriging model at iteration", model@n - n.init, "- optimization stopped \n")
      cat("lastmodel is the model at iteration", model@n - n.init -1, "\n")
      if(control$lightreturn){
        res <- list(par = minpar, value = min(fvalues))
      }else{
        res <- list(par = minpar,
                    value = min(fvalues),
                    low_dim_design = DoE,
                    y = fvalues,
                    A = A_hat,
                    model = model)
      }
      return(res)
    }
    model <- newmodel
    print(min(fvalues))
  }
  
  if(control$lightreturn){
    res <- list(par = minpar, value = min(fvalues))
  }else{
    res <- list(par = minpar,
                value = min(fvalues),
                low_dim_design = DoE,
                y = fvalues,
                A = A_hat,
                model = model,
                homMod = modelD)
  }
  
  return(res)
  
}

#' REMBO for unconstrained problems, based on multi-objective way
#' @title Random Embedding Bayesian Optimization
#' @param par vector whose length is used to define the initial DoE or just the low dimension
#' @param fn function to be minimized over
#' @param ... additional parameters of fn
#' @param lower,upper bounds for optimization
#' @param budget total number of calls to the objective function
#' @param highDimGP should the GP use true design values?
#' @param kmcontrol an optional list of control parameters to be passed to the \code{\link[DiceKriging]{km}} model:
#' \code{iso}, \code{covtype}, \code{formula}. In addition, boolean \code{codereestim} is passed to \code{\link[DiceKriging]{update.km}}
#' @param control an optional list of control parameters. See "Details"
#' @param init optional list with elements \code{Amat} to provide a random matrix, \code{low_dim_design} for an initial design in the low-dimensional space and \code{fvalues} for the corresponding response.
#' When passing initial response values, care should be taken that the mapping with \code{Amat} of the design actually correspond to high-dimensional designs giving \code{fvalues}.
#' @details Options available from \code{control} are:
#' \itemize{
#' \item \code{Atype} see \code{\link[RRembo]{selectA}};
#' \item \code{reverse} if \code{TRUE}, use the new mapping from the zonotope,
#'  otherwise the original mapping with convex projection;
#' \item \code{bxsize} scalar controling the box size in the low-dimensional space;
#' \item \code{testU} with the regular mapping, set to \code{TRUE} to check that points are in U (to avoid non-injectivity);
#' \item \code{standard} for using settings of the original REMBO method;
#' \item \code{maxitOptA} if \code{Atype} is \code{optimized}, number of optimization iterations;
#' \item \code{lightreturn} only returns \code{par} and \code{value};
#' \item \code{warping} either \code{"standard"} for kY, \code{"kX"} or \code{"Psi"};
#' \item \code{designtype} one of "\code{LHS}", "\code{maximin}" and "\code{unif}",
#'  see \code{\link[RRembo]{designZ}} or \code{\link[RRembo]{designU}};
#' \item \code{tcheckP} minimal distance to an existing solution, see \code{\link[GPareto]{checkPredict}}
#' \item \code{roll} to alternate between optimization methods;
#' \item \code{inneroptim} optimization method for EI
#' \item \code{popsize, gen} population size and number of optimization generations of EI
#' }
#' @importFrom rgenoud genoud
#' @importFrom DiceOptim EI
#' @importFrom GPareto checkPredict
#' @import DiceKriging
#' @importFrom pso psoptim
#' @importFrom DEoptim DEoptim
#' @importFrom rgenoud genoud
#' @import OOR
#' @author Mickael Binois
#' @export
#' @references
#' O. Roustant, D. Ginsbourger & Y. Deville, DiceKriging, DiceOptim: Two R Packages for the Analysis of Computer Experiments by Kriging-Based Metamodeling and Optimization Journal of Statistical Software, 2012, 51, 1-55\cr \cr
#' Z. Wang, F. Hutter, M. Zoghi, D. Matheson, N. de Freitas (2016), Bayesian Optimization in a Billion Dimensions via Random Embeddings, JAIR. \cr \cr
#' M. Binois, D. Ginsbourger, O. Roustant (2015), A Warped Kernel Improving Robustness in Bayesian Optimization Via Random Embeddings, Learning and Intelligent Optimization, Springer \cr \cr
#' M. Binois, D. Ginsbourger, O. Roustant (2018), On the choice of the low-dimensional domain for global optimization via random embeddings, arXiv:1704.05318 \cr \cr
#' M. Binois (2015), Uncertainty quantification on Pareto fronts and high-dimensional strategies in Bayesian optimization, with applications in multi-objective automotive design, PhD thesis, Mines Saint-Etienne.
#' @examples
#' \dontrun{
#' set.seed(42)
#' library(rgl)
#' library(DiceKriging)
#'
#' lowd <- 2
#' highD <- 25
#'
#' ntest <- 50
#' maxEval <- 100
#'
#' res <- rep(0, ntest)
#' for(i in 1:ntest){
#'   #ii <- sample(1:highD, 2)
#'   ii <- c(1,2)
#'   branin_mod2 <- function(X){
#'     if(is.null(nrow(X))) X <- matrix(X, nrow = 1)
#'     X <- X[, c(ii[1], ii[2]), drop = FALSE]
#'     return(apply(X, 1, branin))
#'   }
#'   sol <- MactiveREMBO(par = rep(NA, lowd), branin_mod2, lower = rep(0, highD),
#'                    upper = rep(1, highD), budget = maxEval, highDimGP = TRUE)
#'   res[i] <- sol$value
#'   cat(sol$value, " ", i, "\n")
#' }
#'
#' plot(res - 0.397887, type = "b")
#' boxplot(res - 0.397887)
#' }
MactiveREMBO <- function(par, fn, lower, upper, budget, ..., highDimGP, batch_size = 5,
                         homcontrol = list(beta0 = 0, covtype = "Matern5_2"),
                         kmcontrol = list(covtype = "matern5_2", iso = TRUE, covreestim = TRUE, formula =~1),
                         control = list(Atype = 'isotropic', reverse = TRUE, bxsize = NULL, testU = TRUE, standard = FALSE,
                                        maxitOptA = 100, lightreturn = FALSE, warping = 'Psi', designtype = 'unif',
                                        tcheckP = 1e-5, roll = F, plot = TRUE,
                                        inneroptim = "pso", popsize = 80, gen = 40),
                         init = NULL){
  # Initialisation
  d <- length(par)
  D <- length(lower)
  
  # if(is.null(control$Atype)) control$Atype <- 'isotropic'
  # if(is.null(control$Ysetting)) control$Ysetting <- 'max'
  # if(is.null(control$testU)) control$testU <- TRUE
  # if(is.null(control$standard)) control$standard <- FALSE
  # if(is.null(control$maxitOptA)) control$maxitOptA <- 100
  if(is.null(control$lightreturn)) control$lightreturn <- FALSE
  if(is.null(control$warping)) control$warping <- 'Psi'
  if(is.null(control$inneroptim)) control$inneroptim <- 'pso'
  if(is.null(control$popsize)) control$popsize <- 80
  if(is.null(control$gen)) control$gen <- 40
  # if(is.null(control$designtype)) control$designtype <- 'unif'
  if(is.null(control$reverse)) control$reverse <- TRUE
  if(is.null(control$maxf)) control$maxf <- control$popsize * control$gen
  if(is.null(control$tcheckP)) control$tcheckP <- 1e-4
  if(is.null(control$roll)) control$roll <- FALSE
  if(is.null(control$plot)) control$plot <- TRUE
  
  if(is.null(kmcontrol$covtype)) kmcontrol$covtype <- 'matern5_2'
  if(is.null(kmcontrol$iso)) kmcontrol$iso <- TRUE
  if(is.null(kmcontrol$covreestim)) kmcontrol$covreestim <- TRUE
  if(is.null(kmcontrol$formula)) kmcontrol$formula <- ~1
  
  # Add for activeREMBO
  if(is.null(homcontrol$known)) homcontrol$known <- list(beta0 = 0, g = 1e-6)
  if(is.null(homcontrol$covtype)) homcontrol$covtype <- "Matern5_2"
  if(is.null(homcontrol$lower)) homcontrol$lower <- rep(0.05, D)
  if(is.null(homcontrol$upper)) homcontrol$upper <- rep(2, D)
  if(is.null(homcontrol$init)) homcontrol$init <- list(theta = rep(0.5, D))
  if(highDimGP) nuggetestim <- FALSE else nuggetestim <- TRUE
  
  
  # Number of points of the DoE
  if(is.null(init$n)){
    if(!is.null(init$low_dim_design)){
      n.init <- 0 # For update error consistency
    }else{
      n.init <- max(4 * d, round(budget/3))
    }
  }else{
    n.init <- init$n  
  }
  
  DoE <- maximinLHS(n = n.init, k = D)
  
  fvalues <- apply(DoE, 1, fn, ...)
  
  library(hetGP)
  library(activegp)
  
  
  modelD <- mleHomGP(X = DoE, Z = fvalues, lower = homcontrol$lower, upper = homcontrol$upper,
                     known = homcontrol$known, covtype = homcontrol$covtype,
                     init = homcontrol$init)
  
  C_hat <- C_GP(modelD)
  A_hat <- eigen(C_hat$mat)$vectors[, 1:d, drop = F]
  tA <- t(A_hat)
  Amat <- cbind(A_hat, matrix(rep(c(1, 0), times = c(1, D-1)), D, D), matrix(rep(c(-1, 0), times = c(1, D-1)), D, D))
  Aind <- cbind(matrix(c(D, 1:D), D + 1, d), rbind(rep(1, D * 2), c(1:D, 1:D), matrix(0, D-1, D*2)))
  
  # Initial mapping
  if(control$warping == 'kY'){
    map <- function(z, A){
      if(is.null(nrow(z)))
        z <- matrix(z, nrow = 1)
      return(z)
    }
  }
  if(control$warping == 'kX'){
    map <- mapZX
    formals(map)$Amat <- Amat
    formals(map)$Aind <- Aind
  }
  if(control$warping == 'Psi'){
    map <- Psi_Z
    formals(map)$Amat <- Amat
    formals(map)$Aind <- Aind
  }
  
  # now build design for Rembo
  design <- (DoE * 2 - 1) # rescale to [-1,1]^D
  design <- ortProj(design, tA) # now design is in Z
  
  # best value so far, slightly perturbed, to increase exploitation in EI optimization
  boundsEIopt <- rowSums(abs(tA)) ## box surrounding the zonotope Z
  
  
  # change design for GP fitting
  if(!highDimGP) design <- map(design, A_hat) else design <- (DoE * 2 - 1)
  
  model <- km(kmcontrol$formula, design = design, response = fvalues, covtype = kmcontrol$covtype, 
              iso = kmcontrol$iso, nugget.estim = nuggetestim,
              control = list(trace = FALSE))
  
  
  while(model@n < budget){
    
    # Expected Improvement function in the REMBO case
    # Not evaluated and penalized if not in U/Z or
    # mval: max penalty value (negative value)
    
    plugin <- min(predict(model, newdata = model@X, type = "UK")$mean)
    
    # TODO: issue at the border with map(x) not equal to x  
    EI_Rembo <- function(x, model, mval = -10){
      if(is.null(nrow(x)))
        x <- matrix(x, nrow = 1)
      
      inDomain <- rep(TRUE, nrow(x))
      
      inDomain <- testZ(x, tA)
      
      res <- rep(NA, nrow(x))
      if(any(inDomain)){
        
        if(highDimGP) xtmp <- mapZX(x[inDomain,], A_hat, Amat = Amat, Aind = Aind) else xtmp <- map(x[inDomain,], A_hat)
        
        # identify too close points
        tmp <- checkPredict(xtmp, list(model), control$tcheckP, distance = 'euclidean') 
        
        res[inDomain[tmp]] <- 0
        if(any(!tmp)) res[inDomain[!tmp]] <- EI(x = xtmp[!tmp,], model = model, plugin = plugin)
      }
      
      if(any(!inDomain)) res[!inDomain] <- mval * apply(x[!inDomain,,drop = FALSE], 1, distance, x2 = 0)
      return(res)
    }
    
    boundsEIopt <- rowSums(abs(tA)) ## box surrounding the zonotope Z
    
    
    ## Nsga optimization of bi-objective 
    library(mco)
    
    # Second objective: maximizer la variance
    af <- function(x, C, mval = -10){
      if(is.null(nrow(x)))
        x <- matrix(x, nrow = 1)
      
      inDomain <- rep(TRUE, nrow(x))
      
      inDomain <- testZ(x, tA)
      
      res <- rep(NA, nrow(x))
      if(any(inDomain)){
        xtmp <- (mapZX(x[inDomain,, drop = F], A_hat, Amat = Amat, Aind = Aind) + 1)/2
        res[inDomain] <- C_var(C, xtmp, grad = FALSE)
      }
      
      if(any(!inDomain)) res[!inDomain] <- mval * apply(x[!inDomain,,drop = FALSE], 1, distance, x2 = 0)
      
      return(res)
    }
    
    biobjfn <- function(x, C, model, mval = -10){
      f1 <- -EI_Rembo(x = x, model = model, mval = mval)
      f2 <- af(x, C)
      if(f2 > 0) f2 <- sqrt(f2)
      f2 <- -f2
      return(c(f1, f2))
    }
    
    opt <- nsga2(fn = biobjfn, idim = d, odim = 2, lower.bounds = -boundsEIopt, upper.bounds = boundsEIopt, 
                 popsize = control$popsize, generations = control$gen,
                 model = model, C = C_hat)
    
    ## Select subset of solutions
    PF <- normalizeFront(opt$value[opt$pareto.optimal,, drop = FALSE])
    best_batch <- 1:min(batch_size, nrow(PF))
    refPoint <- apply(PF, 2, max) + 1
    best_batch_val <- dominatedHypervolume(PF[best_batch,,drop = F], refPoint)
    
    if(nrow(PF) > batch_size){
      for(i in 1:10000){
        tmp_batch <- sample(1:nrow(PF), batch_size)
        tmp_val <- dominatedHypervolume(PF[tmp_batch,,drop = F], refPoint)
        if(tmp_val > best_batch_val){
          best_batch <- tmp_batch
          best_batch_val <- tmp_val
        }
      }
    }
    
    if(control$plot){
      # plot(model)
      if(d == 2){
        par(mfrow = c(1,2))
        plot(opt$par, xlim  = c(-boundsEIopt[1],boundsEIopt[1]), 
             ylim = c(-boundsEIopt[2], boundsEIopt[2]),
             xlab = expression(z[1]), ylab = expression(z[2]))
        points(opt$par[which(opt$pareto.optimal)[best_batch],,drop = F], pch = 20, col = "red")
        points(ortProj(DoE*2-1, tA), col = "blue", pch = 17)
      }
      
      plot(opt$value, xlab = "-Expected Improvement", ylab = "-C variance")
      points(opt$value[which(opt$pareto.optimal)[best_batch],,drop = F], pch = 20, col = "red")
      if(d == 2) par(mfrow = c(1,1))
      title(paste("n: ", model@n))
    }
    
    newX <- ((mapZX(opt$par[which(opt$pareto.optimal)[best_batch],,drop = F],
                    A_hat, Amat = Amat, Aind = Aind) + 1)/2) %*% diag(upper - lower) + matrix(lower, nrow = length(best_batch), ncol = D)
    newY <- apply(newX, 1, fn, ...)
    
    fvalues <- c(fvalues, newY)
    DoE <- rbind(DoE, newX)
    
    modelD <- mleHomGP(X = DoE, Z = fvalues, lower = homcontrol$lower, upper = homcontrol$upper,
                       known = homcontrol$known, covtype = homcontrol$covtype, init = homcontrol$init)
    
    C_hat <- C_GP(modelD)
    A_hat <- eigen(C_hat$mat)$vectors[, 1:d, drop = F]
    
    # Adapt to change of A matrix
    ## Reverse mode only
    ## precomputations
    tA <- t(A_hat)
    Amat <- cbind(A_hat, matrix(rep(c(1, 0), times = c(1, D-1)), D, D), matrix(rep(c(-1, 0), times = c(1, D-1)), D, D))
    Aind <- cbind(matrix(c(D, 1:D), D + 1, d), rbind(rep(1, D * 2), c(1:D, 1:D), matrix(0, D-1, D*2)))
    
    if(control$warping == 'kY'){
      map <- function(z, A){
        if(is.null(nrow(z)))
          z <- matrix(z, nrow = 1)
        return(z)
      }
    }
    if(control$warping == 'kX'){
      map <- mapZX
      formals(map)$Amat <- Amat
      formals(map)$Aind <- Aind
    }
    if(control$warping == 'Psi'){
      map <- Psi_Z
      formals(map)$Amat <- Amat
      formals(map)$Aind <- Aind
    }
    ##
    # now build design for Rembo
    design <- (DoE * 2 - 1) # rescale to [-1,1]^D
    design <- ortProj(design, tA)
    
    if(!highDimGP) design <- map(design, A_hat) else design <- DoE * 2 - 1
    
    minpar <- DoE[which.min(fvalues),]
    
    # newmodel <- try(update(model, newDesign, newY, cov.reestim = kmcontrol$covreestim,
    #                        kmcontrol = list(control = list(trace = FALSE))))
    
    newmodel <- try(km(kmcontrol$formula, design = design, response = fvalues, covtype = kmcontrol$covtype, 
                       iso = kmcontrol$iso, nugget.estim = nuggetestim,
                       control = list(trace = FALSE)))
    
    # if(typeof(newmodel) == "character"){
    #   cat("Error in hyperparameter estimation - old hyperparameter values used instead for model at iteration", model@n - n.init, "\n")
    #   newmodel <- try(update(object = model, newX = newDesign, newy = newY, newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
    # }
    
    if (typeof(newmodel) == "character") {
      cat("Unable to udpate kriging model at iteration", model@n - n.init, "- optimization stopped \n")
      cat("lastmodel is the model at iteration", model@n - n.init -1, "\n")
      if(control$lightreturn){
        res <- list(par = minpar, value = min(fvalues))
      }else{
        res <- list(par = minpar,
                    value = min(fvalues),
                    low_dim_design = DoE,
                    y = fvalues,
                    A = A_hat,
                    model = model)
      }
      return(res)
    }
    model <- newmodel
    print(min(fvalues))
  }
  
  if(control$lightreturn){
    res <- list(par = minpar, value = min(fvalues))
  }else{
    res <- list(par = minpar,
                value = min(fvalues),
                low_dim_design = DoE,
                y = fvalues,
                A = A_hat,
                model = model,
                homMod = modelD)
  }
  
  return(res)
  
}



