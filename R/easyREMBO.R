##' REMBO for unconstrained problems
##' @title Random Embedding Bayesian Optimization
##' @param par vector whose length is used to define the initial DoE or just the low dimension
##' @param fn function to be minimized over
##' @param ... additional parameters of fn
##' @param lower,upper bounds for optimization
##' @param budget total number of calls to the objective function
##' @param kmcontrol an optional list of control parameters to be passed to the \code{\link[DiceKriging]{km}} model:
##' \code{iso}, \code{covtype}, \code{formula}. In addition, boolean \code{codereestim} is passed to \code{\link[DiceKriging]{update.km}}
##' @param control an optional list of control parameters. See "Details"
##' @param init optional list with elements \code{Amat} to provide a random matrix, \code{low_dim_design} for an initial design in the low-dimensional space and \code{fvalues} for the corresponding response.
##' When passing initial response values, care should be taken that the mapping with \code{Amat} of the design actually correspond to high-dimensional designs giving \code{fvalues}.
##' @details Options available from \code{control} are:
##' \itemize{
##' \item \code{Atype} see \code{\link[RRembo]{selectA}};
##' \item \code{reverse} if \code{TRUE}, use the new mapping from the zonotope,
##'  otherwise the original mapping with convex projection;
##' \item \code{bxsize} scalar controling the box size in the low-dimensional space;
##' \item \code{testU} with the regular mapping, set to \code{TRUE} to check that points are in U (to avoid non-injectivity);
##' \item \code{standard} for using settings of the original REMBO method;
##' \item \code{maxitOptA} if \code{Atype} is \code{optimized}, number of optimization iterations;
##' \item \code{lightreturn} only returns \code{par} and \code{value};
##' \item \code{warping} either \code{"standard"} for kY, \code{"kX"} or \code{"Psi"};
##' \item \code{designtype} one of "\code{LHS}", "\code{maximin}" and "\code{unif}",
##'  see \code{\link[RRembo]{designZ}} or \code{\link[RRembo]{designU}};
##' \item \code{tcheckP} minimal distance to an existing solution, see \code{\link[GPareto]{checkPredict}}
##' \item \code{roll} to alternate between optimization methods;
##' \item \code{inneroptim} optimization method for EI
##' \item \code{popsize, gen} population size and number of optimization generations of EI
##' }
##' @importFrom rgenoud genoud
##' @importFrom DiceOptim EI
##' @importFrom GPareto checkPredict
##' @import DiceKriging
##' @importFrom pso psoptim
##' @importFrom DEoptim DEoptim
##' @importFrom rgenoud genoud
##' @import OOR
##' @author Mickael Binois
##' @export
##' @references
##' O. Roustant, D. Ginsbourger & Y. Deville, DiceKriging, DiceOptim: Two R Packages for the Analysis of Computer Experiments by Kriging-Based Metamodeling and Optimization Journal of Statistical Software, 2012, 51, 1-55\cr \cr
##' Z. Wang, F. Hutter, M. Zoghi, D. Matheson, N. de Freitas (2016), Bayesian Optimization in a Billion Dimensions via Random Embeddings, JAIR. \cr \cr
##' M. Binois, D. Ginsbourger, O. Roustant (2015), A Warped Kernel Improving Robustness in Bayesian Optimization Via Random Embeddings, Learning and Intelligent Optimization, Springer \cr \cr
##' M. Binois, D. Ginsbourger, O. Roustant (2018), On the choice of the low-dimensional domain for global optimization via random embeddings, arXiv:1704.05318 \cr \cr
##' M. Binois (2015), Uncertainty quantification on Pareto fronts and high-dimensional strategies in Bayesian optimization, with applications in multi-objective automotive design, PhD thesis, Mines Saint-Etienne.
##' @examples
##' \dontrun{
##' set.seed(42)
##' library(rgl)
##' library(DiceKriging)
##'
##' lowd <- 2
##' highD <- 25
##'
##' ntest <- 50
##' maxEval <- 100
##'
##' res <- rep(0, ntest)
##' for(i in 1:ntest){
##'   #ii <- sample(1:highD, 2)
##'   ii <- c(1,2)
##'   branin_mod2 <- function(X){
##'     if(is.null(nrow(X))) X <- matrix(X, nrow = 1)
##'     X <- X[, c(ii[1], ii[2]), drop = FALSE]
##'     return(apply(X, 1, branin))
##'   }
##'   sol <- easyREMBO(par = rep(NA, lowd), branin_mod2, lower = rep(0, highD),
##'                    upper = rep(1, highD), budget = maxEval)
##'   res[i] <- sol$value
##'   cat(sol$value, " ", i, "\n")
##' }
##'
##' plot(res - 0.397887, type = "b")
##' boxplot(res - 0.397887)
##' }
easyREMBO <- function(par, fn, lower, upper, budget, ...,
                      kmcontrol = list(covtype = "matern5_2", iso = TRUE, covreestim = TRUE, formula =~1),
                      control = list(Atype = 'isotropic', reverse = TRUE, bxsize = NULL, testU = TRUE, standard = FALSE,
                                     maxitOptA = 100, lightreturn = FALSE, warping = 'Psi', designtype = 'unif',
                                     tcheckP = 1e-4, roll = F,
                                     inneroptim = "pso", popsize = 80, gen = 40),
                      init = NULL){
  # Initialisation
  d <- length(par)
  D <- length(lower)
  
  if(is.null(control$Atype)) control$Atype <- 'isotropic'
  # if(is.null(control$Ysetting)) control$Ysetting <- 'max'
  if(is.null(control$testU)) control$testU <- TRUE
  if(is.null(control$standard)) control$standard <- FALSE
  if(is.null(control$maxitOptA)) control$maxitOptA <- 100
  if(is.null(control$lightreturn)) control$lightreturn <- FALSE
  if(is.null(control$warping)) control$warping <- 'Psi'
  if(is.null(control$inneroptim)) control$inneroptim <- 'pso'
  if(is.null(control$popsize)) control$popsize <- 80
  if(is.null(control$gen)) control$gen <- 40
  if(is.null(control$designtype)) control$designtype <- 'unif'
  if(is.null(control$reverse)) control$reverse <- TRUE
  if(is.null(control$maxf)) control$maxf <- control$popsize * control$gen
  if(is.null(control$tcheckP)) control$tcheckP <- 1e-4
  if(is.null(control$roll)) control$roll <- FALSE
  
  if(is.null(kmcontrol$covtype)) kmcontrol$covtype <- 'matern5_2'
  if(is.null(kmcontrol$iso)) kmcontrol$iso <- TRUE
  if(is.null(kmcontrol$covreestim)) kmcontrol$covreestim <- TRUE
  if(is.null(kmcontrol$formula)) kmcontrol$formula <- ~1
  
  if(!is.null(init$fvalues) && (is.null(init$Amat) || is.null(init$low_dim_design)))
    stop("When providing initial fvalues, both Amat and low_dim_design must be provided.")
  
  # Selection of the random embedding matrix
  if(is.null(init$Amat)){
    A <- selectA(d, D, type = control$Atype, control = list(maxit = control$maxitOptA))
  }else{
    A <- init$Amat
  }
  
  if(d == D){
    A <- diag(D)
    control$bxsize <- 1
    bxsize <- 1
  }
    
  
  tA <- t(A)
  
  # Selecting low dimensional domain ([-bxsize, bxsize]^d)
  if(is.null(control$bxsize)){
    if(control$standard){
      bxsize <- sqrt(d)
    }else{
      bxsize <- sqrt(D)
    }
  }
  
  # Mapping functions given the choices
  if(!control$reverse){
    if(control$warping == 'kY'){
      map <- function(y, A){
        if(is.null(nrow(y)))
          y <- matrix(y, nrow = 1)
        return(y)
      }
    }
    if(control$warping == 'kX'){
      map <- randEmb
    }
    if(control$warping == 'Psi'){
      if(control$Atype == 'standard'){
        map <- Psi_Y_nonort
        formals(map)$pA <- ginv(t(A) %*% A) %*% t(A)
        formals(map)$invA <- ginv(A)
      }else{
        map <- Psi_Y
      }
      
    }
  }else{
    # precomputations
    Amat <- cbind(A, matrix(rep(c(1, 0), times = c(1, D-1)), D, D), matrix(rep(c(-1, 0), times = c(1, D-1)), D, D))
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
  }
  
  # Number of points of the DoE
  if(is.null(init$n) && is.null(init$low_dim_design)){
    n.init <- max(4 * d, round(budget/3))
    ## This is a DiceKriging issue, comment lines 48-50 of "kmStruct.R"
    # if(control$warping == "kX"){
    #   if(n.init < D + 1 && D + 1 < budget){
    #     n.init <- D + 1
    #   }else{
    #     print("Budget too small for warping kX, need at least D+2")
    #   }
    # }
  }else{
    if(!is.null(init$n)) n.init <- init$n
    if(!is.null(init$low_dim_design)) n.init <- 0 # For update error consistency
  }
  
  if(control$reverse){
    # Create design
    if(is.null(init$low_dim_design)){
      DoE <- designZ(n.init, tA, bxsize, type = control$designtype)
    }else{
      # Or use the one provided
      indtest <- testZ(init$low_dim_design, tA) # check that provided designs are actually in the zonotope when using the new mapping
      if(!all(indtest)) warning("Not all initial low dimensional designs belong to Z.")
      DoE <- init$low_dim_design
    }
    
    if(is.null(init$fvalues)){
      fvalues <- fn(((mapZX(DoE, A, Amat = Amat, Aind = Aind) + 1 )/2) %*% diag(upper - lower) + matrix(lower, nrow = nrow(DoE), ncol = length(lower), byrow = T), ...)
    }else{
      fvalues <- init$fvalues
      if(length(fvalues) != nrow(DoE)) stop("Number of rows of design and length of response of provided designs do not match.")
    }
  }else{
    # Create design
    if(is.null(init$low_dim_design)){
      DoE <- designU(n.init, A, bxsize, type = control$designtype, standard = control$standard)
    }else{
      DoE <- init$low_dim_design
    }
    if(is.null(init$fvalues)){
      fvalues <- fn(((randEmb(DoE, A) + 1 )/2) %*% diag(upper - lower) + matrix(lower, nrow = nrow(DoE), ncol = length(lower), byrow = T), ...)
    }else{
      fvalues <- init$fvalues
      if(length(fvalues) != nrow(DoE)) stop("Number of rows of design and length of response of provided designs do not match.")
    }
  }
  cat("Initial best value: ", min(fvalues), "\n")
  design <- map(DoE, A)
  
  model <- try(km(kmcontrol$formula, design = design, response = fvalues, covtype = kmcontrol$covtype, iso = kmcontrol$iso,
              control = list(trace = FALSE)))
  if(class(model) == "try-error"){
    model <- km(kmcontrol$formula, design = design, response = fvalues, covtype = kmcontrol$covtype, iso = kmcontrol$iso,
                control = list(trace = FALSE), nugget.estim = T)
  }
  
  # Expected Improvement function in the REMBO case
  # Not evaluated and penalized if not in U/Z or
  # mval: max penalty value (negative value)
  EI_Rembo <- function(x, model, mval = -10){
    if(is.null(nrow(x)))
      x <- matrix(x, nrow = 1)
    
    inDomain <- rep(TRUE, nrow(x))
    
    if(control$reverse){
      inDomain <- testZ(x, tA)
    }else{
      if(control$testU){
        inDomain <- testU(x, A)
      }
    }
    res <- rep(NA, nrow(x))
    if(any(inDomain)){
      xtmp <- map(x[inDomain,], A)
      # identify too close points
      tmp <- checkPredict(xtmp, list(model), control$tcheckP, distance = 'euclidean') 
      
      res[inDomain[tmp]] <- 0
      if(any(!tmp)) res[inDomain[!tmp]] <- EI(x = xtmp[!tmp,], model = model)
    }
    # else{
    #   return(mval * distance(x, 0)) # penalty
    # }
    if(any(!inDomain)) res[!inDomain] <- mval * apply(x[!inDomain,,drop = FALSE], 1, distance, x2 = 0)
    return(res)
  }
  
  if(control$reverse){
    boundsEIopt <- rowSums(abs(tA)) ## box surrounding the zonotope Z
  }else{
    boundsEIopt <- rep(bxsize, d)
  }
  
  # best value so far, slightly perturbed, to increase exploitation in EI optimization
  spartan <- pmin(boundsEIopt, pmax(-boundsEIopt, DoE[which.min(fvalues),] + rnorm(d, sd = 0.05)))
  
  while(model@n < budget){
    
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
                     control = list(fnscale = -1, maxit = control$gen, s = control$popsize, vectorize = T), model = model)
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
    if(control$reverse){
      newY <- fn(((mapZX(opt$par, A, Amat = Amat, Aind = Aind) + 1)/2) %*% diag(upper - lower) + lower, ...)
    }else{
      newY <- fn(((randEmb(opt$par, A) +  1)/2) %*% diag(upper - lower) + lower, ...)
    }
    if(newY < min(fvalues)) cat(model@n + 1, ": ", newY, "\n")
    
    newX <- opt$par
    newDesign <- map(opt$par, A)
    fvalues <- c(fvalues, newY)
    DoE <- rbind(DoE, newX)
    design <- rbind(design, newDesign)
    
    ind <- which.min(fvalues)
    
    spartan <- pmin(boundsEIopt, pmax(-boundsEIopt, DoE[which.min(fvalues),] + rnorm(d, sd = 0.05)))
    
    if(control$reverse){
      minpar <- ((mapZX(DoE[ind,], A, Amat = Amat, Aind = Aind)+1)/2) %*% diag(upper - lower) + matrix(lower, nrow = 1, ncol = length(lower), byrow = T)
    }else{
      minpar <- ((randEmb(DoE[ind,], A)+1)/2) %*% diag(upper - lower) + matrix(lower, nrow = 1, ncol = length(lower), byrow = T)
    }
    
    newmodel <- try(update(model, newDesign, newY, cov.reestim = kmcontrol$covreestim,
                           kmcontrol = list(control = list(trace = FALSE))))
    
    if(typeof(newmodel) == "character"){
      cat("Error in hyperparameter estimation - old hyperparameter values used instead for model at iteration", model@n - n.init, "\n")
      newmodel <- try(update(object = model, newX = newDesign, newy = newY, newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
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
                    A = A,
                    model = model)
      }
      return(res)
    }
    model <- newmodel
    # print(min(fvalues))
  }
  
  if(control$lightreturn){
    res <- list(par = minpar, value = min(fvalues))
  }else{
    res <- list(par = minpar,
                value = min(fvalues),
                low_dim_design = DoE,
                y = fvalues,
                A = A,
                model = model)
  }
  
  return(res)
  
}

