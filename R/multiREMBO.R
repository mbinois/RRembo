##' REMBO for multi-objective unconstrained problems
##' @title Random Embedding Bayesian Optimization
##' @param par vector whose length is used to define the initial DoE or just the low dimension
##' @param fn function to be minimized over
##' @param ... additional parameters of fn
##' @param critcontrol optional controls for hypervolume computations, see \code{\link[GPareto]{crit_EHI}}
##' @param lower,upper bounds for optimization
##' @param budget total number of calls to the objective function
##' @param control an optional list of control parameters. See "Details"
##' @details Options available from control are:
##' \itemize{
##' \item \code{Atype} see \code{\link[RRembo]{selectA}};
##' \item \code{reverse} if \code{TRUE}, use the new mapping from the zonotope,
##'  otherwise the original mapping with convex projection;
##' \item \code{Amat} matrix defining the random embedding;
##' \item \code{bxsize} scalar controling the box size in the low-dimensional space;
##' \item \code{testU} with the regular mapping, set to \code{TRUE} to check that points are in U (to avoid non-injectivity);
##' \item \code{standard} for using settings of the original REMBO method;
##' \item \code{maxitOptA} if \code{Atype} is \code{optimized}, number of optimization iterations;
##' \item \code{lightreturn} only returns \code{par} and \code{value};
##' \item \code{warping} either \code{"standard"} for kY, \code{"kX"} or \code{"Psi"};
##' \item \code{covtype}, \code{covreestim}, \code{iso} see \code{\link[DiceKriging]{km}};
##' \item \code{designtype} one of "LHS", "maximin" and 'unif',
##'  see \code{\link[RRembo]{designZ}} or \code{\link[RRembo]{designU}};
##' \item \code{tcheckP} minimal distance to an existing solution, see \code{\link[GPareto]{checkPredict}}
##' \item \code{roll} to alternate between optimization methods;
##' \item \code{initdesigns} initial design matrix
##' \item \code{inneroptim} optimization method for EI
##' \item \code{popsize, gen} population size and number of optimization generations of EI
##' }
##' @importFrom GPareto plotParetoEmp
##' @importFrom eaf is.nondominated
##' @importFrom GPareto crit_EHI
##' @export
##' @examples
##' \dontrun{
##' set.seed(42)
##' library(GPareto)
##' library(eaf)
##'
##' lowd <- 2
##' highD <- 25
##'
##' maxEval <- 100
##'
##' n.grid <- 101
##' test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
##' response.grid <- t(apply(test.grid, 1, P1))
##' PFref <- response.grid[is.nondominated(response.grid),]
##' 
##' ii <- c(1,2)
##' P1_mod <- function(X){
##'   if(is.null(nrow(X))) X <- matrix(X, nrow = 1)
##'   return(P1(X[,c(ii[1], ii[2]), drop = FALSE]))
##' }
##' 
##' plot(response.grid, pch = '.', xlab = "f1", ylab = "f2")                   
##' plotParetoEmp(PFref, col = "red")
##' 
##' sol <- multiREMBO(par = rep(NA, lowd), P1_mod, lower = rep(0, highD),
##'                    upper = rep(1, highD), budget = maxEval)
##' points(sol$values, col = "blue", pch = 20)
##'
##' }
multiREMBO <- function(par, fn, lower, upper, budget, ..., critcontrol = NULL, 
                       control = list(Atype = 'isotropic', reverse = TRUE, Amat = NULL, bxsize = NULL, testU = TRUE, standard = FALSE,
                                      maxitOptA = 100, lightreturn = FALSE, warping = 'Psi', covtype = "matern5_2", covreestim = TRUE, designtype = 'unif',
                                      tcheckP = 1e-5, roll = F, initdesigns = NULL, iso = T,
                                      inneroptim = "pso", popsize = 80, gen = 40)){
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
  if(is.null(control$covtype)) control$covtype <- 'matern5_2'
  if(is.null(control$iso)) control$iso <- TRUE
  if(is.null(control$covreestim)) control$covreestim <- TRUE
  if(is.null(control$inneroptim)) control$inneroptim <- 'pso'
  if(is.null(control$popsize)) control$popsize <- 80
  if(is.null(control$gen)) control$gen <- 40
  if(is.null(control$designtype)) control$designtype <- 'unif'
  if(is.null(control$reverse)) control$reverse <- TRUE
  if(is.null(control$maxf)) control$maxf <- control$popsize * control$gen
  if(is.null(control$tcheckP)) control$tcheckP <- 1e-4
  if(is.null(control$roll)) control$roll <- FALSE
  if(is.null(critcontrol$refPoint)) noRef <- TRUE else noRef <- FALSE
  if(is.null(critcontrol$extendper)) critcontrol$extendper <- 0.2
  
  # Selection of the random embedding matrix
  if(is.null(control$Amat)){
    A <- selectA(d, D, type = control$Atype, control = list(maxit = control$maxitOptA))
  }else{
    A <- control$Amat
  }
  
  if(d == D)
    A <- diag(D)
  
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
    if(control$warping == 'kpA'){
      map <- function(y, A){
        return(ortProj(randEmb(y, A), tA))
      }
    }
    if(control$warping == 'Psi'){
      if(control$Atype == 'standard'){
        map <- Psi_Y_nonort
      }else{
        map <- Psi_Y
      }
      
    }
  }else{
    if(control$warping == 'kY' || control$warping == 'pA'){
      map <- function(z, A){
        if(is.null(nrow(z)))
          z <- matrix(z, nrow = 1)
        return(z)
      }
    }
    if(control$warping == 'kX'){
      map <- mapZX
    }
    if(control$warping == 'Psi'){
      map <- Psi_Z
    }
  }
  
  # Number of points of the DoE
  n.init <- max(4 * d, round(budget/3))
  
  if(control$reverse){
    DoE <- designZ(n.init, tA, bxsize, type = control$designtype)
    
    if(!is.null(control$initdesigns)){
      indtest <- rep(TRUE, nrow(control$initdesigns))
      if(control$reverse){
        indtest <- testZ(control$initdesigns, tA)
      }
      # else: not yet implemented
      if(any(indtest))
        DoE[which(indtest),] <- control$initdesigns
    }
    
    fvalues <- fn(((mapZX(DoE, A) + 1 )/2) %*% diag(upper - lower) + matrix(lower, nrow = nrow(DoE), ncol = length(lower), byrow = T), ...)
  }else{
    DoE <- designU(n.init, A, bxsize, type = control$designtype, standard = control$standard)
    fvalues <- fn(((randEmb(DoE, A) + 1 )/2) %*% diag(upper - lower) + matrix(lower, nrow = nrow(DoE), ncol = length(lower), byrow = T), ...)
  }
  
  design <- map(DoE, A)
  
  nobj <- ncol(fvalues)
  
  observations <- fvalues
  
  model <- list()
  for(i in 1:nobj){
    model[[i]] <- km(~1, design = design, response = fvalues[,i], covtype = control$covtype, iso = control$iso,
                     control = list(trace = FALSE))
  }
  
  
  # Expected Improvement function in the REMBO case
  # Not evaluated and penalized if not in U/Z or
  # mval: max penalty value (negative value)
  EHI_Rembo <- function(x, model, mval = -10, PF = NULL, critcontrol = NULL){
    if(is.null(nrow(x)))
      x <- matrix(x, nrow = 1)
    
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
      tmp <- checkPredict(xtmp, model, control$tcheckP, distance = 'euclidean') 
      
      res[inDomain[tmp]] <- 0
      if(any(!tmp)) res[inDomain[!tmp]] <- crit_EHI(x = xtmp[!tmp,], model = model, paretoFront = PF, critcontrol = critcontrol)
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
  ind <- is.nondominated(fvalues)
  PF <- fvalues[ind,]
  spartan <- DoE[ind,] + matrix(rnorm(d*nrow(PF), sd = 0.05), nrow(PF))
  
  while(model[[1]]@n < budget){
    
    if(noRef){
      PF_range <- apply(PF, 2, range)
      critcontrol$refPoint <- matrix(PF_range[2,] + pmax(1, (PF_range[2,] - PF_range[1,]) * critcontrol$extendper), 1, nobj)
    }
    
    # Change optimization method every time
    if(control$roll){
      if(model[[1]]@n %% 4 == 0) control$inneroptim <- 'genoud'
      if(model[[1]]@n %% 4 == 1) control$inneroptim <- 'DEoptim'
      if(model[[1]]@n %% 4 == 2) control$inneroptim <- 'genoud'
      if(model[[1]]@n %% 4 == 3) control$inneroptim <- 'DEoptim'
    }
    
    if(control$inneroptim == 'DEoptim'){
      EHI_2_DE <- function(x, model, mval = -10){
        return(-EHI_Rembo(x, model, mval))
      }
      parinit <- 2*matrix(runif(d*control$popsize), control$popsize) - 1
      parinit[1:nrow(spartan),] <- spartan
      opt <- DEoptim(fn = EHI_2_DE, lower = -boundsEIopt, upper = boundsEIopt, model = model,
                     control = list(initialpop = parinit, NP = control$popsize, itermax = control$gen, trace = F),
                     critcontrol = critcontrol, PF = PF)
      opt$value <- -opt$optim$bestval
      opt$par <- opt$optim$bestmem
    }
    
    if(control$inneroptim == 'pso'){
      parinit <- spartan[sample(1:nrow(spartan), 1),]
      opt <- psoptim(parinit, fn = EHI_Rembo, lower = -boundsEIopt, upper = boundsEIopt,
                     control = list(fnscale = -1, maxit = control$gen, s = control$popsize, vectorize = T), model = model,
                     critcontrol = critcontrol, PF = PF)
      opt$value <- -opt$value
    }
    
    if(control$inneroptim == "genoud"){
      parinit <- 2*matrix(runif(d*control$popsize), control$popsize) - 1
      parinit[1:nrow(spartan),] <- spartan
      opt <- genoud(EHI_Rembo, nvars = d, max = TRUE, starting.values = parinit,
                    Domains = cbind(-boundsEIopt, boundsEIopt), print.level = 0,
                    model = model,
                    critcontrol = critcontrol, PF = PF,
                    boundary.enforcement = 2, pop.size = control$popsize,
                    max.generations = control$gen)
    }
    
    if(control$inneroptim == 'StoSOO'){
      opt <- StoSOO(par = rep(NA, d), fn = EHI_Rembo, lower = -boundsEIopt, upper = boundsEIopt,
                    critcontrol = critcontrol, PF = PF,
                    settings = list(max = TRUE, type = 'det'), nb_iter = control$maxf, model = model)
    }
    
    if(opt$value <= 0){ # no improvement found or not even in domain
      cat("No improvement or not even in Z \n")
    }
    if(control$reverse){
      newY <- fn(((mapZX(opt$par, A) + 1)/2) %*% diag(upper - lower) + lower, ...)
    }else{
      newY <- fn(((randEmb(opt$par, A) +  1)/2) %*% diag(upper - lower) + lower, ...)
    }
    
    newX <- opt$par
    newDesign <- map(opt$par, A)
    fvalues <- rbind(fvalues, newY)
    DoE <- rbind(DoE, newX)
    design <- rbind(design, newDesign)
    
    ind <- is.nondominated(fvalues)
    PF <- fvalues[ind,]
    spartan <- DoE[ind,] + matrix(rnorm(d*nrow(PF), sd = 0.05), nrow(PF))
    
    if(control$reverse){
      minpar <- ((mapZX(DoE[ind,], A)+1)/2) %*% diag(upper - lower) + matrix(lower, nrow = nrow(PF), ncol = length(lower), byrow = T)
    }else{
      minpar <- ((randEmb(DoE[ind,], A)+1)/2) %*% diag(upper - lower) + matrix(lower, nrow = nrow(PF), ncol = length(lower), byrow = T)
    }
    
    ## Update models
    newmodel <- list()
    for(i in 1:nobj){
      newmodel[[i]] <- try(update(model[[i]], newDesign, newY[i], cov.reestim = control$covreestim,
                             kmcontrol = list(control = list(trace = FALSE))))
      
      if(typeof(newmodel[[i]]) == "character"){
        cat("Error in hyperparameter estimation - old hyperparameter values used instead for model at iteration", model[[i]]@n - n.init, "\n")
        newmodel[[i]] <- try(update(object = model[[i]], newX = newDesign, newy = newY[i], newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
      }
      
      if (typeof(newmodel) == "character") {
        cat("Unable to udpate kriging model at iteration", model[[i]]@n - n.init, "- optimization stopped \n")
        cat("lastmodel is the model at iteration", model[[i]]@n - n.init -1, "\n")
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
      
    }
    model <- newmodel
    plotParetoEmp(PF, col = "blue")
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

