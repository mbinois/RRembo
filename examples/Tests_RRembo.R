library("doParallel")
library("foreach")
library("RRembo")
library("DiceKriging")

for(case in 1:5){
  print(case)
  nCores <- detectCores() - 1
  cl <-  makeCluster(nCores)
  registerDoParallel(cl)
  
  set.seed(42)
  save_intermediate <- TRUE
  
  if(case == 1){
    d <- 6
    D <- 50
    nrep <- 25
    budget <- 250
    lower <- rep(0, D)
    upper <- rep(1, D)
    ftest <- hartman6_mod_log
    popsize <- 80
    covtype <- "matern5_2"
    roll <- T
    
    # To ensure fairness among runs
    mat_effective <- matrix(0, nrep, d) # matrix of effective variables (for D)
    for(i in 1:nrep){
      mat_effective[i,] <- sample(1:D, d)
    }
    fstar <- -1.200678
  }
  
  if(case == 2){
    d <- 2
    D <- 25
    nrep <- 25
    budget <- 100
    lower <- rep(0, D)
    upper <- rep(1, D)
    popsize <- 80
    covtype <- "matern3_2"
    roll <- T
    
    ftest <- branin_mod
    
    # To ensure fairness among runs
    mat_effective <- matrix(0, nrep, d) # matrix of effective variables (for D)
    for(i in 1:nrep){
      mat_effective[i,] <- sample(1:D, d)
    }
    
    fstar <- 0.3978874
    
  }
  
  if(case == 3){
    d <- 6
    D <- 17
    budget <- 250
    nrep <- 25
    lower <- rep(0, D)
    upper <- rep(1, D)
    popsize <- 80
    covtype <- "matern5_2"
    roll <- T
    
    cola_mod <- function(x, ii = NULL){
      if(is.null(dim(x))) x <- matrix(x, 1)
      return(apply(x, 1, cola))
    } 
    
    ftest <- cola_mod
    
    # To ensure fairness among runs
    mat_effective <-  NULL # matrix of effective variables (for D)
    
    fstar <- 11.7464
    
  }
  
  if(case == 4){
    d <- 10
    D <- 80
    budget <- 250
    nrep <- 25
    lower <- rep(0, D)
    upper <- rep(1, D)
    popsize <- 80
    covtype <- "matern5_2"
    roll <- T
    ftest <- levy
    
    # To ensure fairness among runs
    mat_effective <- matrix(0, nrep, d) # matrix of effective variables (for D)
    for(i in 1:nrep){
      mat_effective[i,] <- sample(1:D, d)
    }
    
    fstar <- 0
    
  }
  
  if(case == 5){
    d <- 2
    D <- 80
    budget <- 100
    nrep <- 25
    lower <- rep(0, D)
    upper <- rep(1, D)
    popsize <- 80
    covtype <- "matern5_2"
    roll <- T
    ftest <- giunta
    
    # To ensure fairness among runs
    mat_effective <- matrix(0, nrep, d) # matrix of effective variables (for D)
    for(i in 1:nrep){
      mat_effective[i,] <- sample(1:D, d)
    }
    
    fstar <- 0
    
  }
  
  cat('RO \n')
  # Random optimization in the high-dimensional space
  tRO <- Sys.time()
  res_RO <- foreach(i=1:nrep, .packages = c("RRembo", "DiceKriging"), .combine = 'rbind') %dopar% {
    set.seed(i)
    res <- ftest(matrix(runif(D*budget), budget), ii = mat_effective[i,])
  }
  tRO <- difftime(Sys.time(), tRO, units = "sec")
  
  ###
  
  cat('REMBO standard \n')
  # Standard REMBO-like method
  tsY <- Sys.time()
  res_standard_kY <- foreach(i=1:nrep, .packages = c("RRembo", "DiceKriging"), .combine = 'rbind') %dopar% {
    set.seed(i)
    res <- easyREMBO(par = runif(d), fn = ftest, budget = budget, lower = lower, upper = upper,
                     ii = mat_effective[i,],
                     control = list(Atype = 'standard', reverse = FALSE, warping = 'kY', testU = FALSE, standard = TRUE, popsize = popsize, gen = 40,
                                    inneroptim = "StoSOO", covtype = covtype, roll = roll))
    res$y
  }
  tsY <- difftime(Sys.time(), tsY, units = "sec")
  if(save_intermediate) save.image(paste0("Wksp_RRembo_case_", case, ".RData"))
  
  ###
  
  cat('REMBO standard k_X \n')
  # Standard REMBO-like method
  tsX <- Sys.time()
  res_standard_kX <- foreach(i=1:nrep, .packages = c("RRembo", "DiceKriging"), .combine = 'rbind') %dopar% {
    set.seed(i)
    res <- easyREMBO(par = runif(d), fn = ftest, budget = budget, lower = lower, upper = upper,
                     ii = mat_effective[i,],
                     control = list(Atype = 'standard', reverse = FALSE, warping = 'kX', testU = FALSE, standard = TRUE, popsize = popsize, gen = 40,
                                    inneroptim = "StoSOO", covtype = covtype, roll = roll))
    res$y
  }
  tsX <- difftime(Sys.time(), tsX, units = "sec")
  if(save_intermediate) save.image(paste0("Wksp_RRembo_case_", case, ".RData"))
  
  ###
  
  cat('REMBO standard k_Psi \n')
  # Standard REMBO-like method
  tsP <- Sys.time()
  res_standard_kPsi <- foreach(i=1:nrep, .packages = c("RRembo", "DiceKriging"), .combine = 'rbind') %dopar% {
    set.seed(i)
    res <- easyREMBO(par = runif(d), fn = ftest, budget = budget, lower = lower, upper = upper,
                     ii = mat_effective[i,],
                     control = list(Atype = 'standard', reverse = FALSE, warping = 'Psi', testU = FALSE, standard = TRUE, popsize = popsize, gen = 40,
                                    inneroptim = "StoSOO", covtype = covtype, roll = roll))
    res$y
  }
  tsP <- difftime(Sys.time(), tsP, units = "sec")
  if(save_intermediate) save.image(paste0("Wksp_RRembo_case_", case, ".RData"))
  
  ###
  
  cat('REMBO reverse + A Gaussian + Psi \n')
  trP <- Sys.time()
  res_reverse_kPsi <- foreach(i=1:nrep, .packages = c("RRembo", "DiceKriging"), .combine = 'rbind') %dopar% {
    set.seed(i)
    res <- easyREMBO(par = runif(d), fn = ftest, budget = budget, lower = lower, upper = upper,
                     ii = mat_effective[i,],
                     control = list(Atype = 'Gaussian', reverse = TRUE, warping = 'Psi', testU = TRUE, standard = FALSE,  popsize = popsize, gen = 40,
                                    inneroptim = "StoSOO", covtype = covtype, roll = roll))
    res$y
  }
  trP <- difftime(Sys.time(), trP, units = "sec")
  if(save_intermediate) save.image(paste0("Wksp_RRembo_case_", case, ".RData"))
  
  ###
  
  cat('REMBO reverse + A Gaussian + kY \n')
  trY <- Sys.time()
  res_reverse_kY <- foreach(i=1:nrep, .packages = c("RRembo", "DiceKriging"), .combine = 'rbind') %dopar% {
    set.seed(i)
    res <- easyREMBO(par = runif(d), fn = ftest, budget = budget, lower = lower, upper = upper,
                     ii = mat_effective[i,],
                     control = list(Atype = 'Gaussian', reverse = TRUE, warping = 'kY', testU = TRUE, standard = FALSE,  popsize = popsize, gen = 40,
                                    inneroptim = "StoSOO", covtype = covtype, roll = roll))
    res$y
  }
  trY <- difftime(Sys.time(), trY, units = "sec")
  if(save_intermediate) save.image(paste0("Wksp_RRembo_case_", case, ".RData"))
  
  ###
  
  cat('REMBO reverse + Gaussian + kX \n')
  trX <- Sys.time()
  res_reverse_kX <- foreach(i=1:nrep, .packages = c("RRembo", "DiceKriging"), .combine = 'rbind') %dopar% {
    set.seed(i)
    res <- easyREMBO(par = runif(d), fn = ftest, budget = budget, lower = lower, upper = upper,
                     ii = mat_effective[i,],
                     control = list(Atype = 'Gaussian', reverse = TRUE, warping = 'kX', testU = TRUE, standard = FALSE,  popsize = popsize, gen = 40,
                                    inneroptim = "StoSOO", covtype = covtype, roll = roll))
    res$y
  }
  trX <- difftime(Sys.time(), trX, units = "sec")
  if(save_intermediate) save.image(paste0("Wksp_RRembo_case_", case, ".RData"))
  
  ### 
  
  stopCluster(cl)
  
  
  # Post-treatment for progression
  rRO <- res_RO
  rskY <- res_standard_kY; rskX <- res_standard_kX; rskP <- res_standard_kPsi
  rrkY <- res_reverse_kY; rrkX <- res_reverse_kX; rrkP <- res_reverse_kPsi
  
  for(i in 1:nrep){
    for(j in 2:budget){
      if(rRO[i, j] > rRO[i, j-1]) rRO[i, j] <- rRO[i, j-1]
      if(rskY[i, j] > rskY[i, j-1]) rskY[i, j] <- rskY[i, j-1]
      if(rskX[i, j] > rskX[i, j-1]) rskX[i, j] <- rskX[i, j-1]
      if(rskP[i, j] > rskP[i, j-1]) rskP[i, j] <- rskP[i, j-1]
      if(rrkY[i, j] > rrkY[i, j-1]) rrkY[i, j] <- rrkY[i, j-1]
      if(rrkX[i, j] > rrkX[i, j-1]) rrkX[i, j] <- rrkX[i, j-1]
      if(rrkP[i, j] > rrkP[i, j-1]) rrkP[i, j] <- rrkP[i, j-1]
    }
  }
  
  # Plot final results
  boxplot(rRO[,budget] - fstar, 
          rskY[,budget] - fstar, rskX[,budget] - fstar, rskP[,budget] - fstar,
          rrkY[,budget] - fstar, rrkX[,budget] - fstar, rrkP[,budget] - fstar)
  
  
  plot(apply(rRO, 2, mean), type = 'l', lwd = 2)
  lines(apply(rskY, 2, mean), col = 2, lwd = 2)
  lines(apply(rskX, 2, mean), col = 3, lwd = 2)
  lines(apply(rskP, 2, mean), col = 4, lwd = 2)
  lines(apply(rrkY, 2, mean), col = 5, lwd = 2)
  lines(apply(rrkX, 2, mean), col = 6, lwd = 2)
  lines(apply(rrkP, 2, mean), col = 7, lwd = 2)
  
  legend("topright", legend = c("RO", "std kY", "std kX", "std kP",
                                "rev kY", "rev kX", "rev kP"), col = 1:7, lty = 1)
  
}

