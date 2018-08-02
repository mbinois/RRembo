library("doParallel")
library("foreach")

nCores <- detectCores() - 1
cl <-  makeCluster(nCores)
registerDoParallel(cl)
library(RemboIV)
library(DiceKriging)

set.seed(42)

case <- 2
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
  popsize <- 40
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


cat('RO \n')
# Random optimization in the high-dimensional space
res_RO <- foreach(i=1:nrep, .packages = c("RemboIV", "DiceKriging"), .combine = 'c') %dopar% {
  set.seed(i)
  res <- min(ftest(matrix(runif(D*budget), budget), ii = mat_effective[i,]))
}

cat('REMBO standard \n')
# Standard REMBO-like method
res_standard_kY <- foreach(i=1:nrep, .packages = c("RemboIV", "DiceKriging"), .combine = 'c') %dopar% {
  set.seed(i)
  res <- easyREMBO(par = runif(d), fn = ftest, budget = budget, lower = lower, upper = upper,
               ii = mat_effective[i,],
               control = list(Atype = 'standard', reverse = FALSE, warping = 'kY', testU = FALSE, standard = TRUE, popsize = popsize, gen = 40,
                              inneroptim = "StoSOO", covtype = covtype, roll = roll))
  res$value
}

cat('REMBO standard k_X \n')
# Standard REMBO-like method
res_standard_kX <- foreach(i=1:nrep, .packages = c("RemboIV", "DiceKriging"), .combine = 'c') %dopar% {
  set.seed(i)
  res <- easyREMBO(par = runif(d), fn = ftest, budget = budget, lower = lower, upper = upper,
                ii = mat_effective[i,],
                control = list(Atype = 'standard', reverse = FALSE, warping = 'kX', testU = FALSE, standard = TRUE, popsize = popsize, gen = 40,
                               inneroptim = "StoSOO", covtype = covtype, roll = roll))
  res$value
}

cat('REMBO standard k_Psi \n')
# Standard REMBO-like method
res_standard_kPsi <- foreach(i=1:nrep, .packages = c("RemboIV", "DiceKriging"), .combine = 'c') %dopar% {
  set.seed(i)
  res <- easyREMBO(par = runif(d), fn = ftest, budget = budget, lower = lower, upper = upper,
                ii = mat_effective[i,],
                control = list(Atype = 'standard', reverse = FALSE, warping = 'Psi', testU = FALSE, standard = TRUE, popsize = popsize, gen = 40,
                               inneroptim = "StoSOO", covtype = covtype, roll = roll))
  res$value
}

cat('REMBO reverse + A Gaussian + Psi \n')
res_reverse_G_kPsi <- foreach(i=1:nrep, .packages = c("RemboIV", "DiceKriging"), .combine = 'c') %dopar% {
  set.seed(i)
  res <- easyREMBO(par = runif(d), fn = ftest, budget = budget, lower = lower, upper = upper,
                ii = mat_effective[i,],
                control = list(Atype = 'Gaussian', reverse = TRUE, warping = 'Psi', testU = TRUE, standard = FALSE,  popsize = popsize, gen = 40,
                               inneroptim = "StoSOO", covtype = covtype, roll = roll))
  res$value
}

cat('REMBO reverse + A Gaussian + kY \n')
res_reverse_G_kY <- foreach(i=1:nrep, .packages = c("RemboIV", "DiceKriging"), .combine = 'c') %dopar% {
  set.seed(i)
  res <- easyREMBO(par = runif(d), fn = ftest, budget = budget, lower = lower, upper = upper,
                ii = mat_effective[i,],
                control = list(Atype = 'Gaussian', reverse = TRUE, warping = 'kY', testU = TRUE, standard = FALSE,  popsize = popsize, gen = 40,
                               inneroptim = "StoSOO", covtype = covtype, roll = roll))
  res$value
}



cat('REMBO reverse + Gaussian + kX \n')
res_reverse_G_kX <- foreach(i=1:nrep, .packages = c("RemboIV", "DiceKriging"), .combine = 'c') %dopar% {
  set.seed(i)
  res <- easyREMBO(par = runif(d), fn = ftest, budget = budget, lower = lower, upper = upper,
                ii = mat_effective[i,],
                control = list(Atype = 'Gaussian', reverse = TRUE, warping = 'kX', testU = TRUE, standard = FALSE,  popsize = popsize, gen = 40,
                               inneroptim = "StoSOO", covtype = covtype, roll = roll))
  res$value
}

stopCluster(cl)

boxplot(res_RO - fstar, 
        res_standard_kY - fstar, res_standard_kX - fstar, res_standard_kPsi - fstar,
        res_reverse_G_kY - fstar, res_reverse_G_kX - fstar, res_reverse_G_kPsi - fstar)


