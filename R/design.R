##' Select design in U
##' @param p number of points of the design
##' @param A random embedding matrix
##' @param bxsize scaler bounds on the domain, i.e., bxsize x [-1,1]^D
##' @param type design type, one of "LHS", "maximin" and 'unif'
##' @param standard standard REMBO approach
##' @author Mickael Binois
##' @references 
##' M. Binois, D. Ginsbourger, O. Roustant (2018), On the choice of the low-dimensional domain for global optimization via random embeddings, arXiv:1704.05318 \cr \cr
##' M. Binois (2015), Uncertainty quantification on Pareto fronts and high-dimensional strategies in Bayesian optimization, with applications in multi-objective automotive design, PhD thesis, Mines Saint-Etienne.
##' @importFrom lhs maximinLHS randomLHS augmentLHS optAugmentLHS
##' @export
##' @examples
##' ## Example of designs in U
##' set.seed(42)
##' d <- 2; D <- 5
##'
##' A <- selectA(d, D, type = 'optimized')
##' size <- 5 # box size of Y
##' ntest <- 10000
##' Y <- size * (2 * matrix(runif(ntest * d), ntest, d) - 1)
##'
##' inU <- testU(Y, A)
##' colors <- rep('black', ntest)
##' colors[inU] <- 'green'
##' plot(Y, col = colors, pch = 20, cex = 0.5)
##'
##' p <- 20
##' designs1 <- designU(p, A, size)
##' designs2 <- designU(p, A, size, type = 'LHS')
##'
##' points(designs1, col = 'red', pch = 20)
##' points(designs2, col = 'blue', pch = 20)
##' legend("topright", legend = c("LHS", "Maximin LHS"), col = c("blue", "red"), pch = c(20,20))
##'
designU <- function(p, A, bxsize, type = "maximin", standard = FALSE){
  d <- ncol(A)
  n <- 0
  designs <- NULL
  while(n < p){
    #initialisation
    if(is.null(designs)){
      if(type == 'maximin')
        designs <- maximinLHS(p, d)
      if(type == 'LHS')
        designs <- randomLHS(p, d)
      if(type == 'unif')
        designs <- matrix(runif(p*d), p)
    }else{
      if(type == 'maximin')
        designs <- augmentLHS(designs, p - n)
      if(type == 'LHS')
        designs <- optAugmentLHS(designs, p - n)
      if(type == 'unif')
        designs <- rbind(designs, matrix(runif((p - n)*d), p - n))
    }
    if(standard){
      ind <- !duplicated(randEmb(2 * bxsize * designs - bxsize, A))
    }else{
      ind <- testU(2 * bxsize * designs - bxsize, A)
    }
    
    n <- sum(ind)
  }
  return(2 * bxsize * designs[ind,] - bxsize)
  
}


##' Select design in Z
##' @param p number of points of the design
##' @param pA orthogonal projection onto Ran(A) matrix
##' @param bxsize bounds of the domain
##' @param type design type, one of "LHS", "maximin"
##' @author Mickael Binois
##' @references 
##' M. Binois, D. Ginsbourger, O. Roustant (2018), On the choice of the low-dimensional domain for global optimization via random embeddings, arXiv:1704.05318 \cr \cr
##' M. Binois (2015), Uncertainty quantification on Pareto fronts and high-dimensional strategies in Bayesian optimization, with applications in multi-objective automotive design, PhD thesis, Mines Saint-Etienne.
##' @export
##' @examples
##' ## Example of designs in Z
##' set.seed(42)
##' d <- 2; D <- 5
##'
##' A <- selectA(d, D, type = 'optimized')
##' size <- sqrt(D) # box size of Z
##' ntest <- 10000
##' Z <- size * (2 * matrix(runif(ntest * d), ntest, d) - 1)
##'
##' inZ <- testZ(Z, t(A))
##' colors <- rep('black', ntest)
##' colors[inZ] <- 'green'
##' plot(Z, col = colors, pch = 20, cex = 0.5)
##'
##' p <- 20
##' designs1 <- designZ(p, t(A), size)
##' designs2 <- designZ(p, t(A), size, type = 'LHS')
##'
##' points(designs1, col = 'red', pch = 20)
##' points(designs2, col = 'blue', pch = 20)
##' legend("topright", legend = c("LHS", "Maximin LHS"), col = c("blue", "red"), pch = c(20,20))
##'
designZ <- function(p, pA, bxsize, type = "maximin"){
  d <- nrow(pA)
  n <- 0
  designs <- NULL
  ind <- NULL
  while(n < p){
    #initialisation
    if(is.null(designs)){
      if(type == 'maximin')
        designs <- maximinLHS(p, d)
      if(type == 'LHS')
        designs <- randomLHS(p, d)
      if(type == 'unif')
        designs <- matrix(runif(p*d), p)
      designs <- 2 * bxsize * designs - bxsize
      ind <- testZ(designs, pA)
    }else{
      # Then add points uniformly
      ndesigns <- 2 * bxsize * matrix(runif((p - n)*d), p - n) - bxsize
      designs <- rbind(designs, ndesigns)
      ind <- c(ind, testZ(ndesigns, pA))
      bxsize <- bxsize/2 # to avoid taking too much time for this step
    }
    
    n <- sum(ind)
  }
  return(designs[ind,])
}



