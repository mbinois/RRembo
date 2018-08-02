##' Select a random embedding matrix
##' @param d small dimension
##' @param D high dimension
##' @param type method of random sampling of coefficients or selection procedure, one of
##' \itemize{
##' \item 'Gaussian' for standard Gaussian i.i.d. coefficients and orthonormalization
##' \item 'isotropic' for random points on the d-sphere
##' \item 'optimized' for known optimal solutions (e.g., d = 2) or using a potential
##' \item 'standard' for the original REMBO iid random matrix
##' }
##' @param control list to be passed to \code{\link[stats]{optim}} in the \code{optimized} case (d > 2)
##' @return randomly selected matrix with orthogonal columns and normalized rows (except for standard)
##' @export
##' @importFrom far orthonormalization
##' @importFrom stats rnorm
##' @author Mickael Binois
##' @references M. Binois (2015), Uncertainty quantification on Pareto fronts and high-dimensional strategies in Bayesian optimization, with applications in multi-objective automotive design, PhD thesis, Mines Saint-Etienne.
##' @examples
##' ## Example of orthogonal projections
##' d <- 2; D <- 6
##' A1 <- selectA(d, D, type = 'Gaussian')
##' A2 <- selectA(d, D, type = 'isotropic')
##' A3 <- selectA(d, D, type = 'optimized')
##'
##' n <- 10000
##' size <- 10
##' Y <- size * (2 * matrix(runif(n * d), n) - 1)
##'
##' Z1 <- ortProj(randEmb(Y, A1), t(A1))
##' Z2 <- ortProj(randEmb(Y, A2), t(A2))
##' Z3 <- ortProj(randEmb(Y, A3), t(A3))
##'
##' par(mfrow = c(1, 3))
##' plot(Z1, asp = 1)
##' plot(Z2, asp = 1)
##' plot(Z3, asp = 1)
##'
##' par(mfrow = c(1, 1))
##'
selectA <- function(d, D, type = 'isotropic', control = list(n = 30, maxit = 100, maxit2 = 10)){
  if(is.null(control$n))
    control$n <- 10
  if(is.null(control$maxit))
    control$maxit <- 100
  if(is.null(control$maxit2))
    control$maxit2 <- 10

  if(type == 'Gaussian'){
    A <- matrix(rnorm(D * d), D, d)
    A <- orthonormalization(A, basis = FALSE)
  }

  if(type == 'standard'){
    A <- matrix(rnorm(D * d), D, d)
  }

  if(type == 'isotropic'){
    A <- matrix(rnorm(D * d), D, d)
    for(i in 1:control$n){
      A <- A/sqrt(rowSums(A^2))
      A <- orthonormalization(A, basis = FALSE)
    }
  }

  if(type == 'optimized'){
    A <- selectA_II_Sphere(d = d, D = D, display = FALSE, maxit = control$maxit)
    A <- orthonormalization(A, basis = FALSE)
    for(i in 1:control$n){
      A <- A/sqrt(rowSums(A^2))
      A <- selectA_II_Sphere(d = d, D = D, display = FALSE, A = A, maxit = control$maxit2)
      A <- orthonormalization(A, basis = FALSE)
    }
  }
  return(A)

}

##' Select a matrix A based on minimizing a potential function
##' @title SelectA
##' @param d embedding dimension
##' @param D high dimension
##' @param A optional starting matrix (d > 2)
##' @param init_dir optional start vector direction (d = 2)
##' @noRd
##' @importFrom numDeriv grad
##' @importFrom stats runif
##' @importFrom stats optim
##' @importFrom stats dist
selectA_II_Sphere <- function(d, D, A = NULL, init_dir = NULL, display = TRUE, maxit = 1000){
  if(d == 2){
    rr <- pi/D
    
    AA <- init_dir
    
    if(is.null(init_dir))
      AA <- rnorm(2)
    
    V <- runif(D-1)
    V <- as.numeric(V<0.5)*2-1
    
    Mrot <- matrix(c(cos(rr), -sin(rr), sin(rr), cos(rr)), byrow = T, nrow = 2)
    
    tmp <- AA
    
    for(i in 1:(D - 1)){
      tmp <- as.numeric(Mrot %*% tmp)
      AA <- rbind(AA, V[i]*tmp)
    }
    
    # normalisation des colonnes
    AA <- AA/sqrt(rowSums(AA^2))
    
  }else{
    
    AA <- A
    if(is.null(A)){
      AA <- matrix(rnorm(d*D), D, d)
    }
    
    Astart <- AA/sqrt(rowSums(AA^2))
    
    sol <- optim(par = as.vector(Astart), fn = potential_constraint, gr = gr_potential_constraint,
                 lowd = d, highD = D, method = "L-BFGS-B", control = list(maxit = maxit))
    
    AA <- matrix(sol$par, D, d)
  }
  return(AA)
}

potential <- function(A, lowd, highD){
  A <- matrix(A, highD, lowd)
  Ap <- rbind(A,-A)
  return(-min(dist(Ap)))
}

gr_potential <- function(A, lowd, highD){
  return(grad(potential, x = A, lowd = lowd, highD = highD, methods.args = list(r = 4)))
}

constraint <- function(A, lowd, highD){
  A <- matrix(A, highD, lowd)
  return(sum((rowSums(A^2) - 1)^2))
}

gr_constraint <- function(A, lowd, highD){
  A <- matrix(A, highD, lowd)
  return(grad(constraint, x = A, lowd = lowd, highD = highD, methods.args = list(r = 4)))
}

potential_constraint <- function(A, lowd, highD, pena = 1000){
  return(potential(A, lowd = lowd, highD = highD) + pena * constraint(A, lowd = lowd, highD = highD))
}

gr_potential_constraint <- function(A, lowd, highD, pena = 1000){
  return(grad(potential_constraint, x = A, method.args = list(r=4),
              lowd = lowd, highD = highD, pena = pena))
}




