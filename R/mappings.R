##' Random embedding of a low dimensional subspace on the centered hypercube, based on a matrix and convex projection
##' Random embedding mapping
##' @param y matrix of low dimensional coordinates, one point per row
##' @param A random embedding matrix
##' @return Matrix corresponding to convex projection onto [-1, 1]^D of A * y
##' @export
##' @references 
##' Z. Wang, F. Hutter, M. Zoghi, D. Matheson, N. de Freitas (2016), Bayesian Optimization in a Billion Dimensions via Random Embeddings, JAIR. \cr \cr
##' @examples
##' ## Small dimensional examples for illustration
##' d <- 1; D <- 2
##' set.seed(42)
##'
##' ntest <- 1000
##' size <- 1.5 # box size of Y
##' Y <- size * (2 * matrix(runif(ntest * d), ntest, d) - 1)
##' A <- matrix(rnorm(D *d), D, d)
##' X <- randEmb(Y, A)
##' plot(t(A %*% t(Y)), type = 'l', ylim = c(-1,1), xlim = c(-1, 1))
##' points(X, col = 'red', pch = 20)
##'
##' library(rgl)
##' d <- 2; D <- 3
##' Y <- size * (2 * matrix(runif(ntest * d), ntest, d) - 1)
##' A <- matrix(rnorm(D *d), D, d)
##' X <- randEmb(Y, A)
##' plot3d(t(A %*% t(Y)))
##' points3d(X, col = 'red', pch = ".")
randEmb <- function(y, A){
  if(is.null(nrow(y)))
    y <- matrix(y, nrow = 1)
  Xmap <- t(tcrossprod(A, y))
  Xmap <- pmin(Xmap, 1)
  Xmap <- pmax(Xmap, -1)
  return(Xmap)
}

##' Test if a given point is in U
##' @param y matrix corresponding to points to be tested
##' @param A random embedding matrix
##' @export
##' @examples
##' ## Identification of the U set
##' set.seed(42)
##' d <- 2; D <- 5
##'
##' A <- selectA(d, D)
##' size <- 15 # box size of Y
##' ntest <- 10000
##' Y <- size * (2 * matrix(runif(ntest * d), ntest, d) - 1)
##'
##' inU <- testU(Y, A)
##' colors <- rep('black', ntest)
##' colors[inU] <- 'green'
##' plot(Y, col = colors, pch = 20, cex = 0.5)
##'
testU <- function(y, A){
  d <- ncol(A)
  D <- nrow(A)
  if(is.null(nrow(y)))
    y <- matrix(y, nrow = 1)
  X <- randEmb(y, A)
  X[which(abs(X) < 1)] <- 0
  return(!rowSums(abs(X)) > (D - d))
}

##' Orthogonal projection on the subspace spanned by A
##' @title Orthogonal projection on Ran(A)
##' @param x matrix of high-dimensional coordinates, one point per row
##' @param pA matrix giving the coordinates of the orthogonal projection onto Ran(A)
##' @return z matrix of coordinates of the orthogonal projection onto Ran(A)
##' @details It is assumed that rows of pA are orthonormal, such that this linear transformation gives the coordinates on Ran(A) with an orthonormal basis
##' @export
##' @author Mickael Binois
##' @examples
##' ## Example of orthogonal projection
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
ortProj <- function(x, pA){
  if(is.null(nrow(x)))
    x <- matrix(x, nrow = 1)
  return(t(tcrossprod(pA, x)))
}


##' Map from zonotope Z to the convex projection of A Y (E)
##' @title Alternative mapping
##' @param z matrix of low-dimensional coordinates in the zonotope, one point per row
##' @param A random embedding matrix, such that \code{t(A)A = Id}
##' @param eps to avoid some numerical issues with quadratic programming
##' @param Aind,Amat optional matrices to be passed to \code{\link[quadprog]{solve.QP.compact}}
##' @details Numerical problems may occur, then the equality constraints are slightly relaxed.
##' @author Mickael Binois
##' @export
##' @importFrom quadprog solve.QP solve.QP.compact
##' @references 
##' M. Binois, D. Ginsbourger, O. Roustant (2018), On the choice of the low-dimensional domain for global optimization via random embeddings, arXiv:1704.05318 \cr
##' @examples
##' ## Example of forward - backward mapping
##' d <- 5; D <- 26
##' A <- selectA(d, D, type = 'Gaussian')
##'
##' n <- 10000
##' size <- 10
##' Y <- size * (2 * matrix(runif(n * d), n) - 1)
##' X <- randEmb(Y, A)
##' Z <- ortProj(X, t(A))
##' # Errors are catched and the problem slightly relaxed
##' Xback <- mapZX(Z, A)
##'
##' print(max(abs(Xback - X)))
##' print(mean(Xback - X))
mapZX <- function(z, A, eps = 1e-6, Amat = NULL, Aind = NULL){
  if(is.null(nrow(z)))
    z <- matrix(z, nrow = 1)
  return(t(apply(z, 1, mapzx_1, A = A, eps = eps)))
}

# for one vector
mapzx_1 <- function(z, A, eps = 1e-6, Aind = NULL, Amat = NULL){
  
  # Points in I need no optimization
  if(max(abs(A %*% z)) <= 1)
    return(A %*% z)
  
  D <- nrow(A)
  d <- ncol(A)
  
  if(is.null(Aind)) Aind <- cbind(matrix(c(D, 1:D), D + 1, d), rbind(rep(1, D * 2), c(1:D, 1:D), matrix(0, D-1, D*2)))
  if(is.null(Amat)) Amat <- cbind(A, matrix(rep(c(1, 0), times = c(1, D-1)), D, D), matrix(rep(c(-1, 0), times = c(1, D-1)), D, D))
  
  # tmp <- try(solve.QP(Dmat = diag(D), dvec = as.vector(A %*% z),
  #                     meq = d,
  #                     Amat = cbind(A, diag(D), -diag(D)),
  #                     bvec = c(z, rep(-1, D), rep(-1,D)), factorized = T)$solution)
  tmp <- try(solve.QP.compact(Dmat = diag(D), dvec = as.vector(A %*% z),
                      meq = d,
                      Amat = Amat,
                      Aind = Aind,
                      bvec = c(z, rep(-1, D), rep(-1,D)), factorized = T)$solution)

  if(length(tmp) == D){
    x <- tmp
  }else{
    tmp <- try(solve.QP(Dmat = diag(D), dvec = as.vector(A %*% z),
                        meq = 0,
                        Amat = cbind(A, -A, diag(D), -diag(D)),
                        bvec = c(z - eps, -z + eps,
                                 rep(-1, D), rep(-1,D)), factorized = T)$solution)
    if(length(tmp) == D){
      x <- tmp
    }else{
      z <- (1 - eps) * z
      tmp <- try(solve.QP(Dmat = diag(D), dvec = as.vector(A %*% z),
                          meq = d,
                          Amat = cbind(A, diag(D), -diag(D)),
                          bvec = c(z, rep(-1, D), rep(-1,D)), factorized = T)$solution)
    }
    if(length(tmp) == D){
      x <- tmp
    }
  }
  
  return(x)
}

##' Test if a given point is in Z
##' @param z matrix corresponding to points to be tested
##' @param pA matrix of orthogonal projection onto Ran(A) 
##' @param eps numerical precision on the result
##' @param bxsize radius of the sphere enclosing the zonotope
##' @export
##' @return ..
##' @importFrom linprog solveLP
##' @export
##' @examples
##' ## Identification of the Z set
##' set.seed(42)
##' d <- 2; D <- 5
##'
##' A <- selectA(d, D)
##' size <- 2 # box size of Y
##' ntest <- 10000
##' Z <- size * (2 * matrix(runif(ntest * d), ntest, d) - 1)
##'
##' inZ <- testZ(Z, t(A))
##' colors <- rep('black', ntest)
##' colors[inZ] <- 'green'
##' plot(Z, col = colors, pch = 20, cex = 0.5)
##'
##' ## Identification of a zone
##' X <- mapZX(Z[inZ,], A)
##' inZone1 <- which(abs(X[,1]) < 1)
##' points((Z[inZ,])[inZone1,], col ='red')
##'
testZ <- function(z, pA, eps = 1e-10, bxsize = sqrt(ncol(pA))){
  if(is.null(nrow(z)))
    z <- matrix(z, nrow = 1)
  tmp <- apply(z, 1, testz_1, pA = pA, bxsize = bxsize)
  
  return(tmp > (1 - eps))
}


## @param z potential point
## @param pA matrix with generators (ncol = nb of generators)
## @value linear programming  solution (if tmp$opt == 1 : c'est un sommet, sinon non)
## Solution from linear programming. If sol$opt <= 1, the point is in the zonotope
testz_1 <- function(x, pA, bxsize){
  n <- ncol(pA) # number of generators
  m <- nrow(pA)
  
  if(sqrt(sum(x^2)) > bxsize)
    return(0)
  
  if(all(abs(x %*% pA) <= 1)) return(1)
  
  tmp <- solveLP(cvec = c(1, rep(0, n)),
                 bvec = c(pA %*% rep(-1, n), rep(1, n)),
                 Amat = rbind(cbind(x, -2*pA), cbind(0, diag(n))),
                 const.dir = c(rep("==", m), rep("<=", n)),
                 maximum = TRUE,
                 lpSolve = TRUE)
  # tmp2 <- lp(direction = "max", objective.in = c(1, rep(0, n)),
  #            const.mat = rbind(cbind(x, -2*pA), cbind(0, diag(n))),
  #            const.rhs = c(pA %*% rep(-1, n), rep(1, n)),
  #            const.dir = c(rep("==", m), rep("<=", n)))
  
  return(tmp$opt)
}






