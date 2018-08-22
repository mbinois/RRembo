##' Warping psi (original low dimensional space) 
##' @title Warping psi for points in Y
##' @param y matrix of low dimensional coordinates, one point per row
##' @param A random embedding matrix
##' @export
##' @references 
##' M. Binois, D. Ginsbourger, O. Roustant (2015), A Warped Kernel Improving Robustness in Bayesian Optimization Via Random Embeddings, Learning and Intelligent Optimization, Springer \cr \cr
##' @examples
##' set.seed(42)
##' d <- 2; D <- 5
##'
##' A <- selectA(d, D, type = 'optimized')
##' size <- 5 # box size of Y
##' ntest <- 10000
##' Y <- size * (2 * matrix(runif(ntest * d), ntest, d) - 1)
##' Z <- Psi_Y(Y, A)
##' plot(Z)
Psi_Y <- function(y, A){
  if(is.null(dim(y)))
    y <- matrix(y, nrow = 1)
  
  # pA <- tcrossprod(A)
  
  px <- randEmb(y, A)
  
  # papx <- tcrossprod(pA, px)
  papx <- tcrossprod(A, px %*% A)
  
  pivot <- Az <- papx
  
  for(i in 1:nrow(y)){
    if(max(abs(papx[,i])) < 1){
    }else{
      pivot[,i] <- papx[,i]/max(abs(papx[,i]))
      tmp <- distance(pivot[,i], px[i,])
      tmp2 <- sqrt(sum(pivot[,i]^2))
      Az[,i] <- pivot[,i]* (tmp2 + tmp)/tmp2
    }
  }
  
  return(crossprod(Az, A))
}

##' Warping psi for the alternative mapping procedure
##' @title Warping psi for points in the zonotope Z
##' @param z matrix of low dimensional coordinates, one point per row. These points must belong to Z, that can be checked by testZ
##' @param A random embedding matrix
##' @param eps to avoid some numerical issues with quadratic programming
##' @references 
##' M. Binois, D. Ginsbourger, O. Roustant (2018), On the choice of the low-dimensional domain for global optimization via random embeddings, arXiv:1704.05318 \cr
##' @export
##' @examples
##' # Comparison with between results starting from Y and from Z
##' set.seed(42)
##' d <- 2; D <- 5
##'
##' A <- selectA(d, D, type = 'optimized')
##' sizeY <- 5 # box size of Y
##' sizeZ <- sqrt(D) # box size of Z
##' ntest <- 10000
##' Y <- sizeY * (2 * matrix(runif(ntest * d), ntest, d) - 1)
##' X <- randEmb(Y, A)
##' Z <- ortProj(X, t(A))
##'
##' PsisY <- Psi_Y(Y, A)
##' ind <- testZ(Z, t(A), eps = 1e-8) ## Should be all TRUE, up to numerical error
##' PsisZ <- Psi_Z(Z[ind,], A)
##' plot(PsisY)
##' points(PsisZ, col = 'red', pch = 20)
##' print(max(abs(PsisY - PsisZ)))
Psi_Z <- function(z, A, eps = 1e-6){
  if(is.null(dim(z)))
    z <- matrix(z, nrow = 1)
  
  papx <- tcrossprod(A, z)
  
  px <- mapZX(z, A, eps = eps)
  
  pivot <- Az <- papx
  
  for(i in 1:nrow(z)){
    if(max(abs(papx[,i])) < 1){
    }else{
      pivot[,i] <- papx[,i]/max(abs(papx[,i]))
      tmp <- distance(pivot[,i], px[i,])
      tmp2 <- sqrt(sum(pivot[,i]^2))
      Az[,i] <- pivot[,i]* (tmp2 + tmp)/tmp2
    }
  }
  
  return(crossprod(Az, A))
}


##' Warping psi for points in Y, non normalized A
##' @param y matrix of low dimensional coordinates, one point per row
##' @param A random embedding matrix with non-othogonal columns
##' @param pA optional projection matrix onto Ran(A)
##' @param invA optional pseudo inverse of A
##' @export
##' @examples
##' set.seed(42)
##' d <- 2; D <- 5
##' library(MASS)
##' A <- selectA(d, D, type = 'standard')
##' size <- 5 # box size of Y
##' ntest <- 10000
##' Y <- size * (2 * matrix(runif(ntest * d), ntest, d) - 1)
##' Z <- Psi_Y(Y, A)
##' plot(Z)
##' Z2 <- Psi_Y_nonort(Y, A)
##' points(Z2, col = 'red')
##' library(far)
##' B <- orthonormalization(A, basis = FALSE)
##' Z3 <- Psi_Y(Y, B)
##' points(Z3, pch = 20)
##' @importFrom MASS ginv
##' @seealso  \code{\link[RemboIV]{Psi_Y}} if \code{A} has orthogonal columns
Psi_Y_nonort <- function(y, A, pA = NULL, invA = NULL){
  if(is.null(dim(y)))
    y <- matrix(y, nrow = 1)
  
  if(is.null(pA)) 
    pA <- A %*% ginv(t(A) %*% A) %*% t(A)
  
  if(is.null(invA)) invA <- ginv(A)

  px <- randEmb(y, A)

  papx <- t(pA %*% t(px))

  pivot <- Az <- papx

  for(i in 1:nrow(y)){
    if(max(abs(papx[i,])) < 1){
    }else{
      pivot[i,] <- papx[i,]/max(abs(papx[i,]))
      tmp <- distance(pivot[i,], px[i,])
      tmp2 <- sqrt(sum(pivot[i,]^2))
      Az[i,] <- pivot[i, ]* (tmp2 + tmp)/tmp2
    }
  }

  return(tcrossprod(Az, invA))
}

distance <- function(x1,x2){
  return(sqrt(sum((x1-x2)^2)))
}