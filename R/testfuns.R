##' Test functions adapted for the Rembo problem, i.e., possibly with non-active variables. All are defined in the [0,1]^d hypercube.\cr \cr
##' Implemented functions:\cr
##' Sphere function (Euclidean distance to given optimum) (nD)
##' @title Test function
##' @param x vector of input location, or, for \code{sphere} and \code{giunta}, matrix specifying locations where the function is to be evaluated, one point per row.
##' @param x0 optimum in the effective space (same length as \code{ii})
##' @param ii effective dimensions indices
## ' @return euclidean distance between x0[,ii] and x[,ii]
##' @rdname Test_functions
##' @export
##' @examples
##' ## Sphere function example
##' set.seed(35)
##' n <- 101
##' d <- 2
##' D <- 3
##'
##' ygrid <- seq(-3,3, length.out = n)
##' Y <- as.matrix(expand.grid(ygrid, ygrid))
##' A <- selectA(d, D)
##' opt <- c(0.5,0.5)
##' fgrid <- sphere(randEmb(Y, A), x0 = opt, ii = c(1,2))
##' filled.contour(ygrid, ygrid, matrix(fgrid, n), color = terrain.colors)
##' 
sphere <- function(x, x0, ii){
  if(is.null(nrow(x)))
    x <- matrix(x, nrow = 1)
  return(rowSums((x0 - x[,ii])^2))
}


##' Branin function (2D)
## ' @title Branin function for high dim tests
## ' @param x where to test
## ' @param ii indices of x to consider
##' @export
##' @rdname Test_functions
##' @importFrom DiceKriging branin
##' @seealso \code{\link[DiceKriging]{branin}} and \code{\link[DiceKriging]{hartman6}} for the original \code{branin} and \code{hartman6} functions.
branin_mod <- function(x, ii = c(1,2)){
  if(is.null(nrow(x))){
    x <- matrix(x, nrow = 1)
  }
  return(apply(x[, c(ii[1], ii[2]), drop  = F], 1, branin))
}

##' Hartman6 function (6D)
## ' @title Hartman6 function for high dim tests
## ' @param x where to test
## ' @param ii indices of x to consider
##' @export
##' @importFrom DiceKriging hartman6
##' @rdname Test_functions
hartman6_mod <- function(x, ii = c(1,2,3,4,5,6)){
  if(is.null(nrow(x)))
    x <- matrix(x, nrow = 1)
  return(apply(x[,c(ii[1], ii[2], ii[3], ii[4], ii[5], ii[6]), drop = F], 1, hartman6))
}


##' Hartman6 function (Jones version, 6D)
## ' @title Hartman6 function for high dim tests (Jones Version)
## ' @param x where to test
## ' @param ii indices of x to consider
##' @export
##' @references 
##' D. R. Jones, M. Schonlau, W. Welch (1998),  Efficient global optimization of expensive black-box functions
##' @rdname Test_functions
hartman6_mod_log <- function(x, ii = c(1,2,3,4,5,6)){
  if(is.null(nrow(x)))
    x <- matrix(x, nrow = 1)
  return(-log(-apply(x[,c(ii[1], ii[2], ii[3], ii[4], ii[5], ii[6]), drop = FALSE], 1, hartman6)))
}


##' Cola function (17D)
##' @title Cola function for high dim tests
## ' @param x input vector of dimension 17 in [0, 1]
##' @export
##' @rdname Test_functions
##' @references
##' Madsen, Kaj, and Julius Zilinskas. Testing branch-and-bound methods for global optimization. IMM, Department of Mathematical Modelling, Technical Universityof Denmark, 2000.
##' @examples
##' # cola function
##' xstar <- c(0.1629765, 0.6627425, 0.51240525, 0.38952613, 0.39005, 0.52558138, 0.0894825,
##'            0.6063985, 0.06719375, 0.81655625, 0.38809425, 0.67624, 0.11579125,
##'            0.74532125, 0.12766, 0.39901887, 0.2887775)
##' cola(xstar) # 11.7464
##' 
cola <- function(x){

  # Solution in original space
  # xstar <- c(0.651906, 1.30194, 0.099242, -0.883791, -0.8796, 0.204651, -3.28414, 0.851188, -3.46245,
  #          2.53245, -0.895246, 1.40992, -3.07367, 1.96257, -2.97872, -0.807849, -1.68978)

  if(length(x) != 17) warning("Input dimension of x is not 17\n")

  u <- x
  x <- y <- rep(0, 10)
  x[2] <- u[1]*4
  x[3:10] <- u[seq(2, 17, by = 2)]*8 - 4
  y[3:10] <- u[seq(3, 17, by = 2)]*8 - 4

  r <- sqrt(outer(x, x, "-")^2 + outer(y, y, "-")^2)

  dmat <- matrix(c(1.27, 1.69, 2.04, 3.09, 3.20, 2.86, 3.17, 3.21, 2.38,
                   1.69, 1.43, 2.35, 3.18, 3.22, 2.56, 3.18, 3.18, 2.31,
                   2.04, 2.35, 2.43, 3.26, 3.27, 2.58, 3.18, 3.18, 2.42,
                   3.09, 3.18, 3.26, 2.85, 2.88, 2.59, 3.12, 3.17, 1.94,
                   3.20, 3.22, 3.27, 2.88, 1.55, 3.12, 1.31, 1.70, 2.85,
                   2.86, 2.56, 2.58, 2.59, 3.12, 3.06, 1.64, 1.36, 2.81,
                   3.17, 3.18, 3.18, 3.12, 1.31, 1.64, 3.00, 2.95, 2.56,
                   3.21, 3.18, 3.18, 3.17, 1.70, 1.36, 2.95, 1.32, 2.91,
                   2.38, 2.31, 2.42, 1.94, 2.85, 2.81, 2.56, 2.91, 2.97), 9)

  res <- sum((dmat[lower.tri(dmat, diag = T)] - r[lower.tri(r)])^2)

  return(res)
}

##' Giunta function (nD)
##' @title Giunta function
## ' @param x where to test in [0, 1]^D
## ' @param ii indices of x to consider
##' @export
##' @rdname Test_functions
##' @examples 
##' # Giunta function
##' xstar <- c(0.73366, 0.73366)
##' giunta(xstar) # 0.06447042
##' 
giunta <- function(x, ii = NULL){
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(ii)) ii <- 1:ncol(x)
  x <- x*2 - 1

  tmp <- sin(16/15*x[,ii, drop = F] - 1) + sin(16/15 * x[,ii, drop = F] - 1)^2 + 1/50 * sin(4*(16/15 * x[,ii, drop = F] -1))

  return(0.6 + rowSums(tmp))
}

##' Levy function (nD)
##' @title Levy function
##' @export
##' @rdname Test_functions
##' @examples 
##' # Levy function
##' xstar <- rep(0.55, 10)
##' levy(xstar) # 0
##' 
##' @details Levy function adapted from the code of Sonja Surjanovic and Derek Bingham
##' and available at https://www.sfu.ca/~ssurjano/levy.html.
levy <- function(x, ii = NULL){
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(ii)) ii <- 1:ncol(x)
  
  x <- x * 20 - 10
  
  d <- length(ii)
  w <- 1 + (x[,ii, drop = F] - 1)/4
  
  term1 <- (sin(pi*w[,1]))^2 
  term3 <- (w[,d]-1)^2 * (1+1*(sin(2*pi*w[,d]))^2)
  
  wi <- w[,1:(d-1), drop = F]
  sum <- rowSums((wi-1)^2 * (1+10*(sin(pi*wi+1))^2))
  
  y <- term1 + sum + term3
  return(y)
}
