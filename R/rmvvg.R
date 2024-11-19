#' Generate Matrix-Variate Variance Gamma Samples
#'
#' Generates random samples from the matrix-variate variance gamma (MVVG) distribution, under the identifiability constraint set by [].
#'
#' MVVG samples are formulated through the normal variance-mean mixture \eqn{M + WA + \sqrt{W}Z}, where \eqn{W \sim Gamma(\gamma, \gamma)}.
#'
#' @param n number of observations
#' @param M \eqn{p \times q} location matrix
#' @param A \eqn{p \times q} skewness matrix
#' @param Sigma \eqn{p \times p} covariance matrix
#' @param Psi \eqn{q \times q} covariance matrix
#' @param gamma scalar mixing parameter
#' @details
#' Gamma must be \eqn{>0}. Sigma and  Psi must be positive definite covariance matrices.
#'
#' @return rmvvg returns a list of random samples.
#' @import MixMatrix
#' @import stats
#' @export
#' @seealso \code{\link{dmvvg}}
#' @author Samuel Soon
#' @examples
#' M <- cbind(rep(1, 5), c(1, 0, 1, 0, 1))
#' A <- matrix(c(1,2), 5, 2, byrow = TRUE)
#' Sigma <- diag(5)
#' Psi <- matrix(c(4,2,2,3), 2, 2)
#' gamma <- 3
#'
#' rmvvg(2, M, A, Sigma, Psi, gamma)
rmvvg <- function(n, M, A, Sigma, Psi, gamma) {
  if(gamma <= 0){
    stop("Gamma must be > 0 ")
  }
  # get mixing parameter, dimensions
  w <- rgamma(n,gamma, gamma)
  mrow <- nrow(M)
  mcol <- ncol(M)

  # generate n mvn samples
  z <-
    rmatrixnorm(n,
                matrix(0, mrow, mcol),
                U = Sigma,
                V = Psi ,
                list = TRUE)

  # mix the variables together for a single mvvg obs, add rownames to help denote them later
  help <- function(z, w, A, M) {

    ret <- M + w * A + sqrt(w) * z
    rownames(ret) <- 1:nrow(ret)
    ret
  }

  # generate n mvvg observations
  ret <-
    mapply(help,
           z,
           w,
           MoreArgs = list(A = A, M = M),
           SIMPLIFY = FALSE)


  ret

}
