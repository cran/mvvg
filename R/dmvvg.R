#' Calculate Matrix-Variate Variance Gamma Density
#'
#' Determines density of observations from a Matrix-variate variance gamma (MVVG) distribution, under the identifiability constraint set by [].
#'
#' MVVG samples are formulated through the normal variance-mean mixture \eqn{M + WA + \sqrt{W}Z}, where \eqn{W \sim Gamma(\gamma, \gamma)}.
#'
#' @param X \eqn{p \times q} observed matrix value
#' @param M \eqn{p \times q} location matrix
#' @param A \eqn{p \times q} skewness matrix
#' @param Sigma \eqn{p \times p} covariance matrix
#' @param Psi \eqn{q \times q} covariance matrix
#' @param gamma scalar mixing parameter
#' @param log returns log-likelihood if TRUE, default is FALSE.
#' @details
#' Gamma must be \eqn{>0}. Sigma and Psi must be positive definite covariance matrices.
#'
#' @return dmvvg returns the probability density corresponding to the inputted values and parameters.
#' @import nlme
#' @import psych
#' @export
#' @seealso \code{\link{rmvvg}}
#' @author Samuel Soon
#' @examples
#' M <- cbind(rep(1, 5), c(1, 0, 1, 0, 1))
#' A <- matrix(c(1,2), 5, 2, byrow = TRUE)
#' Sigma <- diag(5)
#' Psi <- matrix(c(4,2,2,3), 2, 2)
#' gamma <- 3
#'
#' X <- rmvvg(1, M, A, Sigma, Psi, gamma)[[1]]
#' dmvvg(X, M, A, Sigma, Psi, gamma)
dmvvg <- function(X, M, A, Sigma, Psi, gamma, log = FALSE) {
  if(gamma <= 0){
    stop("Gamma must be > 0 ")
  }
  Sigma_inv <- solve(Sigma)
  Psi_inv <- solve(Psi)
  rhofun <- tr(Sigma_inv %*% A %*% Psi_inv %*% t(A))
  deltafun <- tr(Sigma_inv %*% (X - M) %*% Psi_inv %*% t(X - M))
  n <- nrow(X)
  p <- ncol(X)

  t1 <- log(2) +
    gamma * log(gamma) +
    tr(Sigma_inv %*% (X - M) %*% Psi_inv %*% t(A))
  t2 <- (n*p/2) * log(2*pi) +
    (p/2) * log(det(Sigma)) +
    (n/2) * log(det(Psi)) +
    lgamma(gamma)
  t3pow <- (gamma/2 - n*p/4)
  t3base <- (log(deltafun) - log(rhofun + 2*gamma))
  t3 <- ifelse(t3pow == 0, 0,
               ifelse(t3base == -Inf, -Inf,
                      t3pow*t3base))
  t4 <- log(besselK(sqrt((rhofun + 2*gamma) * deltafun), abs(gamma - n*p/2)))
  ret <- t1 - t2 + t3 + t4
  if(!log){
    ret <- exp(ret)
  }
  ret
}
