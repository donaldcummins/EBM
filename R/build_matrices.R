# Copyright (C) 2020  Donald Cummins
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

BuildABQ <- function(gamma, C, kappa, epsilon, sigma_eta, sigma_xi) {
  k <- length(C)
  if (k == 2) {
    A <- rbind(
      c(-gamma, 0, 0),
      c(1/C[1], -(kappa[1] + epsilon*kappa[2])/C[1], epsilon*kappa[2]/C[1]),
      c(0, kappa[2]/C[2], -kappa[2]/C[2])
    )
  } else if (k == 3) {
    A <- rbind(
      c(-gamma, 0, 0, 0),
      c(1/C[1], -(kappa[1] + kappa[2])/C[1], kappa[2]/C[1], 0),
      c(0, kappa[2]/C[2], -(kappa[2] + epsilon*kappa[3])/C[2], epsilon*kappa[3]/C[2]),
      c(0, 0, kappa[3]/C[3], -kappa[3]/C[3])
    )
  } else {
    stop("number of boxes k must be two or three")
  }
  B <- matrix(0, k + 1)
  B[1] <- gamma
  Q <- matrix(0, k + 1, k + 1)
  Q[1, 1] <- sigma_eta^2
  Q[2, 2] <- (sigma_xi/C[1])^2
  return(list(A = A, B = B, Q = Q))
}

BuildAdBdQd <- function(A, B, Q) {
  k <- nrow(A)
  Ad <- expm::expm(A)
  Bd <- solve(A, (Ad - diag(k)) %*% B)
  H <- rbind(
    cbind(-A, Q),
    cbind(matrix(0, k, k), t(A))
  )
  G <- expm::expm(H)
  Qd <- t(G[(k + 1):(2*k), (k + 1):(2*k)]) %*% G[1:k, (k + 1):(2*k)]
  return(list(Ad = Ad, Bd = Bd, Qd = Qd))
}

BuildGamma0 <- function(Ad, Qd) {
  k <- nrow(Ad)
  Gamma0 <- solve(diag(k^2) - kronecker(Ad, Ad), as.vector(Qd))
  Gamma0 <- matrix(Gamma0, k)
  return(list(Gamma0 = Gamma0))
}

BuildCd <- function(kappa, epsilon) {
  k <- length(kappa)
  if (k == 2) {
    Cd <- rbind(
      c(0, 1, 0),
      c(1, -kappa[1] + (1 - epsilon)*kappa[2], -(1 - epsilon)*kappa[2])
    )
  } else if (k == 3) {
    Cd <- rbind(
      c(0, 1, 0, 0),
      c(1, -kappa[1], (1 - epsilon)*kappa[3], -(1 - epsilon)*kappa[3])
    )
  } else {
    stop("k must be two or three")
  }
  return(list(Cd = Cd))
}

#' Building k-box model matrices
#'
#' \code{BuildMatrices} takes physical parameters of the k-box model
#' differential equations and constructs continuous- and discrete-time model
#' matrices.
#'
#' @param gamma stochastic forcing correlation parameter.
#' @param C vector of box heat capacities.
#' @param kappa vector of heat transfer coefficients.
#' @param epsilon deep ocean heat uptake efficacy factor.
#' @param sigma_eta stochastic forcing standard deviation parameter.
#' @param sigma_xi standard deviation of stochastic temperature disturbances.
#'
#' @return \code{BuildMatrices} returns a list containing k-box model matrices:
#' A, B, Q, Ad, Bd, Qd, Gamma0, Cd.
#' @export
#'
#' @examples
#' # set physical parameters
#' parameters <- list(
#'   gamma = 2.2,
#'   C = c(7.0, 80.0),
#'   kappa = c(1.2, 0.75),
#'   epsilon = 1.2,
#'   sigma_eta = 0.54,
#'   sigma_xi = 0.72
#' )
#'
#' # build matrices
#' matrices <- with(parameters, {
#'   BuildMatrices(gamma, C, kappa, epsilon, sigma_eta, sigma_xi)
#' })
#'
#' # print output
#' print(matrices)
BuildMatrices <- function(gamma, C, kappa, epsilon, sigma_eta, sigma_xi) {
  ABQ <- BuildABQ(gamma, C, kappa, epsilon, sigma_eta, sigma_xi)
  AdBdQd <- with(ABQ, BuildAdBdQd(A, B, Q))
  Gamma0 <- with(AdBdQd, BuildGamma0(Ad, Qd))
  Cd <- BuildCd(kappa, epsilon)
  return(c(ABQ, AdBdQd, Gamma0, Cd))
}


