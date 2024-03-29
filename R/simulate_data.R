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

StepResponse <- function(Ad, Bd, F_4xCO2, n) {
  k <- nrow(Ad)
  x <- matrix(0, k, n + 1)
  x[1, 1] <- F_4xCO2
  for (i in 1:n) {
    x[, i + 1] <- Ad %*% x[, i] + Bd * F_4xCO2
  }
  return(x[, -1])
}

StepResponseAnalytic <- function(A, kappa, F_4xCO2, n) {
  A <- A[-1, -1]
  ones <- rep(1, nrow(A))
  func <- function(t) {
    F_4xCO2/kappa[1] * (ones - expm::expAtv(A, ones, t)$eAtv)
  }
  return(sapply(1:n, func))
}

SimNoise <- function(Ad, Qd, Gamma0, n) {
  k <- nrow(Ad)
  y <- matrix(0, k, n)
  y[, 1] <- t(chol(Gamma0)) %*% stats::rnorm(k)
  for (i in 1:(n - 1)) {
    y[, i + 1] <- Ad %*% y[, i] + t(chol(Qd)) %*% stats::rnorm(k)
  }
  return(y)
}

SimStep <- function(Ad, Bd, Qd, Gamma0, F_4xCO2, n) {
  signal <- StepResponse(Ad, Bd, F_4xCO2, n)
  noise <- SimNoise(Ad, Qd, Gamma0, n)
  return(signal + noise)
}

#' Simulating k-box model abrupt 4xCO2 step responses
#'
#' \code{SimStepData} takes physical parameters of a k-box model and simulates
#' the first n years of its response to an abrupt 4xCO2 forcing.
#'
#' @param gamma stochastic forcing correlation parameter.
#' @param C vector of box heat capacities.
#' @param kappa vector of heat transfer coefficients.
#' @param epsilon deep ocean heat uptake efficacy factor.
#' @param sigma_eta stochastic forcing standard deviation parameter.
#' @param sigma_xi standard deviation of stochastic temperature disturbances.
#' @param F_4xCO2 effective radiative forcing due to CO2 quadrupling.
#' @param n length of simulation in years.
#'
#' @return \code{SimStepData} returns a matrix with two rows containing time
#' series of global mean surface temperature and top-of-the-atmosphere net
#' downward radiative flux respectively.
#' @export
#' @seealso \code{\link{BuildMatrices}}, \code{\link{FitKalman}}.
#'
#' @examples
#' # set physical parameters
#' parameters <- list(
#'   gamma = 2.2,
#'   C = c(7.0, 80.0),
#'   kappa = c(1.2, 0.75),
#'   epsilon = 1.2,
#'   sigma_eta = 0.54,
#'   sigma_xi = 0.72,
#'   F_4xCO2 = 7.4,
#'   n = 150
#' )
#'
#' # simulate step response
#' step <- with(parameters, {
#'   SimStepData(gamma, C, kappa, epsilon, sigma_eta, sigma_xi, F_4xCO2, n)
#' })
#'
#' # plot results
#' plot.ts(t(step), main = "Simulated Step Response")
SimStepData <- function(gamma, C, kappa, epsilon, sigma_eta, sigma_xi, F_4xCO2, n) {

  # build matrices
  m <- BuildMatrices(gamma, C, kappa, epsilon, sigma_eta, sigma_xi)

  # generate data
  x <- with(m, SimStep(Ad, Bd, Qd, Gamma0, F_4xCO2, n))
  y <- matrix(0, 2, n)
  for (i in 1:n) {
    y[, i] <- with(m, Cd %*% x[, i])
  }
  rownames(y) <- c("T1", "N")

  # return simulated dataset
  return(y)

}

TransientResponseAnalytic <- function(A, kappa, F_4xCO2, n) {
  A <- A[-1, -1]
  k <- nrow(A)
  ones <- rep(1, k)
  func <- function(t) {
    ts <- rep(t, k)
    log(1.01)/log(4) * F_4xCO2/kappa[1] *
      (ts - solve(A) %*% (expm::expm(A*t) - diag(k)) %*% ones)
  }
  return(sapply(1:n, func))
}

ImpulseResponseAnalytic <- function(A, kappa, t) {
  A <- A[-1, -1]
  ones <- rep(1, nrow(A))
  func <- function(t) {
    -1/kappa[1] * A %*% expm::expAtv(A, ones, t)$eAtv
  }
  return(sapply(t, func))
}



