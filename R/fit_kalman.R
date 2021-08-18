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

#' Fitting k-box models to abrupt 4xCO2 experiments
#'
#' \code{FitKalman} fits a k-box model to time series of global mean surface
#' temperature and top-of-the-atmosphere net downward radiative flux from an
#' abrupt 4xCO2 climate model experiment. Parameters are estimated using the
#' Kalman filter and maximum likelihood estimation. Estimated parameters are
#' returned in a list along with a suite of useful derived quantities.
#'
#' @param inits list of parameter starting values containing the following
#' elements:
#' \itemize{
#'   \item \code{gamma} - stochastic forcing correlation parameter;
#'   \item \code{C} - vector of box heat capacities;
#'   \item \code{kappa} - vector of heat transfer coefficients;
#'   \item \code{epsilon} - deep ocean heat uptake efficacy factor;
#'   \item \code{sigma_eta} - stochastic forcing standard deviation parameter;
#'   \item \code{sigma_xi} - standard deviation of stochastic temperature
#'   disturbances;
#'   \item \code{F_4xCO2} - effective radiative forcing due to CO2 quadrupling.
#' }
#' @param T1 time series of global mean surface temperature.
#' @param N time series of top-of-the-atmosphere net downward radiative flux.
#' @param alpha quadratic penalty applied to heat capacity of deep-ocean box.
#' @param maxeval maximum number of iterations in BOBYQA optimization algorithm.
#'
#' @return \code{FitKalman} returns a list containing the following elements:
#' \itemize{
#'   \item \code{mle} - vector of logged maximum likelihood parameter estimates;
#'   \item \code{vcov} - approximate covariance matrix of estimated parameters;
#'   \item \code{se} - vector of approximate standard errors of logged parameter
#'   estimates;
#'   \item \code{confint} - matrix with two rows containing lower and upper
#'   bounds of approximate 95 percent confidence intervals for logged
#'   parameters;
#'   \item \code{AIC} - Akaike's information criterion score;
#'   \item \code{p} - list of maximum likelihood parameter estimates;
#'   \item \code{m} - list containing matrix representation of fitted model;
#'   \item \code{tau} - vector containing characteristic timescales of fitted
#'   model;
#'   \item \code{T1} - time series of global mean surface temperature used in
#'   model fitting;
#'   \item \code{N} - time series of top-of-the-atmosphere net downward
#'   radiative flux used in model fitting;
#'   \item \code{kf} - fitted Kalman filter object (see \code{\link[FKF]{fkf}});
#'   \item \code{step} - matrix with k + 1 rows containing expected trajectories
#'   of fitted model state variables F, T1, ..., Tk under abrupt 4xCO2 forcing;
#'   \item \code{transient} - matrix with k rows containing expected
#'   trajectories of box temperatures under a one-percent-per-year increasing
#'   CO2 scenario;
#'   \item \code{impulse} - matrix with k rows containing expected
#'   trajectories of box temperatures under a unit-impulse forcing scenario;
#'   \item \code{ECS} - equilibrium climate sensitivity (equilibrium
#'   temperature increase after doubling CO2);
#'   \item \code{TCR} - transient climate response (temperature increase after
#'   70 years of one-percent-per-year increasing CO2).
#' }
#' @export
#' @seealso \code{\link{BuildMatrices}}, \code{\link{SimStepData}}.
#'
#' @examples
#' # load CMIP5 abrupt 4xCO2 runs
#' data("CMIP5")
#'
#' # set three-box model starting values
#' inits3 <- list(
#'   gamma = 2.0,
#'   C = c(5.0, 20.0, 100.0),
#'   kappa = c(1.0, 2.0, 1.0),
#'   epsilon = 1.0,
#'   sigma_eta = 0.5,
#'   sigma_xi = 0.5,
#'   F_4xCO2 = 5.0
#' )
#'
#' # fit three-box model to HadGEM2-ES
#' HadGEM2 <- with(CMIP5$MOHC, {
#'   FitKalman(inits = inits3, T1 = temp, N = flux)
#' })
#'
#' # print parameter estimates
#' print(HadGEM2$p)
FitKalman <- function(inits, T1, N, alpha = 0, maxeval = 1e+05) {

  # transform to optimization domain
  inits <- Transform(inits)
  dataset <- rbind(T1, N)

  # minimize negative log-likelihood
  opt <- nloptr::bobyqa(
    x0 = inits,
    fn = KalmanNegLogLik,
    dataset = dataset,
    alpha = alpha,
    control = list(maxeval = maxeval))

  # check convergence
  if (opt$iter == maxeval) {
    stop("Iterations exhausted.")
  }
  if (opt$convergence < 0) {
    stop("Convergence failure!")
  }
  cat("Success! Convergence achieved in", opt$iter, "iterations.\n")

  # extract parameter estimates
  mle <- opt$par

  # calculate covariance matrix
  vcov <- solve(numDeriv::hessian(
    func = KalmanNegLogLik,
    x = mle,
    dataset = dataset,
    alpha = alpha
  ))

  # calculate standard errors
  se <- sqrt(diag(vcov))

  # calculate confidence intervals
  confint <- rbind(
    mle - stats::qnorm(0.975)*se,
    mle + stats::qnorm(0.975)*se
  )

  # calculate AIC
  AIC <- 2*length(opt$par) + 2*opt$value

  # back-transform parameters
  p <- BackTransform(mle)

  # build matrices at mle
  m <- with(p, BuildMatrices(gamma, C, kappa, epsilon, sigma_eta, sigma_xi))

  # calculate characteristic timescales
  tau <- with(m, -1/eigen(A[-1, -1])$values)

  # run Kalman filter at mle
  kf <- with(c(p, m), KalmanFilter(Ad, Bd, Qd, Gamma0, Cd, F_4xCO2, dataset))

  # step response of fitted model
  step <- with(c(p, m), StepResponse(Ad, Bd, F_4xCO2, 150))

  # transient response
  transient <- with(c(p, m), TransientResponseAnalytic(A, kappa, F_4xCO2, 150))

  # impulse response
  impulse <- with(c(p, m), ImpulseResponseAnalytic(A, kappa, 0:150))

  # ECS and TCR
  ECS <- with(p, 0.5*F_4xCO2/kappa[1])
  TCR <- transient[1, 70]

  # return output
  return(list(
    mle = mle,
    vcov = vcov,
    se = se,
    confint = confint,
    AIC = AIC,
    p = p,
    m = m,
    tau = tau,
    T1 = T1,
    N = N,
    kf = kf,
    step = step,
    transient = transient,
    impulse = impulse,
    ECS = ECS,
    TCR = TCR
  ))

}






