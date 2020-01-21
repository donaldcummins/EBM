KalmanNegLogLik <- function(par, dataset, k=2) {

  # back-transform parameters
  p <- BackTransform(par)

  # build matrices
  m <- with(p, BuildMatrices(gamma, C, kappa, epsilon, sigma_eta, sigma_xi))

  # run kalman filter
  kf <- with(c(p, m), KalmanFilter(Ad, Bd, Qd, Gamma0, Cd, F_4xCO2, dataset))

  # return negative log-likelihood
  return(-kf$logLik)

}
