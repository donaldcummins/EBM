KalmanNegLogLik <- function(par, temp, flux) {

  # back-transform parameters
  par <- exp(par)

  # extract parameters
  nboxes <- as.integer((length(par) - 3)/2)
  n <- length(temp)
  C <- par[1:nboxes]
  kappa <- par[(nboxes + 1):(2*nboxes)]
  sigma <- par[2*nboxes + 1]
  G <- par[2*nboxes + 2]
  epsilon <- par[2*nboxes + 3]

  # build matrices
  m <- BuildMatrices(C, kappa, sigma)

  # run kalman filter
  kf <- KalmanFilter(m$Ad, m$Bd, m$Qd, m$Gamma0, G, temp)

  # log-likelihood for fluxes
  fluxLogLik <- sum(dnorm(flux, G - kappa[1]*temp, epsilon, log = TRUE))

  # return negative log-likelihood
  return(-(kf$logLik + fluxLogLik))

}
