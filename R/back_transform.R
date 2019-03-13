BackTransform <- function(par) {

  # back-transform parameters
  par <- exp(par)

  # extract parameters
  nboxes <- as.integer((length(par) - 3)/2)
  C <- par[1:nboxes]
  kappa <- par[(nboxes + 1):(2*nboxes)]
  sigma <- par[2*nboxes + 1]
  G <- par[2*nboxes + 2]
  epsilon <- par[2*nboxes + 3]

  # build matrices
  m <- BuildMatrices(C, kappa, sigma)

  # return list of output
  return(list(C = C, kappa = kappa, sigma = sigma, G = G, epsilon = epsilon,
              A = m$A, B = m$B, Q = m$Q, Ad = m$Ad, Bd = m$Bd, Qd = m$Qd,
              Gamma0 = m$Gamma0))

}




