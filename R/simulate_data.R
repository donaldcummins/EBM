StepResponse <- function(Ad, Bd, G, n) {
  nboxes <- nrow(Ad)
  x <- matrix(0, nboxes, n + 1)
  for (i in 1:n) {
    x[, i + 1] <- Ad %*% x[, i] + Bd * G
  }
  return(x[, -1])
}

SimNoise <- function(Ad, Qd, Gamma0, n) {
  nboxes <- nrow(Ad)
  y <- matrix(0, nboxes, n)
  y[, 1] <- t(chol(Gamma0)) %*% rnorm(nboxes)
  for (i in 1:(n - 1)) {
    y[, i + 1] <- Ad %*% y[, i] + t(chol(Qd)) %*% rnorm(nboxes)
  }
  return(y)
}

SimStep <- function(Ad, Bd, Qd, Gamma0, G, n) {
  signal <- StepResponse(Ad, Bd, G, n)
  noise <- SimNoise(Ad, Qd, Gamma0, n)
  return(signal + noise)
}

SimData <- function(C, kappa, sigma, G, epsilon, n) {

  # build matrices
  mats <- BuildABQ(C, kappa, sigma)
  dmats <- BuildAdBdQd(mats$A, mats$B, mats$Q)
  Ad <- dmats$Ad
  Bd <- dmats$Bd
  Qd <- dmats$Qd
  Gamma0 <- BuildGamma0(Ad, Qd)

  # generate data
  temp <- SimStep(Ad, Bd, Qd, Gamma0, G, n)[1, ]
  flux <- G - kappa[1]*temp + rnorm(n, sd = epsilon)

  # return simulated dataset
  return(list(temp = temp, flux = flux))

}

TransientResponse <- function(Ad, Bd, G, n) {
  nboxes <- nrow(Ad)
  alpha <- G/log(4)
  G <- alpha*log(1.01^(1:n))
  x <- matrix(0, nboxes, n + 1)
  for (i in 1:n) {
    x[, i + 1] <- Ad %*% x[, i] + Bd * G[i]
  }
  return(x[, -1])
}



