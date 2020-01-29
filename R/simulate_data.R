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
  y[, 1] <- t(chol(Gamma0)) %*% rnorm(k)
  for (i in 1:(n - 1)) {
    y[, i + 1] <- Ad %*% y[, i] + t(chol(Qd)) %*% rnorm(k)
  }
  return(y)
}

SimStep <- function(Ad, Bd, Qd, Gamma0, F_4xCO2, n) {
  signal <- StepResponse(Ad, Bd, F_4xCO2, n)
  noise <- SimNoise(Ad, Qd, Gamma0, n)
  return(signal + noise)
}

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



