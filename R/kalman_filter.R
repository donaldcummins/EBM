KalmanFilter <- function(Ad, Bd, Qd, Gamma0, Cd, F_4xCO2, dataset) {

  # get dimension of state vector
  k <- nrow(Ad)

  # set initial conditions (t=0)
  x0 <- c(F_4xCO2, rep(0, k - 1))

  # Kalman filter parameters
  a0 <- as.vector(Ad %*% x0 + Bd*F_4xCO2) # first observation one year in
  P0 <- Gamma0
  dt <- Bd*F_4xCO2
  ct <- matrix(0, 2, 1)
  Tt <- array(Ad, c(k, k, 1))
  Zt <- array(Cd, c(2, k, 1))
  HHt <- array(Qd, c(k, k, 1))
  GGt <- array(diag(1e-12, 2), c(2, 2, 1))
  yt <- dataset

  # run filter
  kf <- FKF::fkf(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt)

  # return FKF object
  return(kf)

}




