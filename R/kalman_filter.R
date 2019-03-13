KalmanFilter <- function(Ad, Bd, Qd, Gamma0, G, temp) {

  # get number of boxes
  nboxes <- nrow(Ad)

  # Kalman filter parameters
  a0 <- as.vector(Bd*G)
  P0 <- Gamma0
  dt <- Bd*G
  ct <- matrix(0)
  Tt <- array(Ad, c(nboxes, nboxes, 1))
  Zt <- array(c(1, rep(0, (nboxes - 1))), c(1, nboxes, 1))
  HHt <- array(Qd, c(nboxes, nboxes, 1))
  GGt <- array(1e-12, c(1, 1, 1))
  yt <- matrix(temp, 1)

  # run filter
  kf <- FKF::fkf(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt)

  # return FKF object
  return(kf)

}

KalmanFilterTCR <- function(Ad, Bd, Qd, Gamma0, G, temp) {

  # get number of boxes
  nboxes <- nrow(Ad)
  alpha <- G/log(4)
  G <- alpha*log(1.01^(0:(length(temp) - 1)))

  # Kalman filter parameters
  a0 <- as.vector(Bd*G[1])
  P0 <- Gamma0
  dt <- Bd %*% G
  ct <- matrix(0)
  Tt <- array(Ad, c(nboxes, nboxes, 1))
  Zt <- array(c(1, rep(0, (nboxes - 1))), c(1, nboxes, 1))
  HHt <- array(Qd, c(nboxes, nboxes, 1))
  GGt <- array(1e-12, c(1, 1, 1))
  yt <- matrix(temp, 1)

  # run filter
  kf <- FKF::fkf(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt)

  # return FKF object
  return(kf)

}

KalmanFilterObs <- function(Ad, Bd, Qd, Gamma0, G, temp, se) {

  # get number of boxes and length of dataset
  n <- length(temp)
  nboxes <- nrow(Ad)

  # Kalman filter parameters
  a0 <- as.vector(Bd*G[1])
  P0 <- Gamma0
  dt <- Bd %*% G
  ct <- matrix(0)
  Tt <- array(Ad, c(nboxes, nboxes, 1))
  Zt <- array(c(1, rep(0, (nboxes - 1))), c(1, nboxes, 1))
  HHt <- array(Qd, c(nboxes, nboxes, 1))
  GGt <- array(se^2, c(1, 1, n))
  yt <- matrix(temp, 1)

  # run filter
  kf <- FKF::fkf(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt)

  # return FKF object
  return(kf)

}



