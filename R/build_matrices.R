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

BuildMatrices <- function(gamma, C, kappa, epsilon, sigma_eta, sigma_xi) {
  ABQ <- BuildABQ(gamma, C, kappa, epsilon, sigma_eta, sigma_xi)
  AdBdQd <- with(ABQ, BuildAdBdQd(A, B, Q))
  Gamma0 <- with(AdBdQd, BuildGamma0(Ad, Qd))
  Cd <- BuildCd(kappa, epsilon)
  return(c(ABQ, AdBdQd, Gamma0, Cd))
}


