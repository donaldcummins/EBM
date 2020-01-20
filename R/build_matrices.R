BuildABQ <- function(C, kappa, sigma_xi) {
  k <- length(C)
  A <- matrix(0, k, k)
  for (i in 1:(k - 1)) {
    A[i, i] <- -(kappa[i] + kappa[i + 1])/C[i]
    A[i, i + 1] <- kappa[i + 1]/C[i]
    A[i + 1, i] <- kappa[i + 1]/C[i + 1]
  }
  A[k, k] <- -kappa[k]/C[k]
  B <- matrix(0, k)
  B[1] <- 1/C[1]
  Q <- matrix(0, k, k)
  Q[1, 1] <- (sigma_xi/C[1])^2
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

BuildMatrices <- function(C, kappa, sigma_xi) {
  ABQ <- BuildABQ(C, kappa, sigma_xi)
  AdBdQd <- with(ABQ, BuildAdBdQd(A, B, Q))
  Gamma0 <- with(AdBdQd, BuildGamma0(Ad, Qd))
  return(c(ABQ, AdBdQd, Gamma0))
}


