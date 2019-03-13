BuildABQ <- function(C, kappa, sigma) {
  nboxes <- length(C)
  A <- matrix(0, nboxes, nboxes)
  for (i in 1:(nboxes - 1)) {
    A[i, i] <- -(kappa[i] + kappa[i + 1])/C[i]
    A[i, i + 1] <- kappa[i + 1]/C[i]
    A[i + 1, i] <- kappa[i + 1]/C[i + 1]
  }
  A[nboxes, nboxes] <- -kappa[nboxes]/C[nboxes]
  B <- matrix(0, nboxes)
  B[1] <- 1/C[1]
  Q <- matrix(0, nboxes, nboxes)
  Q[1, 1] <- (sigma/C[1])^2
  return(list(A = A, B = B, Q = Q))
}

BuildAdBdQd <- function(A, B, Q) {
  nboxes <- nrow(A)
  Ad <- expm::expm(A)
  Bd <- solve(A, (Ad - diag(nboxes)) %*% B)
  H <- rbind(cbind(-A, Q),
             cbind(matrix(0, nboxes, nboxes), t(A)))
  G <- expm::expm(H)
  Qd <- t(G[(nboxes + 1):(2*nboxes), (nboxes + 1):(2*nboxes)]) %*%
    G[1:nboxes, (nboxes + 1):(2*nboxes)]
  return(list(Ad = Ad, Bd = Bd, Qd = Qd))
}

BuildGamma0 <- function(Ad, Qd) {
  nboxes <- nrow(Ad)
  Gamma0 <- solve(diag(nboxes^2) - kronecker(Ad, Ad), as.vector(Qd))
  return(matrix(Gamma0, nboxes))
}

BuildMatrices <- function(C, kappa, sigma) {
  mats <- BuildABQ(C, kappa, sigma)
  dmats <- BuildAdBdQd(mats$A, mats$B, mats$Q)
  Gamma0 <- BuildGamma0(dmats$Ad, dmats$Qd)
  return(list(A = mats$A, B = mats$B, Q = mats$Q,
              Ad = dmats$Ad, Bd = dmats$Bd, Qd = dmats$Qd,
              Gamma0 = Gamma0))
}


