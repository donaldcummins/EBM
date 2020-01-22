Transform <- function(p) {
  log(unlist(p))
}

BackTransform <- function(par) {
  par <- exp(unname(par))
  if (length(par) == 9) {
    p <- list(
      gamma = par[1],
      C = par[2:3],
      kappa = par[4:5],
      epsilon = par[6],
      sigma_eta = par[7],
      sigma_xi = par[8],
      F_4xCO2 = par[9]
    )
    return(p)
  } else if (length(par) == 11) {
    p <- list(
      gamma = par[1],
      C = par[2:4],
      kappa = par[5:7],
      epsilon = par[8],
      sigma_eta = par[9],
      sigma_xi = par[10],
      F_4xCO2 = par[11]
    )
  } else {
    stop("k must be two or three")
  }
}





