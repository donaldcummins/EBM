FitAllTwoBox <- function(CMIP5) {

  # extract model names
  models <- names(CMIP5)

  # list to store fits
  twoBoxFits <- list()

  # two-box starting values
  inits <- log(c(7.3, 106, 1.13, 0.73, 0.7, 7.41, 0.2))

  # fit models
  for (model in models) {
    cat("Fitting model", model, "...\n")
    twoBoxFits[[model]] <-
      FitKalman(inits, CMIP5[[model]]$temp, CMIP5[[model]]$flux)
  }

  # return list of fits
  return(twoBoxFits)

}

FitAllGeoffroy <- function(CMIP5) {

  # extract model names
  models <- names(CMIP5)

  # list to store fits
  geoffroyFits <- list()

  # fit models
  for (model in models) {
    cat("Fitting model", model, "...\n")
    geoffroyFits[[model]] <-
      FitGeoffroy(CMIP5[[model]]$temp, CMIP5[[model]]$flux)
  }

  # return list of fits
  return(geoffroyFits)

}

FitAllBest <- function(CMIP5) {

  # extract model names
  models <- names(CMIP5)

  # list to store fits
  bestFits <- list()

  # starting values
  inits2 <- log(c(7.3, 106, 1.13, 0.73, 0.7, 7.41, 0.2))
  inits3 <- log(c(5, 30, 150, 1.13, 3.5, 1, 0.7, 7.41, 0.2))
  initsCCCMA <- log(c(3, 10, 20, 40, 1, 2, 1, 0.5, 0.6, 7.6, 0.3))
  initsCNRM <- log(c(3, 10, 30, 60, 1, 3.5, 3.5, 1, 0.4, 7.3, 0.2))
  initsIPSL <- log(c(1, 3, 16, 94, 0.79, 2.5, 2.2, 0.65, 0.5, 6.39, 0.25))

  # fit models
  for (model in models) {
    cat("Fitting model", model, "...\n")
    if (model == "GISS") {
      bestFits[[model]] <-
        FitKalman(inits2, CMIP5[[model]]$temp, CMIP5[[model]]$flux)
    } else if (model == "CCCMA") {
      bestFits[[model]] <-
        FitKalman(initsCCCMA, CMIP5[[model]]$temp, CMIP5[[model]]$flux)
    } else if (model == "CNRM") {
      bestFits[[model]] <-
        FitKalman(initsCNRM, CMIP5[[model]]$temp, CMIP5[[model]]$flux)
    } else if (model == "IPSL") {
      bestFits[[model]] <-
        FitKalman(initsIPSL, CMIP5[[model]]$temp, CMIP5[[model]]$flux)
    } else {
      bestFits[[model]] <-
        FitKalman(inits3, CMIP5[[model]]$temp, CMIP5[[model]]$flux)
    }
  }

  # return list of fits
  return(bestFits)

}

TwoBoxTable <- function(twoBoxFits) {

  # extract model names
  models <- names(twoBoxFits)

  # dataframe to store fits
  twoBoxTable <- as.data.frame(matrix(nrow = 16, ncol = 10))
  names(twoBoxTable) <- c("model", "C1", "C2", "kappa1", "kappa2", "sigma",
                          "G", "epsilon", "tauf", "taus")
  twoBoxTable$model <- models

  # fill table
  for (model in models) {
    twoBoxTable[twoBoxTable$model == model, 2:8] <-
      exp(twoBoxFits[[model]]$mle)
    twoBoxTable[twoBoxTable$model == model, 9:10] <-
      -1/eigen(BackTransform(twoBoxFits[[model]]$mle)$A)$values
  }

  # return dataframe of fits
  return(twoBoxTable)

}

GeoffroyTable <- function(geoffroyFits) {

  # extract model names
  models <- names(geoffroyFits)

  # dataframe to store fits
  geoffroyTable <- as.data.frame(matrix(nrow = 16, ncol = 8))
  names(geoffroyTable) <- c("model", "C1", "C2", "kappa1", "kappa2", "G",
                            "tauf", "taus")
  geoffroyTable$model <- models

  # fill table
  for (model in models) {
    geoffroyTable[geoffroyTable$model == model, -1] <-
      exp(geoffroyFits[[model]])
  }

  # return dataframe of fits
  return(geoffroyTable)

}

BestTable <- function(bestFits) {

  # extract model names
  models <- names(bestFits)

  # dataframe to store fits
  bestTable <- as.data.frame(matrix(nrow = 16, ncol = 12))
  names(bestTable) <- c("model", "C1", "C2", "C3", "C4", "kappa1", "kappa2",
                       "kappa3", "kappa4", "sigma", "G", "epsilon")
  bestTable$model <- models

  # fill table
  for (model in models) {
    if (model == "GISS") {
      bestTable[bestTable$model == model, -c(1, 4, 5, 8, 9)] <-
        exp(bestFits[[model]]$mle)
    } else if (model == "CCCMA") {
      bestTable[bestTable$model == model, -1] <-
        exp(bestFits[[model]]$mle)
    } else if (model == "CNRM") {
      bestTable[bestTable$model == model, -1] <-
        exp(bestFits[[model]]$mle)
    } else if (model == "IPSL") {
      bestTable[bestTable$model == model, -1] <-
        exp(bestFits[[model]]$mle)
    } else {
      bestTable[bestTable$model == model, -c(1, 5, 9)] <-
        exp(bestFits[[model]]$mle)
    }
  }

  # return dataframe of fits
  return(bestTable)

}

# predict transient response given vector of parameters on log scale
FilterTCR <- function(par, temp) {
  par <- exp(par[!is.na(par)])
  if (length(par) == 7) {
    C <- par[1:2]
    kappa <- par[3:4]
    sigma <- par[5]
    G <- par[6]
  } else if (length(par) == 9) {
    C <- par[1:3]
    kappa <- par[4:6]
    sigma <- par[7]
    G <- par[8]
  } else if (length(par) == 11) {
    C <- par[1:4]
    kappa <- par[5:8]
    sigma <- par[9]
    G <- par[10]
  }
  ABQ <- BuildABQ(C, kappa, sigma) # sigma not needed here
  AdBdQd <- BuildAdBdQd(ABQ$A, ABQ$B, ABQ$Q)
  Gamma0 <- BuildGamma0(AdBdQd$Ad, AdBdQd$Qd)
  kf <- KalmanFilterTCR(AdBdQd$Ad, AdBdQd$Bd, AdBdQd$Qd, Gamma0, G, temp)
  return(kf)
}

TCRTable <- function(TCR, twoBoxFits, bestFits) {

  # extract model names
  models <- names(twoBoxFits)

  # dataframe for negative log-likelihoods
  logLikTable <- as.data.frame(matrix(nrow = 16, ncol = 3))
  names(logLikTable) <- c("model", "two-Box", "best")
  logLikTable$model <- models

  # calculate TCR log-likelihood
  for (model in models) {
    temp <- TCR[[model]]$temp
    twoBoxNegLogLik <- -FilterTCR(twoBoxFits[[model]]$mle, temp)$logLik
    bestNegLogLik <- -FilterTCR(bestFits[[model]]$mle, temp)$logLik
    logLikTable[logLikTable$model == model, -1] <-
      c(twoBoxNegLogLik, bestNegLogLik)
  }

  # return dataframe of negative log-likelihoods
  return(logLikTable)

}



























