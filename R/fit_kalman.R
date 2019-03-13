FitKalman <- function(inits, temp, flux, maxeval = 100000) {

  # minimize negative log-likelihood
  opt <- nloptr::bobyqa(inits, KalmanNegLogLik,
                        temp = temp, flux = flux,
                        control = list(maxeval = maxeval))

  # check convergence
  if (opt$iter == maxeval) {
    stop("Iterations exhausted.")
  }
  if (opt$convergence < 0) {
    stop("Convergence failure!")
  }

  # extract parameter estimates
  mle <- opt$par

  # calculate covariance matrix
  vcov <- solve(numDeriv::hessian(KalmanNegLogLik, opt$par,
                                  temp = temp, flux = flux))

  # calculate standard errors
  se <- sqrt(diag(vcov))

  # calculate confidence intervals
  confint <- rbind(mle - qnorm(0.975)*se,
                   mle + qnorm(0.975)*se)

  # calculate AIC
  AIC <- 2*length(opt$par) + 2*opt$value

  # back-transform parameters
  b <- BackTransform(mle)

  # calculate fitted values
  fitted.values <- list(temp = StepResponse(b$Ad, b$Bd, b$G, length(temp)),
                       flux = b$G - b$kappa[1]*temp)

  # calculate residuals
  residuals <- list(temp = temp - fitted.values$temp[1, ],
                    flux = flux - fitted.values$flux)

  # run Kalman filter at mle
  kf <- KalmanFilter(b$Ad, b$Bd, b$Qd, b$Gamma0, b$G, temp)

  # return list of output on log scale
  return(list(mle = mle, se = se, confint = confint, vcov = vcov, AIC = AIC,
              temp = temp, flux = flux, kf = kf,
              fitted.values = fitted.values, residuals = residuals))

}






