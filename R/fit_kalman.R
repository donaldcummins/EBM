FitKalman <- function(inits, T1, N, maxeval = 100000) {

  # transform to optimization domain
  inits <- Transform(inits)
  dataset <- rbind(T1, N)

  # minimize negative log-likelihood
  opt <- nloptr::bobyqa(
    x0 = inits,
    fn = KalmanNegLogLik,
    dataset = dataset,
    control = list(maxeval = maxeval))

  # check convergence
  if (opt$iter == maxeval) {
    stop("Iterations exhausted.")
  }
  if (opt$convergence < 0) {
    stop("Convergence failure!")
  }
  cat("Success! Convergence achieved in", opt$iter, "iterations.\n")

  # extract parameter estimates
  mle <- opt$par

  # calculate covariance matrix
  vcov <- solve(numDeriv::hessian(
    func = KalmanNegLogLik,
    x = mle,
    dataset = dataset
  ))

  # calculate standard errors
  se <- sqrt(diag(vcov))

  # calculate confidence intervals
  confint <- rbind(
    mle - qnorm(0.975)*se,
    mle + qnorm(0.975)*se
  )

  # calculate AIC
  AIC <- 2*length(opt$par) + 2*opt$value

  # back-transform parameters
  p <- BackTransform(mle)

  # build matrices at mle
  m <- with(p, BuildMatrices(gamma, C, kappa, epsilon, sigma_eta, sigma_xi))

  # calculate characteristic timescales
  tau <- with(m, -1/eigen(A[-1, -1])$values)

  # run Kalman filter at mle
  kf <- with(c(p, m), KalmanFilter(Ad, Bd, Qd, Gamma0, Cd, F_4xCO2, dataset))

  # step response of fitted model
  step <- with(c(p, m), StepResponse(Ad, Bd, F_4xCO2, 150))

  # transient response
  transient <- with(c(p, m), TransientResponseAnalytic(A, kappa, F_4xCO2, 150))

  # impulse response
  impulse <- with(c(p, m), ImpulseResponseAnalytic(A, kappa, 0:150))

  # ECS and TCR
  ECS <- with(p, 0.5*F_4xCO2/kappa[1])
  TCR <- transient[1, 70]

  # return output
  return(list(
    mle = mle,
    vcov = vcov,
    se = se,
    confint = confint,
    AIC = AIC,
    p = p,
    m = m,
    tau = tau,
    T1 = T1,
    N = N,
    kf = kf,
    step = step,
    transient = transient,
    impulse = impulse,
    ECS = ECS,
    TCR = TCR
  ))

}






