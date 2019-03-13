FitGeoffroy <- function(temp, flux) {

  # regress flux on temperature
  reg1 <- lm(flux ~ temp)
  forcing <- coef(reg1)[1]
  lambda <- -coef(reg1)[2]
  Teq <- forcing/lambda
  epsilon <- sigma(reg1)

  # estimate as and taus
  reg2 <- lm(log(1-temp/Teq)[30:150] ~ c(30:150))
  as <- exp(coef(reg2)[1])
  taus <- -1/coef(reg2)[2]
  af <- 1 - as

  # estimate tauf
  # i <- c(1:(min(which(1-temp[1:10]/Teq-as*exp(-c(1:10)/taus) <= 0)-1, 10)))
  # tauf <- mean(i/(log(af)-log(1-temp[i]/Teq-as*exp(-i/taus))))
  # reg3 <- nls(temp ~ Teq*(af*(1-exp(-c(1:150)/tauf))+as*(1-exp(-c(1:150/taus)))),
  #             start = list(tauf = tauf))
  # tauf <- coef(reg3)
  lsobj <- function(tauf) {
    sum((temp - Teq*(af*(1-exp(-c(1:150)/tauf))+as*(1-exp(-c(1:150/taus)))))^2)
  }
  tauf <- optimize(lsobj, c(0,taus))$minimum

  # retrieve C, C0, kappa
  C <- lambda/(af/tauf + as/taus)
  C0 <- lambda*(tauf*af + taus*as) - C
  kappa <- C0/(tauf*as + taus*af)

  # return vector of estimated parameters on log scale
  return(unname(log(c(C, C0, lambda, kappa, forcing, epsilon))))

}



