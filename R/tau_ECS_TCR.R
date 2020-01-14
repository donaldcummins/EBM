Timescales <- function(mle) {
  p <- BackTransform(mle)
  -1/eigen(p$A)$values
}

ECS <- function(mle) {
  p <- BackTransform(mle)
  with(p, G/kappa[1])/2
}

TCR <- function(mle) {
  p <- BackTransform(mle)
  with(p, TransientResponseAnalytic(A, kappa, G, 70))[1,70]
}
