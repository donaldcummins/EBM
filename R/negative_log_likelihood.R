# Copyright (C) 2020  Donald Cummins
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

KalmanNegLogLik <- function(par, dataset, alpha) {

  # back-transform parameters
  p <- BackTransform(par)

  # build matrices
  m <- with(p, BuildMatrices(gamma, C, kappa, epsilon, sigma_eta, sigma_xi))

  # run kalman filter
  kf <- with(c(p, m), KalmanFilter(Ad, Bd, Qd, Gamma0, Cd, F_4xCO2, dataset))

  # calculate penalty
  penalty <- alpha*p$C[length(p$C)]^2

  # return (penalized) negative log-likelihood
  return(-kf$logLik + penalty)

}
