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





