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

# function to read in temperatures
ReadTemps <- function(models) {

  # initialize dataframe
  temps <- data.frame(matrix(nrow = 150, ncol = length(models)))
  names(temps) <- models

  # extract temperatures from files
  for (model in models) {
    temps[[model]] <- read.table(paste0("AB4_tas_",model,".txt"))[3:152,]
  }

  # return dataframe of temperatures
  return(temps)

}

# function to read in fluxes
ReadFluxes <- function(models) {

  # initialize dataframe
  fluxes <- data.frame(matrix(nrow = 150, ncol = length(models)))
  names(fluxes) <- models

  # extract fluxes from files
  for (model in models) {
    fluxes[[model]] <- read.table(paste0("AB4_rtot_",model,".txt"))[3:152,]
  }

  # return dataframe of temperatures
  return(fluxes)

}

# function to build CMIP5 list of dataframes
BuildCMIP5 <- function(models, temps, fluxes) {

  # build list of dataframes
  CMIP5 <- lapply(models, function(model) {data.frame(temp = temps[[model]],
                                                      flux = fluxes[[model]])})

  # name each dataframe in list
  names(CMIP5) <- models

  # return list of dataframes
  return(CMIP5)

}

# names of models
models <- c("BCC", "BNU", "CCCMA", "CNRM",
            "CSIRO", "GFDL", "GISS", "IAP",
            "INM", "IPSL", "MIROC", "MOHC",
            "MPIM", "MRI", "NCAR", "NCC")

# read in data
temps <- ReadTemps(models)
fluxes <- ReadFluxes(models)

# build CMIP5 list of dataframes
CMIP5 <- BuildCMIP5(models, temps, fluxes)





