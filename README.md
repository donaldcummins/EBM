# EBM

EBM is an R package for maximum likelihood estimation of k-box stochastic energy balance models.

## Getting started

EBM can be installed directly from GitHub using the [devtools](https://cran.r-project.org/package=devtools) package for R.

### Installing devtools

To install devtools from [CRAN](https://cran.r-project.org/), open an R session and run

```
install.packages("devtools")
```

If devtools requires compilation on your system this step may take a few minutes.

### Installing EBM

To install EBM using devtools, run

```
devtools::install_github("donaldcummins/EBM")
```

Once installed, EBM can be loaded into an R session using

```
library(EBM)
```

## Licence

EBM is licenced under the GNU General Public License v3.0 - see the [LICENCE](LICENCE) file for details.

## Acknowledgements

EBM depends on the following R packages:

* [expm](https://cran.r-project.org/package=expm)
* [FKF](https://cran.r-project.org/package=FKF)
* [nloptr](https://cran.r-project.org/package=nloptr)
* [numDeriv](https://cran.r-project.org/package=numDeriv)

We thank the respective authors for making their software available.





