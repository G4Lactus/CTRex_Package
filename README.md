CT-Rex Selector
================

## CT-Rex Selector

The CT-Rex selector performs complex-valued terminated-random
experiments (CT-Rex) using the CT-LARS algorithm and fuses the selected
sets of all random experiments to obtain a final set of selected
variables. It provably controls the false discovery rate (FDR) at a
user-defined target level while maximizing the number of selected
variables. Thereby it achieves a high true positive rate (TPR) (i.e.,
power). The CT-Rex selector can be applied in various fields, such as in
engineering, communications, physics, and all disciplines being
presented with the problem of high-dimensional FDR control in the
complex number domain.

Note: A Python and C++ Implementation is coming soon.

# Installation

Before installing the `CT-Rex` package, you need to install the required
`CT-LARS` package from
[Github](https://github.com/G4Lactus/CTLARS_Package) (developer
version). Afterwards you can install the `CT-Rex` package
[Github](https://github.com/G4Lactus/CTrex_Package).

``` r
# install.packages("devtools")
# devtools::install_github("https://github.com/G4Lactus/CTLARS_Package")
# devtools::install_github("https://github.com/G4Lactus/CTrex_Package")
```

# Quick Start

In the following we present a short demo how to use it.

``` r
set.seed(42)
n <- 100
p <- 250

beta_cardinality <- 5
beta_support <- sort(sample(p, beta_cardinality))
beta <- rep(0, times = p)
beta[beta_support] <- exp(1i * stats::runif(beta_cardinality, 0, 2*pi))

X <- matrix(data = complex(
                    real = stats::rnorm(n * p),
                    imaginary = stats::rnorm(n * p)),
            nrow = n,
            ncol = p)

y <- X[, beta_support] %*% beta[beta_support]
noise <- complex(real = stats::rnorm(n, sd = 1/sqrt(2)),
                 imaginary = stats::rnorm(n, sd = 1/sqrt(2)))
                 
y <- y + noise
                 
ctrex_obj <- CTrex::ctrex(X = X, y = y, tFDR = 0.1)

cat("Selected variables:\n", which(ctrex_obj$selected_var > 0))
```

# Links

If you interested in FDR controlled variable selection in the real
number space have a look at:

T-Rex paper: <https://doi.org/10.48550/arXiv.2110.06048>

TRexSelector package (stable version):
[CRAN-TRexSelector](https://CRAN.R-project.org/package=TRexSelector).

TRexSelector package (developer version):
[GitHub-TRexSelector](https://github.com/jasinmachkour/TRexSelector).

README file:
[GitHub-readme](https://htmlpreview.github.io/?https://github.com/jasinmachkour/TRexSelector/blob/main/README.html).

Vignette:
[GitHub-vignette](https://htmlpreview.github.io/?https://github.com/jasinmachkour/TRexSelector/blob/main/vignettes/TRexSelector_usage_and_simulations.html).

tlars package: [CRAN-tlars](https://CRAN.R-project.org/package=tlars)
and [GitHub-tlars](https://github.com/jasinmachkour/tlars).
