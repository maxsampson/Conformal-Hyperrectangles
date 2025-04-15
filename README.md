# Conformal-Hyperrectangles
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

*Conformal-Hyperrectangles* is an R package that can be used to compute conformalized hyperrectangles for multi-target responses.

Sampson, M. and Chan, K.S. (2024). [Conformal Multi-Target Hyperrectangles](https://onlinelibrary.wiley.com/doi/full/10.1002/sam.11710)

To install the package, run

```R
# install.packages("devtools")
devtools::install_github("maxsampson/Conformal-Hyperrectangles")
```

The functions p_dim_point_conformal and p_dim_quantile_conformal require 6 and 5 inputs, Y, X, alpha, train_prop, cal1_prop (only for point regression), and x_test.

- Y is an n by p matrix of the observed responses
- X is an n by d matrix of the observed covariates
- alpha is a numeric; the desired miscoverage level
- train_prop is the proportion of data used to train the initial models
- cal1_prop is the proportion of the data used to form in the initial hyperrectangles when using a point estimator
- x_test is an m by d matrix of the covariates we are interested in creating hyperrectangular prediction regions for

A simple example to be written later

```R
set.seed(1)
```

## License

This package is free and open source software, licensed under GPL 3.
