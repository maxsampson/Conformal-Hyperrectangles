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
- cal1_prop (only for point regression) is the proportion of the data used to form in the initial hyperrectangles when using a point estimator
- x_test is an m by d matrix of the covariates we are interested in creating hyperrectangular prediction regions for

A simple example of both functions is given below

```R
set.seed(1)
n <- 10000; d <- 1; p <- 3
X <- matrix(rnorm(n * d), n, d)
Y <- matrix(7 * X ^ 2, nrow = n, ncol = p, byrow = FALSE) + matrix(rnorm(n * p), ncol = p)
x_test <- matrix(rnorm(10 * d), nrow = 10)  # 10 test points
y_test <- matrix(7 * x_test ^ 2, nrow = 10, ncol = p) + rnorm(10 * p)
res <- p_dim_quantile_conformal(Y, X, alpha = 0.1, train_prop = 0.7, x_test = x_test)
mean(rowMeans(res$lower < y_test & res$upper > y_test)) ## out of sample coverage with the quantile approach

res_point <- p_dim_point_conformal(Y, X, alpha = 0.1, train_prop = 0.5, cal1_prop = 0.2, x_test = x_test)
mean((rowMeans(res_point$lower < y_test & res_point$upper > y_test)) ## out of sample coverage with the point approach

```

The models used are quantile random forests and traditional random forests for the quantile and point approach, respectively. These can be changed by the user by editing the functions found in the `R` folder. 

## License

This package is free and open source software, licensed under GPL 3.
