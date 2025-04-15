#' Conformalized Quantile Hyperrectangles
#'
#' Creates conformalized hyperrectanles using univariate quantile regression
#'
#' @param Y the n by p matrix of the response
#' @param X the n by d matrix of the covariates
#' @param alpha numeric; the miscoverage rate
#' @param train_prop numeric; the proportion of data used to train the quantile models
#' @param x_test the matrix containing the covariates that are used to create the prediction regions
#' @return A list with the following components: \item{lower}{matrix of m by p with marginal
#' lower bounds} \item{upper}{matrix of m by p with marginal upper bounds}
#' \item{adj}{matrix of m by p with conformal adjustments used for marginal quantile regression}
#' @examples
#' set.seed(42)
#' n <- 1000; d <- 1; p <- 2
#' X <- matrix(rnorm(n * d), n, d)
#' Y <- matrix(7 * X ^ 2, nrow = n, ncol = p, byrow = FALSE) + matrix(rnorm(n * p), ncol = p)
#' x_test <- matrix(rnorm(10 * d), nrow = 10)  # 10 test points
#' y_test <- matrix(7 * x_test ^ 2, nrow = 10, ncol = p) + rnorm(10 * p)
#' res <- p_dim_quantile_conformal(Y, X, alpha = 0.1, train_prop = 0.7, x_test = x_test)
#' cov <- mean(rowMeans(res$lower < y_test & res$upper > y_test) == 1) ## out of sample coverage

#' # Check shape
#' dim(res$lower)  # should be 10 x 2
#' c(res$lower[1, ], res$upper[1, ])  # marginal intervals for the initial test point
#' @export p_dim_quantile_conformal
#' @import quantregRanger
#' @importFrom stats predict

p_dim_quantile_conformal <- function(Y, X, alpha = 0.1, train_prop = 0.7, x_test) {
    stopifnot(nrow(Y) == nrow(X), is.numeric(alpha), alpha > 0, alpha < 1)

    if (is.vector(x_test)) {
        x_test <- matrix(x_test, nrow = 1)
    }

    n <- nrow(X)
    p <- ncol(Y)
    m <- nrow(x_test)  # number of test points

    # Split into training and calibration sets
    idx <- sample(1:n)
    n_train <- floor(train_prop * n)
    idx_train <- idx[1:n_train]
    idx_cal <- idx[(n_train + 1):n]

    X_train <- X[idx_train, , drop = FALSE]
    Y_train <- Y[idx_train, , drop = FALSE]
    X_cal <- X[idx_cal, , drop = FALSE]
    Y_cal <- Y[idx_cal, , drop = FALSE]

    # Fit quantile regression forest for each dimension of Y
    models <- lapply(1:p, function(j) {
        quantregRanger::quantregRanger(y ~ ., data = data.frame(y = Y_train[, j], x = X_train))
    })

    # Predict lower and upper quantiles for calibration data
    q_lo <- sapply(1:p, function(j) predict(models[[j]], data.frame(x = X_cal), quantiles = alpha / 2))
    q_hi <- sapply(1:p, function(j) predict(models[[j]], data.frame(x = X_cal), quantiles = 1 - alpha / 2))
    side_lengths <- q_hi - q_lo

    # E_{k,j}
    E_kj <- pmax(q_lo - Y_cal, Y_cal - q_hi)

    # W_{k,j}
    W_kj <- E_kj
    for (k in 1:nrow(E_kj)) {
        for (j in 2:p) {
            W_kj[k, j] <- E_kj[k, j] * (side_lengths[k, 1] / side_lengths[k, j])
        }
    }

    # W_k = max_j W_kj
    W_k <- apply(W_kj, 1, max)

    # Adjustment for dimension 1
    adj1 <- sort(W_k)[ceiling((1 - alpha) * (1 + length(W_k)))]

    # Predict for all test points
    q_lo_test <- sapply(1:p, function(j) predict(models[[j]], data.frame(x = x_test), quantiles = alpha / 2))
    q_hi_test <- sapply(1:p, function(j) predict(models[[j]], data.frame(x = x_test), quantiles = 1 - alpha / 2))

    # Ensure shape: m × p
    q_lo_test <- matrix(q_lo_test, nrow = m, ncol = p)
    q_hi_test <- matrix(q_hi_test, nrow = m, ncol = p)
    side_lengths_test <- q_hi_test - q_lo_test

    # Final adjustments for all test points
    Adj <- matrix(NA, nrow = m, ncol = p)
    Adj[, 1] <- adj1
    for (j in 2:p) {
        Adj[, j] <- adj1 * (side_lengths_test[, j] / side_lengths_test[, 1])
    }

    # Final prediction intervals
    lower <- q_lo_test - Adj
    upper <- q_hi_test + Adj

    list(
        lower = lower,        # m × p matrix
        upper = upper,        # m × p matrix
        adj = Adj
    )
}
