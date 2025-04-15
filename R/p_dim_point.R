#' Conformalized Point Hyperrectangles
#'
#' Creates conformalized hyperrectanles using univariate point (e.g. conditional mean) regression
#'
#' @param Y the n by p matrix of the response
#' @param X the n by d matrix of the covariates
#' @param alpha numeric; the miscoverage rate
#' @param train_prop numeric; the proportion of data used to train the conditional point models
#' @param cal1_prop numeric; the proportion of data used to form the initial (unconformalized) hyper rectangles
#' @param x_test the matrix containing the covariates that are used to create the prediction regions
#' @return A list with the following components: \item{lower}{matrix of m by p with marginal
#' lower bounds} \item{upper}{matrix of m by p with marginal upper bounds}
#' \item{adj}{vector of length p, the total adjustment for the point predictions }
#' @examples
#' set.seed(42)
#' n <- 1000; d <- 1; p <- 2
#' X <- matrix(rnorm(n * d), n, d)
#' Y <- matrix(7 * X ^ 2, nrow = n, ncol = p, byrow = FALSE) + matrix(rnorm(n * p), ncol = p)
#' x_test <- matrix(rnorm(10 * d), nrow = 10)  # 10 test points
#' y_test <- matrix(7 * x_test ^ 2, nrow = 10, ncol = p) + rnorm(10 * p)
#' res <- p_dim_point_conformal(Y, X, alpha = 0.1, train_prop = 0.5, cal1_prop = 0.2, x_test = x_test)
#' cov <- mean(rowMeans(res$lower < y_test & res$upper > y_test) == 1) ## out of sample coverage

#' # Check shape
#' dim(res$lower)  # should be 10 x 2
#' c(res$lower[1, ], res$upper[1, ])  # marginal intervals for the initial test point

#' @export p_dim_point_conformal
#' @import ranger
#' @importFrom stats predict

p_dim_point_conformal <- function(Y, X, alpha = 0.1, train_prop = 0.5, cal1_prop = 0.25, x_test) {
    stopifnot(nrow(Y) == nrow(X), is.numeric(alpha), alpha > 0, alpha < 1)

    if (is.vector(x_test)) {
        x_test <- matrix(x_test, nrow = 1)
    }

    n <- nrow(X)
    p <- ncol(Y)
    m <- nrow(x_test)

    # Indices for splitting the data
    idx <- sample(1:n)
    n_train <- floor(train_prop * n)
    n_cal1 <- floor(cal1_prop * n)
    n_cal2 <- n - n_train - n_cal1

    idx_train <- idx[1:n_train]
    idx_cal1 <- idx[(n_train + 1):(n_train + n_cal1)]
    idx_cal2 <- idx[(n_train + n_cal1 + 1):n]

    # Data splits
    X_train <- X[idx_train, , drop = FALSE]
    Y_train <- Y[idx_train, , drop = FALSE]
    X_cal1 <- X[idx_cal1, , drop = FALSE]
    Y_cal1 <- Y[idx_cal1, , drop = FALSE]
    X_cal2 <- X[idx_cal2, , drop = FALSE]
    Y_cal2 <- Y[idx_cal2, , drop = FALSE]

    # Step 2: Train point estimate models using ranger
    models_fhat <- lapply(1:p, function(j) {
        ranger::ranger(y ~ ., data = data.frame(y = Y_train[, j], x = X_train))
    })

    # Step 3: Compute absolute residuals |Y - fhat(X)| on cal1 set
    preds_cal1 <- sapply(1:p, function(j) {
        predict(models_fhat[[j]], data = data.frame(x = X_cal1))$predictions
    })

    scores <- abs(Y_cal1 - preds_cal1)  # columns of scores for each dimension

    # Step 4: Get conformal intervals on cal2 points using scores from cal1
    quantile_scores <- apply(scores, 2, function(v) {
        sort(v)[ceiling((1 - alpha) * (1 + length(v)))]
    })  # length-p vector

    # Step 5: Get point predictions on cal2
    preds_cal2 <- sapply(1:p, function(j) {
        predict(models_fhat[[j]], data = data.frame(x = X_cal2))$predictions
    })

    # Step 6: Side lengths (upper - lower) for cal2
    side_lengths_cal2 <- 2 * quantile_scores

    # Step 7: Compute ratios relative to dimension 1
    ratios <- side_lengths_cal2[1] / side_lengths_cal2[1:p]

    # Step 8: Compute E_{k,j} on cal2
    E_kj <- matrix(NA, nrow = length(preds_cal2[, 1]), ncol = p)
    for (j in 1:p) {
        q_lo_cal2 <- preds_cal2[, j] - quantile_scores[j]
        q_hi_cal2 <- preds_cal2[, j] + quantile_scores[j]
        E_kj[, j] <- pmax(q_lo_cal2 - Y_cal2[, j], Y_cal2[, j] - q_hi_cal2)
    }

    # Step 9: Convert to first dimension using ratios
    A_kj <- E_kj * ratios

    # Step 10: W_k = max_j A_{k,j}
    W_k <- apply(A_kj, 1, max)

    # Step 11: Adjustment in first dimension
    adj1 <- sort(W_k)[ceiling((1 - alpha) * (1 + length(W_k)))]

    # Step 12: Final adjustments for each test point dimension using average ratios
    Adj <- adj1 * 1 / ratios ##since the ratios are constant, we don't need to loop over m

    # Step 13: Get point predictions on x_test
    preds_test <- sapply(1:p, function(j) {
        predict(models_fhat[[j]], data = data.frame(x = x_test))$predictions
    })
    preds_test <- matrix(preds_test, nrow = m, ncol = p)

    Adj <- Adj + quantile_scores

    lower <- preds_test - Adj
    upper <- preds_test + Adj

    list(
        lower = lower,        # m × p matrix
        upper = upper,        # m × p matrix
        adj = Adj             # adjustment matrix
    )
}
