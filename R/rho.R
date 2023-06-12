#' @title rho
#'
#' @description Estimator of the A-R statistic.
#'
#' @param x a numeric vector.
#' @param y a numeric vector with compatible length to \code{x}.
#' @param return_dist return the simulated distribution of the A-R statistic. Default \code{FALSE}.
#' @param m the number of times the A-R algorithm is run.
#' @param ad_hoc_check check that the estimate of the statistic is not larger than one. Default \code{TRUE}.
#' @param scale_monotone return the standardized coefficient. Default is \code{TRUE}.
#'
#' @return Estimate of the expected value or simulated values of the A-R statistic.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' x <- rnorm(n)
#' y <- rnorm(n)
#' rho(x, y)
#' d <- rho(x, y, return_dist = TRUE)
#' hist(d, probability = TRUE)
#'
#' y <- (x - mean(x))^2 + rnorm(n, sd = 0.1)
#' rho(x, y)
#' d <- rho(x, y, return_dist = TRUE)
#' hist(d, probability = TRUE)
#' @export
#' @importFrom stats runif

rho <- function(x, y = NULL, return_dist = FALSE, m = 10^4, ad_hoc_check = TRUE, scale_monotone = TRUE){

  if(!(is.vector(x) || is.vector(y)))
    stop("'x' and 'y' must be numeric vectors")
  if(!is.vector(x) || is.null(y))
    stop("supply both 'x' and 'y' ")
  if(!(is.numeric(x) || is.logical(x)))
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  if(!(is.numeric(y) || is.logical(y)))
    stop("'y' must be numeric")
  stopifnot(is.atomic(y))
  if(length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  if (length(x) == 0L || length(y) == 0L)
    stop("both 'x' and 'y' must be non-empty")

  n <- length(x)

  d_cond <- matrix(0, n, 4)

  d_marg <- matrix(0, n, 4)

  for(i in 1:n){

    d_cond[i, 1] <- mean(y[x <= x[i]] >= y[i])
    d_cond[i, 2] <- mean(y[x >= x[i]] >= y[i])
    d_cond[i, 3] <- mean(y[x <= x[i]] <= y[i])
    d_cond[i, 4] <- mean(y[x >= x[i]] <= y[i])

  }

  y_rank_rev <- rank(-y, ties.method = "max")
  y_rank <- rank(y, ties.method = "max")

  d_marg[, c(1, 2)] <- y_rank_rev/n
  d_marg[, c(3, 4)] <- y_rank/n

  a <- d_cond/d_marg

  a <- as.vector(a)

  i <- 2:(n-1)

  d <- 1/2 + 1/n + (1/2)*sum(1/(i*(n + 1 - i)))

  d <- scale_monotone*d

  if(return_dist){

    p <- rep(0, m)

    for(i in 1:m){

      P <- (1 - mean(a > runif(length(a)), na.rm = TRUE))/(1 - d)

      p[i] <- P - ad_hoc_check*((P > 1)*(P - 1))

    }

  }

  if(!return_dist){

    b <- a < 1

    b <- a*b

    b[b == 0] <- NA

    p <- (1 - (sum(a >= 1) + sum(b, na.rm = TRUE))/(4*n))/(1 - d)

    p <- p - ad_hoc_check*((p > 1)*(p - 1))

  }

  return(p)

}
