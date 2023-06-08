#' @title rho_rank
#'
#' @description Estimator of the A-R statistic using ranks. Faster than \code{rho()} when variables are large. The estimator is the same as the one available in the function \code{rho()} when there are no ties in the data.
#'
#' @param x a numeric vector.
#' @param y a numeric vector with compatible length to \code{x}.
#' @param return_dist return the simulated distribution of the A-R statistic. The Default value is \code{FALSE}.
#' @param m the number of times the A-R algorithm is run.
#' @param ad_hoc_check check that the estimate of the statistic is not larger than one. The Default value is \code{TRUE}.
#' @param scale_monotone return the standardized coefficient. The Default value is \code{TRUE}.
#'
#' @return Estimate of the expected value or simulated values of the A-R statistic
#'
#' @examples
#' set.seed(1)
#' n = 1000
#' x = rnorm(n)
#' y = rnorm(n)
#' system.time(r <- rho(x, y))
#' system.time(r_rank <- rho_rank(x, y))
#' r; r_rank
#' # When there are ties,
#' x = sample(x, replace = TRUE)
#' rho(x, y); rho_rank(x, y)
#' @export
#' @importFrom stats runif

rho_rank = function(x, y = NULL, return_dist = FALSE, m = 10^4, ad_hoc_check = TRUE, scale_monotone = TRUE){

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

  n = length(x)

  d_cond = matrix(0, n, 4)

  d_marg = matrix(0, n, 4)

  y_rank_rev = rank(-y, ties.method = "max")
  y_rank = rank(y, ties.method = "max")

  x_rank_rev = rank(-x, ties.method = "max")
  x_rank = rank(x, ties.method = "max")

  for(i in 1:n){

    d_cond[i, 3] = sum(y[x <= x[i]] <= y[i])

  }

  d_cond[, 1] = x_rank + 1 - d_cond[, 3]
  d_cond[, 2] = y_rank_rev + 1 - d_cond[, 1]
  d_cond[, 4] = y_rank + 1 - d_cond[, 3]

  d_cond[, 1] = d_cond[, 1]/x_rank
  d_cond[, 2] = d_cond[, 2]/x_rank_rev
  d_cond[, 3] = d_cond[, 3]/x_rank
  d_cond[, 4] = d_cond[, 4]/x_rank_rev

  d_marg[, c(1, 2)] = y_rank_rev/n
  d_marg[, c(3, 4)] = y_rank/n

  a = d_cond/d_marg

  a = as.vector(a)

  i = 2:(n-1)

  d = 1/2 + 1/n + (1/2)*sum(1/(i*(n + 1 - i)))

  d = scale_monotone*d

  if(return_dist){

    p = rep(0, m)

    for(i in 1:m){

      P = (1 - mean(a > runif(length(a)), na.rm = TRUE))/(1 - d) # scaled proportion

      p[i] = P - ad_hoc_check*((P > 1)*(P - 1))

    }

  }

  if(!return_dist){

    b = a < 1

    b = a*b

    b[b == 0] = NA

    p = (1 - (sum(a >= 1) + sum(b, na.rm = TRUE))/(4*n))/(1 - d) # scaled proportion

    p = p - ad_hoc_check*((p > 1)*(p - 1))

  }

  return(p)

}

