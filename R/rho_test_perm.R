#' @title rho_test_perm
#'
#' @description Computes the estimate of the A-R statistic and independence inference using a permutation test.
#'
#' @param x a numeric vector.
#' @param y a numeric vector with compatible length to \code{x}.
#' @param R the number of permutations.
#' @param verbose print tracing information. The Default value is \code{FALSE}.
#' @param return_dist should the simulated distribution of the A-R statistic be returned. The Default value is \code{TRUE}.
#'
#' @return A list object.
#' \itemize{
#'    \item{\code{rho_minus_one_mean}}{ estimate of the expected value of the A-R statistic.}
#'    \item{\code{rho_minus_one_sample}}{ if \code{return_dist = TRUE} returns simulated values of the A-R statistic.}
#'    \item{\code{p_value}}{ the test p-value.}
#'  }
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' x <- rnorm(n)
#' y <- rnorm(n)
#' res <- rho_test_perm(x, y)
#' res$rho_minus_one_mean
#' res$p_value
#' hist(res$rho_minus_one_sample, probability = TRUE)
#' @export

rho_test_perm <- function(x, y = NULL, R = 1000, verbose = FALSE, return_dist = TRUE){

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

  rho_minus_ones <- NA

  if(return_dist){
    rho_minus_ones <- rho(x, y, return_dist = TRUE)
  }

  rho_minus_one_mean <- rho(x, y)

  rho_minus_one_perm <- rep(0, R)

  for(i in 1:R){

    rho_minus_one_perm[i] <- rho(sample(x), sample(y))

    if(verbose) cat("\r", round(100*i/R, 2), "%")

  }

  p_value <- mean(rho_minus_one_perm > rho_minus_one_mean)

  return(list(rho_minus_one_mean = rho_minus_one_mean,
              rho_minus_one_sample = rho_minus_ones,
              p_value = p_value)
         )

}
