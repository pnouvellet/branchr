#' Obtain the maximum likelihood from a vector of parameter and its likelihood profile
#'
#' Also gives confidence interval
#'
#'
#' @author Pierre Nouvellet (\email{p.nouvellet@imperial.ac.uk})
#'
#' @export
#'
#' @param theta is a vector containing parameter's values for which the
#' likelihood has been evaluated
#'
#' @param likelihood at the parameter's value (log-likelihood)
#'
#' @param CI is the level to estimate the confidence interval, i.e. the default,
#' CI = 0.95, return the 95percent confidence interval.
#'
#'
#' @return
#'  The function returns a list including:
#' \itemize{
#'
#' \item theta_max_likelihood: the maximum likelihood estimated with the range
#' of theta values evaluated.
#'
#' \item max_likelihood: the log of the maximum likelihood.
#'
#' \item lower_theta,upper_theta: lower and upper bound of the 'CI'percent
#' confidence interval.
#'
#'
#' }
#'
#'
#'
#'
theta_max_likelihood <- function(theta,likelihood,CI = 0.95){

  max_likelihood <- list( likelihood = max(likelihood) )
  max_likelihood$index <- which(likelihood %in% max_likelihood$likelihood)
  max_likelihood$theta <- theta[max_likelihood$index]

  limit <- qchisq(CI, df=1)/2

  CI_profile<-(likelihood-max_likelihood$likelihood+limit)^2
  lower_CI <- theta[ which( CI_profile[1:max_likelihood$index]
                          %in% min(CI_profile[1:max_likelihood$index])) ]

  upper_CI <- theta[ max_likelihood$index - 1 +
                     which( CI_profile[max_likelihood$index:length(CI_profile)]
                            %in% min(CI_profile[max_likelihood$index:length(CI_profile)])) ]

  return ( list (theta_max_likelihood = max_likelihood$theta,
                      max_likelihood = max_likelihood$likelihood,
                      lower_theta = lower_CI, upper_theta = upper_CI) )
}
