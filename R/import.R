#' # Estimate the number of clusters, including those unobserved
#'
#'
#' @author Pierre Nouvellet (\email{p.nouvellet@imperial.ac.uk})
#'
#' @export
#'
#' @param y_obs is the reproduction number, i.e. the average number of secondary cases due to a single case.
#' This can be any positive number.
#'
#' @param rho is the true potential outbreak sizes. It is precomputed by proba_observation.
#'
#' @param profile: the probability of observing size y given z,rho; with z the true number of
#'  cases (reported and un-reported cases) and rho the reporting probability. It is precomputed by proba_observation.
#'
#' @param threshold_z: maximum size of unobserved cluster.
#'
#' @param threshold_import: maximum number of importations evaluated.
#'
#' @param CI optional parameter for the level of the confidence interval.Default
#' is 0.95 for 95% confidence interval
#'
#' @param over is the overdispersion in the offspring distribution.
#'
#' @return
#'  The function returns the maximum likelihood of the total number of clusters
#'  (observed and unobserved), associated confidence intervals and maximum likelohood
#'
#'
#'
#' @examples
#'
#' x <- import(y_obs = y_obs_test, rho = .5,
#'                 profile = prof_test,
#'                 threshold_z = 1e3,
#'                 threshold_import = 1e3)
#' x
#'
import <- function(y_obs,rho,profile,threshold_z,threshold_import, CI = .95, k = NULL){

  z <- matrix(seq(1,threshold_z), nrow = threshold_z, ncol = 1)
  g0 <- dbinom(0,matrix(z,nrow = 1,ncol = threshold_z,byrow=TRUE),rho)

  profile$import <- 0:threshold_import
  profile$Lk_import <- rep(0,threshold_import+1)

  A <- matrix(NA,nrow = threshold_import+1, ncol = length(profile$theta))
  if (is.null(k)){
    for (i in 1:length(profile$theta)){
      R_eff <- R_eff_poisson(R = profile$theta[i])
      f <- f_z_poiss(z = z, s = 1, R = R_eff$R_effective)*R_eff$P_extinction
      p_obs<- (1-g0 %*% f)


      A[,i] <- dnbinom(profile$import, length(y_obs), p_obs, log = TRUE) +
        profile$Likelihood[i]

      # temp <- dnbinom(profile$import, length(y_obs), p_obs, log = TRUE) +
      #   profile$Likelihood[i]
      # profile$Lk_import <- profile$Lk_import + exp(temp)
    }
  } else if (k>0){
    for (i in 1:length(profile$theta)){
      R_eff <- R_eff_NB(R = profile$theta[i], k = k)
      f <- f_z_NB(z = z,s = 1, R = R_eff$R_effective, k = k)*R_eff$P_extinction
      p_obs<- (1-g0 %*% f)


      A[,i] <- dnbinom(profile$import, length(y_obs), p_obs, log = TRUE) +
        profile$Likelihood[i]

      # temp <- dnbinom(profile$import, length(y_obs), p_obs, log = TRUE) +
      #   profile$Likelihood[i]
      # profile$Lk_import <- profile$Lk_import + exp(temp)
    }
  }
  b <- max(c(A))
  profile$Lk_import <- b + log(rowSums( exp(A - b) ) )

  max_likelihood <- theta_max_likelihood(profile$import,profile$Lk_import,CI)

  return ( list (theta_max_likelihood = max_likelihood$theta,
                 max_likelihood = max_likelihood$likelihood,
                 lower_theta =  max_likelihood$lower_theta,
                 upper_theta = max_likelihood$upper_theta) )
}
