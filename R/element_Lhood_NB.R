#' # Likelihood of observing outbreaks of size y given parameters
#'
#'
#'  The function internally relies on pre-computing:
#'  1) g: the probability of observing size y given z,rho; with z the true number of
#'  cases (reported and un-reported cases) and rho the reporting probability
#'
#' 2) g0: the probability of not observing an outbreak (i.e. y=0 or no cases) given z and rho.
#'
#' g and g0 are obtained using the function proba_observation.
#' 
#' adapted from element_Lhood_poisson but accounting for a NB offspring distribution with overdispersion 'over'
#'
#' @author Pierre Nouvellet (\email{p.nouvellet@imperial.ac.uk})
#'
#' @export
#'
#' @param R is the reproduction number, i.e. the average number of secondary cases due to a single case.
#' This can be any positive number.
#' 
#' @param over is the overdispersion in the offspring distribution.
#'
#' @param z is the true potential outbreak sizes. It is precomputed by proba_observation.
#'
#' @param g: the probability of observing size y given z,rho; with z the true number of
#'  cases (reported and un-reported cases) and rho the reporting probability. It is precomputed by proba_observation.
#'
#' @param g0: the probability of not observing an outbreak (i.e. y=0 or no cases) given z and rho.
#' It is precomputed by proba_observation.
#'
#'
#' @return
#'  The function returns the likelihood of the observations for a given {R,rho,overdispersion}.
#'
#'
#' @examples
#'
#' y <- element_Lhood_NB(R = .5, over = 1, z = x$possible_size, g = x$p_y_z, g0 = x$p_0_z)
#' y
#'


element_Lhood_NB<-function(R,over,z,g,g0){

  R_eff <- R_eff_NB(R = R, over = over)
  # f <- f_z_poiss(z,1,R_eff$R_effective)*R_eff$P_extinction
  f <- f_z_NB(z = z, s = 1, R = R_eff$R_effective, over = over)*R_eff$P_extinction
  # get the likelihood
  Likelihood <- (sum(log((g %*% f)))-nrow(g)*log(1-g0 %*% f))
  # correction for threshold?
  # L2 <- -log( sapply(yobs,
  #                    function(x) correct.tail(x,R,k,p,threshold))  )
  return( Likelihood )
}
