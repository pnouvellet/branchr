#' Compute the probability of observing an outbreak of size z given an R<1 and assuming a NB offspring distribution
#'
#' adapted from f_z_poisson but accounting for a NB offspring distribution with overdispersion 'over'
#'
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
#' @return
#'  The function returns a list including:
#' \itemize{
#'
#' \item the probability of observing an outbreak of size z given R/overdispersion
#' 
#'}
#'
#' 
#'
#'
#' @examples
#'
#' x <- f_z_NB(z = 1:1e2, s = 1, R = .5, over = 1)
#' 
#' 
#
f_z_NB <- function(z,s,R,over){

  logf<-lgamma(z*over+z-1)-lgamma(z*over)-lgamma(z+1) + (z-1)*log(R/over) - (z*over+z-1)*log(R/over+1)
  f<-exp(logf)
  
  # logf<-(z-s-1)*log( s*z)+(z-s)*log(R)+ (-z*R) - lgamma(z-s +1)
  # f<-exp(logf)

  # plot(Re(Psize))
  # sum(Psize,na.rm=TRUE)

  return(proba_z=f)



}
