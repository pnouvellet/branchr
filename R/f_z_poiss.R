#' Compute the probability of observing an outbreak of size z given an R<1 and assuming a Poisson offspring distribution
#'
#'
#' @author Pierre Nouvellet (\email{p.nouvellet@imperial.ac.uk})
#'
#' @export
#'
#' @param R is the reproduction number, i.e. the average number of secondary cases due to a single case.
#' This can be any positive number.
#'
#' @return
#'  The function returns a list including:
#' \itemize{
#'
#'  \item the probability of observing an outbreak of size z given R.
#'  
#'  }
#' 
#'
#'
#' @examples
#'
#' x <- f_z_poiss(z = 1:1e2, s = 1, R = .5)
#' x
#'
#
f_z_poiss <- function(z,s,R){
  logf<-(z-s-1)*log( s*z)+(z-s)*log(R)+ (-z*R) - lgamma(z-s +1)
  f<-exp(logf)
  return(proba_z = f)
}
