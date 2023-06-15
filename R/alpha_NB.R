#' Obtain alpha for a Poisson offspring distribution
#'
#' Same as alpha_Poisson but applicable for NB offspring distribution with overdispersion 'over'
#'
#'
#' @author Pierre Nouvellet (\email{p.nouvellet@imperial.ac.uk})
#'
#' @export
#'
#' @param R is the reproduction number, i.e. the average number of secondary cases due to a single case.
#' This can be any positive number. if R is a vector, then the length of R must be 'n' (see below).
#'
#' @param over is the overdispersion in the offspring distribution.
#'
#'

#' @return
#'  The function returns alpha assuming a NB offspring distribution.
#'
#'
#' @examples
#'
#' x <- alpha_NB(R = 1.5, over = 1)
#' x
#'

alpha_NB <- function(R,over){
  
  alpha <- rep(NA,length(R))
  f_sub <- which(R<=1)
  f_super <- which(R>1)
  for(i in 1:length(f_sub)){
    alpha[f_sub[i]] <- 0
  }
  for(i in 1:length(f_super)){
    alpha0 <- R+lamW::lambertW0(-R*exp(-R))
    temp <- optimise(f = f1nb, interval = c(1e-10,1e2), R = R, over = over)
    alpha[f_super[i]] <- unlist(temp[1])
  }

  return(alpha=alpha)
}

f1nb <- function(alpha,R,over){
  p <- R/ (R+over)
  temp <- exp(alpha) * ((1-p)/(1-p*exp(-alpha)))^over -1
  res <- ( temp )^2
  return(res)
}

