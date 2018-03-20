
##########################

#' The Generalized logistic distribution
#'
#' Density, distribution function, quantile function and random generation a generalized logistic distribution.
#' @param x,q vector of quantiles.
#' @param k vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required
#' @param a,b,p  parameters \eqn{\ge 0}, with restrictions.*
#' @param mu mu parameter
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @keywords genlogis
#' 
#' @export
#' @examples 
#' pgenlog(0.5) 
#' curve(dgenlog(x), xlim = c(-3,3)) 
#' 
#' rgenlog(100) 
#' 
#' qgenlog(0.95)
#' 
#' 
#' @name distrib
#'  
#' @return 
#' \code{dgenlog} gives the density, \code{pgenlog} gives the distribution function, 
#' \code{qgenlog} gives the quantile function, and \code{rgenlog} generates random deviates.\cr
#'  
#'  The length of the result is determined by \code{n} for \code{rgenlog}, and is the maximum of the lengths 
#'  of the numerical arguments for the other functions.
#' 
#' @details 
#' 
#' The used distribution for this package is given by: 
#' \deqn{f(x) = ((a + b*(1+p)*(|x-mu|^p))*exp(-(x-mu)*(a+b*(|x-mu|^p)))) / ((exp(-(x-mu)*(a + b* (|x-mu|^p)))+1)^2)}
#' 
#' The default values for \code{a, b, p and mu} produces a function with mean 0 and variance close to 1.\cr 
#' 
#' *Restrictions:\cr 
#' 
#' If \code{p} equals to 0, \code{b} or \code{a} must be 0 otherwise there is identifiability problem.\cr 
#' 
#' The distribution is not defined for \code{a} and \code{b} equal to 0 simultaneously.\cr 
#' 
#' @references 
#' Rathie, P. N. and Swamee, P. K (2006) \emph{On a new invertible generalized logistic distribution
#' approximation to normal distribution}, Technical Research Report in Statistics, 07/2006,
#' Dept. of Statistics, Univ. of Brasilia, Brasilia, Brazil.

#' @rdname distrib
#' @export

pgenlog <- function(q, a = sqrt(2/pi), b = 0.5, p = 2, mu = 0, lower.tail = TRUE){

  if(!missing(a)){
    if(a < 0){
      stop('The argument "a" must be positive.')
    }
  }
  if(!missing(b)){
    if(b < 0){
      
      stop('The argument "b" must be positive.')
    }
  }
  if(!missing(p)){
    if(p < 0){
      stop('The argument "p" must be positive.')
    }
  }
  
  if(p == 0 && b > 0 && a > 0){
    stop('If "p" equals to 0, "b" or "a" must be 0 otherwise there is identifiability problem.')
  }  
  if(b == 0 && a == 0){
    stop('The distribution is not defined for "a" and "b" equal to 0 simultaneously.')
  }
  
  z <- (exp(-(q-mu)*(a+b*(abs(q-mu)^p)))+1)^(-1)
  z <- ifelse(lower.tail == TRUE, z, 1-z)
  
  return(z)
}


#' @rdname distrib
#' @export

dgenlog <- function(x, a = sqrt(2/pi), b = 0.5, p = 2, mu = 0){
  
  if(!missing(a)){
    if(a < 0){
      stop('The argument "a" must be positive.')
    }
  }
  if(!missing(b)){
    if(b < 0){
      
      stop('The argument "b" must be positive.')
    }
  }
  if(!missing(p)){
    if(p < 0){
      stop('The argument "p" must be positive.')
    }
  }

  if(p == 0 && b > 0 && a > 0){
      stop('If "p" equals to 0, "b" or "a" must be 0 otherwise there is identifiability problem.')
  }  
  if(b == 0 && a == 0){
    stop('The distribution is not defined for "a" and "b" equal to 0 simultaneously.')
  } 

  d <- ((a + b*(1+p)*(abs(x-mu)^p))*exp(-(x-mu)*(a+b*(abs(x-mu)^p)))) / ((exp(-(x-mu)*(a + b* (abs(x-mu)^p)))+1)^2) 
  
  d <- ifelse(is.nan(d), 0, d)
  
  return(d)
}

#' @rdname distrib
#' @export

qgenlog <- function(k, a = sqrt(2/pi), b = 0.5, p = 2, mu = 0, lower.tail = TRUE){
  
  if(!missing(a)){
    if(a < 0){
      stop('The argument "a" must be positive.')
    }
  }
  if(!missing(b)){
    if(b < 0){
      
      stop('The argument "b" must be positive.')
    }
  }
  if(!missing(p)){
    if(p < 0){
      stop('The argument "p" must be positive.')
    }
  }
  
  if(p == 0 && b > 0 && a > 0){
    stop('If "p" equals to 0, "b" or "a" must be 0 otherwise there is identifiability problem.')
  }  
  if(b == 0 && a == 0){
    stop('The distribution is not defined for "a" and "b" equal to 0 simultaneously.')
  }
  
  dgen_log <- function(x, a1 = a, b1 = b, p1 = p){
    
    d <- ((a1 + b1*(1+p1)*(abs(x-mu)^p1))*exp(-(x-mu)*(a1+b1*(abs(x-mu)^p1)))) / 
      ((exp(-(x-mu)*(a1 + b1* (abs(x-mu)^p1)))+1)^2) 
    d <- ifelse(is.nan(d), 0, d)

    return(d)
  }
  
  cont_dist <- distr::AbscontDistribution(d = dgen_log)
  
  qdist <- distr::q(cont_dist)
  
  ret <- ifelse(lower.tail == TRUE, qdist(k), qdist(1-k))
  
  return(ret)
  
}

#' @rdname distrib
#' @export

rgenlog <- function(n, a = sqrt(2/pi), b = 0.5, p = 2, mu = 0){
  
  if(!missing(a)){
    if(a < 0){
      stop('The argument "a" must be positive.')
    }
  }
  if(!missing(b)){
    if(b < 0){
      
      stop('The argument "b" must be positive.')
    }
  }
  if(!missing(p)){
    if(p < 0){
      stop('The argument "p" must be positive.')
    }
  }
  
  if(p == 0 && b > 0 && a > 0){
    stop('If "p" equals to 0, "b" or "a" must be 0 otherwise there is identifiability problem.')
  }  
  if(b == 0 && a == 0){
    stop('The distribution is not defined for "a" and "b" equal to 0 simultaneously.')
  }
  
  
  dgen_log <- function(x, a1 = a, b1 = b, p1 = p){
    
    d <- ((a1 + b1*(1+p1)*(abs(x-mu)^p1))*exp(-(x-mu)*(a1+b1*(abs(x-mu)^p1)))) / ((exp(-(x-mu)*(a1 + b1* (abs(x-mu)^p1)))+1)^2) 
    
    d <- ifelse(is.nan(d), 0, d)
    
    return(d)
  }
  
  cont_dist <- distr::AbscontDistribution(d = dgen_log)
  
  rdist <- distr::r(cont_dist)
  
  ret <- rdist(n)
  
  return(ret)
  
}
