
##########################

#' The Generalized logistic distribution
#'
#' Density, distribution function, quantile function and random generation a generalized logistic distribution.
#' @param x,q vector of quantiles.
#' @param k vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required
#' @param a,b,p  parameters >= 0, with restrictions.*
#' @param location location parameter
#' @keywords d, p, q, r, genlogis
#' 
#' @export
#' @examples 
#' pgenlog(0.5) 
#' curve(dgenlog(x), xlim = c(-3,3)) 
#' 
#' @usage 
#' dgenlog(x, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#' pgenlog(q, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#' qgenlog(k, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#' rgenlog(n, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#'
#' 
#' @return 
#' \code{dgenlogis} gives the density, \code{pgenlogis} gives the distribution function, 
#' \code{qgenlogis} gives the quantile function, and \code{rgenlogis} generates random deviates.\cr
#'  
#'  The length of the result is determined by \code{n} for \code{rgenlogis}, and is the maximum of the lengths 
#'  of the numerical arguments for the other functions.
#' 
#' @details 
#' 
#' The used distribution for this package is given by: \deqn{f(x) = ((a + b*(1+p)*(abs(x-location)^p))*exp(-(x-location)*(a+b*(|x-location|^p)))) / ((exp(-(x-location)*(a + b* (|x-location|^p)))+1)^2)}
#'  
#' The \code{qgenlog()} returns values for P(X < x).\cr 
#' 
#' The default values for \code{a, b, p and location} produces a function with mean 0 and variance close to 1.\cr 
#' 
#' *Restrictions:\cr 
#' 
#' If \code{p} equals to 0, \code{b} or \code{a} must be 0 otherwise there is identifiability problem.\cr 
#' 
#' The distribution is not defined for \code{a} and \code{b} equal to 0 simultaneously.\cr 
#' @references 
#' RATHIE, P. N., SWAMEE, P. K. \emph{On a new invertible generalized logistic distribution
#' approximation to normal distribution}, Technical Research Report in Statistics, 07/2006,
#' Dept. of Statistics, Univ. of Brasilia, Brasilia, Brazil. 2006.

pgenlog <- function(q, a = sqrt(2/pi), b = 0.5, p = 2, location = 0){

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
  
  z <- (exp(-(q-location)*(a+b*(abs(q-location)^p)))+1)^(-1)
  
  return(z)
}


#' The Generalized logistic distribution
#'
#' Density, distribution function, quantile function and random generation a generalized logistic distribution.
#' @param x,q vector of quantiles.
#' @param k vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required
#' @param a,b,p  parameters >= 0, with restrictions.*
#' @param location location parameter
#' @keywords d, p, q, r, genlogis
#' 
#' @export
#' @examples 
#' pgenlog(0.5) 
#' curve(dgenlog(x), xlim = c(-3,3)) 
#' 
#' @usage 
#' dgenlog(x, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#' pgenlog(q, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#' qgenlog(k, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#' rgenlog(n, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#'
#' 
#' @return 
#' \code{dgenlogis} gives the density, \code{pgenlogis} gives the distribution function, 
#' \code{qgenlogis} gives the quantile function, and \code{rgenlogis} generates random deviates.\cr
#'  
#'  The length of the result is determined by \code{n} for \code{rgenlogis}, and is the maximum of the lengths 
#'  of the numerical arguments for the other functions.
#' 
#' @details 
#' 
#' The used distribution for this package is given by: \deqn{f(x) = ((a + b*(1+p)*(abs(x-location)^p))*exp(-(x-location)*(a+b*(|x-location|^p)))) / ((exp(-(x-location)*(a + b* (|x-location|^p)))+1)^2)}
#'  
#' The \code{qgenlog()} returns values for P(X < x).\cr 
#' 
#' The default values for \code{a, b, p and location} produces a function with mean 0 and variance close to 1.\cr 
#' 
#' *Restrictions:\cr 
#' 
#' If \code{p} equals to 0, \code{b} or \code{a} must be 0 otherwise there is identifiability problem.\cr 
#' 
#' The distribution is not defined for \code{a} and \code{b} equal to 0 simultaneously.\cr 
#' @references 
#' RATHIE, P. N., SWAMEE, P. K. \emph{On a new invertible generalized logistic distribution
#' approximation to normal distribution}, Technical Research Report in Statistics, 07/2006,
#' Dept. of Statistics, Univ. of Brasilia, Brasilia, Brazil. 2006.


dgenlog <- function(x, a = sqrt(2/pi), b = 0.5, p = 2, location = 0){
  
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
  
  d <- ((a + b*(1+p)*(abs(x-location)^p))*exp(-(x-location)*(a+b*(abs(x-location)^p)))) / ((exp(-(x-location)*(a + b* (abs(x-location)^p)))+1)^2) 
  
  return(d)
}

#' The Generalized logistic distribution
#'
#' Density, distribution function, quantile function and random generation a generalized logistic distribution.
#' @param x,q vector of quantiles.
#' @param k vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required
#' @param a,b,p  parameters >= 0, with restrictions.*
#' @param location location parameter
#' @keywords d, p, q, r, genlogis
#' 
#' @export
#' @examples 
#' pgenlog(0.5) 
#' curve(dgenlog(x), xlim = c(-3,3)) 
#' 
#' @usage 
#' dgenlog(x, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#' pgenlog(q, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#' qgenlog(k, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#' rgenlog(n, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#'
#' 
#' @return 
#' \code{dgenlogis} gives the density, \code{pgenlogis} gives the distribution function, 
#' \code{qgenlogis} gives the quantile function, and \code{rgenlogis} generates random deviates.\cr
#'  
#'  The length of the result is determined by \code{n} for \code{rgenlogis}, and is the maximum of the lengths 
#'  of the numerical arguments for the other functions.
#' 
#' @details 
#' 
#' The used distribution for this package is given by: \deqn{f(x) = ((a + b*(1+p)*(abs(x-location)^p))*exp(-(x-location)*(a+b*(|x-location|^p)))) / ((exp(-(x-location)*(a + b* (|x-location|^p)))+1)^2)}
#'  
#' The \code{qgenlog()} returns values for P(X < x).\cr 
#' 
#' The default values for \code{a, b, p and location} produces a function with mean 0 and variance close to 1.\cr 
#' 
#' *Restrictions:\cr 
#' 
#' If \code{p} equals to 0, \code{b} or \code{a} must be 0 otherwise there is identifiability problem.\cr 
#' 
#' The distribution is not defined for \code{a} and \code{b} equal to 0 simultaneously.\cr 
#' @references 
#' RATHIE, P. N., SWAMEE, P. K. \emph{On a new invertible generalized logistic distribution
#' approximation to normal distribution}, Technical Research Report in Statistics, 07/2006,
#' Dept. of Statistics, Univ. of Brasilia, Brasilia, Brazil. 2006.

qgenlog <- function(k, a = sqrt(2/pi), b = 0.5, p = 2, location = 0){
  
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
    
    d <- ((a1 + b1*(1+p1)*(abs(x-location)^p1))*exp(-(x-location)*(a1+b1*(abs(x-location)^p1)))) / ((exp(-(x-location)*(a1 + b1* (abs(x-location)^p1)))+1)^2) 
    
    return(d)
  }
  
  cont_dist <- distr::AbscontDistribution(d = dgen_log)
  
  qdist <- distr::q(cont_dist)
  
  return(qdist(k))
  
}

#' The Generalized logistic distribution
#'
#' Density, distribution function, quantile function and random generation a generalized logistic distribution.
#' @param x,q vector of quantiles.
#' @param k vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required
#' @param a,b,p  parameters >= 0, with restrictions.*
#' @param location location parameter
#' @keywords d, p, q, r, genlogis
#' 
#' @export
#' @examples 
#' pgenlog(0.5) 
#' curve(dgenlog(x), xlim = c(-3,3)) 
#' 
#' @usage 
#' dgenlog(x, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#' pgenlog(q, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#' qgenlog(k, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#' rgenlog(n, a = sqrt(2/pi), b = 0.5, p = 2, location = 0)
#'
#' 
#' @return 
#' \code{dgenlogis} gives the density, \code{pgenlogis} gives the distribution function, 
#' \code{qgenlogis} gives the quantile function, and \code{rgenlogis} generates random deviates.\cr
#'  
#'  The length of the result is determined by \code{n} for \code{rgenlogis}, and is the maximum of the lengths 
#'  of the numerical arguments for the other functions.
#' 
#' @details 
#' 
#' The used distribution for this package is given by: \deqn{f(x) = ((a + b*(1+p)*(abs(x-location)^p))*exp(-(x-location)*(a+b*(|x-location|^p)))) / ((exp(-(x-location)*(a + b* (|x-location|^p)))+1)^2)}
#'  
#' The \code{qgenlog()} returns values for P(X < x).\cr 
#' 
#' The default values for \code{a, b, p and location} produces a function with mean 0 and variance close to 1.\cr 
#' 
#' *Restrictions:\cr 
#' 
#' If \code{p} equals to 0, \code{b} or \code{a} must be 0 otherwise there is identifiability problem.\cr 
#' 
#' The distribution is not defined for \code{a} and \code{b} equal to 0 simultaneously.\cr 
#' @references 
#' RATHIE, P. N., SWAMEE, P. K. \emph{On a new invertible generalized logistic distribution
#' approximation to normal distribution}, Technical Research Report in Statistics, 07/2006,
#' Dept. of Statistics, Univ. of Brasilia, Brasilia, Brazil. 2006.

rgenlog <- function(n, a = sqrt(2/pi), b = 0.5, p = 2, location = 0){
  
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
    
    d <- ((a1 + b1*(1+p1)*(abs(x-location)^p1))*exp(-(x-location)*(a1+b1*(abs(x-location)^p1)))) / ((exp(-(x-location)*(a1 + b1* (abs(x-location)^p1)))+1)^2) 
    
    return(d)
  }
  
  cont_dist <- distr::AbscontDistribution(d = dgen_log)
  
  rdist <- distr::r(cont_dist)
  
  return(rdist(n))
  
}
