
##########################

#' The Generalized logistic distribution
#'
#' Density, distribution function, quantile function and random generation a generalized logistic distribution.
#' @param x,q vector of quantiles.
#' @param k vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required
#' @param a,b,p  parameters > 0.
#' 
#' @keywords d, p, q, r
#' @export
#' @examples
#' pdgenlog(0.5)
#' curve(dgenlog(x), c(-3,3))
#' 
#' @usage 
#' dgenlog(x, a, b, p)
#' pgenlog(q, a, b, p)
#' qgenlog(k, a, b, p)
#' rgenlog(n, a, b, p)
#' 
#' @details 
#' 
#' The used distribution for this package is given by: \deqn{((a + b*(1+p)*(abs(x)^p))*exp(-x*(a+b*(|x|^p)))) / ((exp(-x*(a + b* (|x|^p)))+1)^2)}
#' 
#' The \code{qgenlog()} returns values for P(X < x).
#' 
#' The default values for \code{a, b and p} are \code{sqrt(2/pi), 0.5 and 2} which produces a functions with mean 0 and variance close to 1.

pgenlog <- function(q, a = sqrt(2/pi), b = 0.5, p = 2){
  
  if(!missing(a)){
    if(a <= 0){
      stop('The argument "a" must be positive')
    }
  }
  if(!missing(a)){
    if(b <= 0){
      
      stop('The argument "b" must be positive')
    }
  }
  if(!missing(a)){
    if(p <= 0){
      stop('The argument "p" must be positive')
    }
  }
  
  z <- (exp(-q*(a+b*(abs(q)^p)))+1)^(-1)
  
  return(z)
}


#' The Generalized logistic distribution
#'
#' Density, distribution function, quantile function and random generation a generalized logistic distribution.
#' @param x,q vector of quantiles.
#' @param k vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required
#' @param a,b,p  parameters > 0.
#' 
#' @keywords d, p, q, r
#' @export
#' @examples 
#' pdgenlog(0.5)
#' curve(dgenlog(x), c(-3,3))
#' 
#' @usage 
#' dgenlog(x, a, b, p)
#' pgenlog(q, a, b, p)
#' qgenlog(k, a, b, p)
#' rgenlog(n, a, b, p)
#' 
#' @details 
#' 
#' The used distribution for this package is given by: \deqn{((a + b*(1+p)*(abs(x)^p))*exp(-x*(a+b*(|x|^p)))) / ((exp(-x*(a + b* (|x|^p)))+1)^2)}
#' 
#' The \code{qgenlog()} returns values for P(X < x).
#' 
#' The default values for \code{a, b and p} are \code{sqrt(2/pi), 0.5 and 2} which produces a functions with mean 0 and variance close to 1.


dgenlog <- function(x, a = sqrt(2/pi), b = 0.5, p = 2){
  
  if(!missing(a)){
    if(a <= 0){
      stop('The argument "a" must be positive')
    }
  }
  if(!missing(a)){
    if(b <= 0){
      
      stop('The argument "b" must be positive')
    }
  }
  if(!missing(a)){
    if(p <= 0){
      stop('The argument "p" must be positive')
    }
  }
  
  d <- ((a + b*(1+p)*(abs(x)^p))*exp(-x*(a+b*(abs(x)^p)))) / ((exp(-x*(a + b* (abs(x)^p)))+1)^2) 
  
  return(d)
}


#' The Generalized logistic distribution
#'
#' Density, distribution function, quantile function and random generation a generalized logistic distribution.
#' @param x,q vector of quantiles.
#' @param k vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required
#' @param a,b,p  parameters > 0.
#' 
#' @keywords d, p, q, r
#' @export
#' @examples 
#' pdgenlog(0.5) 
#' curve(dgenlog(x), c(-3,3))
#' 
#' @usage 
#' dgenlog(x, a, b, p)
#' pgenlog(q, a, b, p)
#' qgenlog(k, a, b, p)
#' rgenlog(n, a, b, p)
#' 
#' @details 
#' 
#' The used distribution for this package is given by: \deqn{((a + b*(1+p)*(abs(x)^p))*exp(-x*(a+b*(|x|^p)))) / ((exp(-x*(a + b* (|x|^p)))+1)^2)}
#' 
#' The \code{qgenlog()} returns values for P(X < x).
#' 
#' The default values for \code{a, b and p} are \code{sqrt(2/pi), 0.5 and 2} which produces a functions with mean 0 and variance close to 1.



qgenlog <- function(k, a = sqrt(2/pi), b = 0.5, p = 2){
  
  if(!missing(a)){
    if(a <= 0){
      stop('The argument "a" must be positive')
    }
  }
  if(!missing(a)){
    if(b <= 0){
      
      stop('The argument "b" must be positive')
    }
  }
  if(!missing(a)){
    if(p <= 0){
      stop('The argument "p" must be positive')
    }
  }
  
  dgen_log <- function(x, a1 = a, b1 = b, p1 = p){
    
    d <- ((a1 + b1*(1+p1)*(abs(x)^p1))*exp(-x*(a1+b1*(abs(x)^p1)))) / ((exp(-x*(a1 + b1* (abs(x)^p1)))+1)^2) 
    
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
#' @param a,b,p  parameters > 0.
#' 
#' @keywords d, p, q, r
#' @export
#' @examples 
#' pdgenlog(0.5) 
#' curve(dgenlog(x), c(-3,3)) 
#' 
#' @usage 
#' dgenlog(x, a, b, p)
#' pgenlog(q, a, b, p)
#' qgenlog(k, a, b, p)
#' rgenlog(n, a, b, p)
#' 
#' @details 
#' 
#' The used distribution for this package is given by: \deqn{((a + b*(1+p)*(abs(x)^p))*exp(-x*(a+b*(|x|^p)))) / ((exp(-x*(a + b* (|x|^p)))+1)^2)}
#' 
#' The \code{qgenlog()} returns values for P(X < x).
#' 
#' The default values for \code{a, b and p} are \code{sqrt(2/pi), 0.5 and 2} which produces a functions with mean 0 and variance close to 1.


rgenlog <- function(n, a = sqrt(2/pi), b = 0.5, p = 2){
  
  if(!missing(a)){
    if(a <= 0){
      stop('The argument "a" must be positive')
    }
  }
  if(!missing(a)){
    if(b <= 0){
      
      stop('The argument "b" must be positive')
    }
  }
  if(!missing(a)){
    if(p <= 0){
      stop('The argument "p" must be positive')
    }
  }
  
  
  dgen_log <- function(x, a1 = a, b1 = b, p1 = p){
    
    d <- ((a1 + b1*(1+p1)*(abs(x)^p1))*exp(-x*(a1+b1*(abs(x)^p1)))) / ((exp(-x*(a1 + b1* (abs(x)^p1)))+1)^2) 
    
    return(d)
  }
  
  cont_dist <- distr::AbscontDistribution(d = dgen_log)
  
  rdist <- distr::r(cont_dist)
  
  return(rdist(n))
  
}
