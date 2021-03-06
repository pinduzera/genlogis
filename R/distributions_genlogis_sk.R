
##########################

#' The Generalized logistic distribution with skewness
#'
#' Density, distribution function, quantile function and random generation a generalized logistic distribution with skewness.
#' @param x,q vector of quantiles.
#' @param k vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required
#' @param a,b,p  parameters \eqn{\le 0}, with restrictions.*
#' @param mu mu parameter
#' @param skew skewness parameter limited to the interval (-1, 1)
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @keywords genlogis
#' 
#' @export
#' @examples 
#' pgenlog_sk(0.5) 
#' curve(dgenlog_sk(x), xlim = c(-3,3)) 
#' 
#' rgenlog_sk(100) 
#' 
#' qgenlog_sk(0.95)
#' 
#' 
#' @name distrib_sk
#'  
#' @return 
#' \code{dgenlog_sk} gives the density, \code{pgenlog_sk} gives the distribution function, 
#' \code{qgenlog_sk} gives the quantile function, and \code{rgenlog_sk} generates random deviates.\cr
#'  
#'  The length of the result is determined by \code{n} for \code{rgenlog_sk}, and is the maximum of the lengths 
#'  of the numerical arguments for the other functions.
#' 
#' @details 
#' 
#' The used distribution for this package is given by: \deqn{f(x) = 2*((a + b*(1+p)*(abs(x-mu)^p))*exp(-(x-mu)*(a+b*(abs(x-mu)^p))))/ 
#'    ((exp(-(x-mu)*(a + b* (abs(x-mu)^p)))+1)^2) *
#'    ((exp(-(skew*(x-mu))*(a+b*(abs(skew*(x-mu))^p)))+1)^(-1)) }
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
#' 
#' Azzalini, A. (1985) \emph{A class of distributions which includes the normal ones}. Scandinavian Journal of Statistics.

#' @rdname distrib_sk
#' @export

pgenlog_sk <- function(q, a = sqrt(2/pi), b = 0.5, p = 2, mu = 0, skew = .5, lower.tail = TRUE){
  
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
  
  if(skew > 1 | skew < -1){
    stop('The skew parameter must be in the interval (-1,1).')
  }
  
  dgen_log <- function(x){
    
    d <- 2*((a + b*(1+p)*(abs(x-mu)^p))*exp(-(x-mu)*(a+b*(abs(x-mu)^p)))) / 
      ((exp(-(x-mu)*(a + b* (abs(x-mu)^p)))+1)^2) *
      ((exp(-(skew*(x-mu))*(a+b*(abs(skew*(x-mu))^p)))+1)^(-1))  
    
    d <- ifelse(is.nan(d), 0, d)
    
    return(d)
  }
  
  cont_dist <- distr::AbscontDistribution(d = dgen_log)
  
  pdist <- distr::p(cont_dist)
  
  ret <- ifelse(lower.tail == TRUE, pdist(q), 1-pdist(q))
  
  return(ret)
}


#' @rdname distrib_sk
#' @export

dgenlog_sk <- function(x, a = sqrt(2/pi), b = 0.5, p = 2, mu = 0, skew = .5){
  
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
  if(skew > 1 | skew < -1){
    stop('The skew parameter must be in the interval (-1,1).')
  }

  d <- 2*((a + b*(1+p)*(abs(x-mu)^p))*exp(-(x-mu)*(a+b*(abs(x-mu)^p)))) / 
    ((exp(-(x-mu)*(a + b* (abs(x-mu)^p)))+1)^2) *
    ((exp(-(skew*(x-mu))*(a+b*(abs(skew*(x-mu))^p)))+1)^(-1))
  
  d <- ifelse(is.nan(d), 0, d)
  
  return(d)
  
}

#' @rdname distrib_sk
#' @export

qgenlog_sk <- function(k, a = sqrt(2/pi), b = 0.5, p = 2, mu = 0, skew = .5, lower.tail = TRUE){
  
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
  if(skew > 1 | skew < -1){
    stop('The skew parameter must be in the interval (-1,1).')
  }
  
  dgen_log <- function(x){
    
    d <- 2*((a + b*(1+p)*(abs(x-mu)^p))*exp(-(x-mu)*(a+b*(abs(x-mu)^p)))) / 
      ((exp(-(x-mu)*(a + b* (abs(x-mu)^p)))+1)^2) *
      ((exp(-(skew*(x-mu))*(a+b*(abs(skew*(x-mu))^p)))+1)^(-1))  
    
    d <- ifelse(is.nan(d), 0, d)
    
    return(d)
  }
  
  cont_dist <- distr::AbscontDistribution(d = dgen_log)
  
  qdist <- distr::q(cont_dist)
  
  ret <- ifelse(lower.tail == TRUE, qdist(k), qdist(1-k))
  
  return(ret)
  
}

#' @rdname distrib_sk
#' @export

rgenlog_sk <- function(n, a = sqrt(2/pi), b = 0.5, p = 2, mu = 0, skew = 0.5){
  
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
  if(skew > 1 | skew < -1){
    stop('The skew parameter must be in the interval (-1,1).')
  }
  
  dgen_log <- function(x){
    
    d <- 2*((a + b*(1+p)*(abs(x-mu)^p))*exp(-(x-mu)*(a+b*(abs(x-mu)^p)))) / 
      ((exp(-(x-mu)*(a + b* (abs(x-mu)^p)))+1)^2) *
      ((exp(-(skew*(x-mu))*(a+b*(abs(skew*(x-mu))^p)))+1)^(-1))  
    
    d <- ifelse(is.nan(d), 0, d)
    
    return(d)
  }
  
  cont_dist <- distr::AbscontDistribution(d = dgen_log)
  
  rdist <- distr::r(cont_dist)
  
  ret <- rdist(n)
  
  return(ret)
  
}

