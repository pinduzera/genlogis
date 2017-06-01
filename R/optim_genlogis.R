#' Optimization for a generalized logistic distribution
#'
#' Maximum likehood estimation of parameters for a generalized logistic distribution.
#' 
#' @param parameters Initial values for the parameters to be optimized over in the following order \code{c(a, b, p, location)},
#'  \code{location} can be omitted and will be equaled to 0.
#' @param data This is the the data to be utilized for the estimation
#' @keywords d, p, q, r
#' 
#' @export
#' @examples 
#' datas <- rgenlog(10000, 1,2,1, 0)
#' genlog_mle(c(.7,1.6, .8, 0),datas) 
#' 
#' @usage 
#' genlog_mle(parameters, data)
#' 
#' @details 
#' Maximum likehood estimation of parameters for the distribution proposed in this package.\cr
#' This function is an application of \code{optim} as a particular case needed for this distribution using the method "\code{L-BFGS-B}".  \cr
#' For more information of the output check \code{help(optim)}.\cr
#' 
#' The used distribution for this package is given by: 
#' \deqn{f(x) = ((a + b*(1+p)*(abs(x-location)^p))*exp(-(x-location)*(a+b*(|x-location|^p)))) / ((exp(-(x-location)*(a + b* (|x-location|^p)))+1)^2)}
#' 
#' \code{help(dgenlog)} for parameters restrictions.\cr


genlog_mle <- function(parameters, data){

genlogis.loglikehood <- function(param = c(sqrt(2/pi),0.5, 2, 0), x){
  
  if(length(param) < 3 | length(param) > 4 ){
    stop('Incorrect number of parameters: param = c(a,b,p,location)')
  }
  
  if(length(param) == 3){
    warning('Location parameter is set to 0')
    location = 0
  }
  
  if(length(param) == 4){
    location = param[4]
  }
  
  a = param[1]
  b = param[2]
  p = param[3]
  
  if(!missing(a)){
    if(a < 0){
      stop('The argument "a" must be positive.')
    }
  }
  if(!missing(a)){
    if(b < 0){
      
      stop('The argument "b" must be positive.')
    }
  }
  if(!missing(a)){
    if(p < 0){
      stop('The argument "p" must be positive.')
    }
  }
  
  if(p == 0 && b > 0 && a > 0){
    stop('If "p" equals to 0, "b" or "a" must be 
         0 otherwise there is identifiability problem.')
  }  
  if(b == 0 && a == 0){
    stop('The distribution is not defined for "a" 
         and "b" equal to 0 simultaneously.')
  }

  z <- sum(log((a+b*(1+p)*abs((x-location))^p ) * exp(-((x-location)*(a+b*abs((x-location))^p)))) -
     log((exp(-((x-location)*(a+b*abs((x-location))^p))) + 1)^2))
    
  return(-z)
}


op <- optim(par=parameters, fn = genlogis.loglikehood, x=data,
      lower = c(0,0,0,-Inf), upper = c(Inf,Inf,Inf,Inf),
      method = 'L-BFGS-B') 

return(op)

}
