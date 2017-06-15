### this is a function of internal use for MLE in optim_mle() function

genlogis.loglikehood <- function(param = c(sqrt(2/pi),0.5, 2, 0), x){
  
  if(length(param) < 3 | length(param) > 4 ){
    stop('Incorrect number of parameters: param = c(a,b,p,location)')
  }
  
  if(length(param) == 3){
    #warning('Location parameter is set to 0')
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
    stop('If "p" equals to 0, "b" or "a" must be 
         0 otherwise there is identifiability problem.')
  }  
  if(b == 0 && a == 0){
    stop('The distribution is not defined for "a" 
         and "b" equal to 0 simultaneously.')
  }
  
  z <- sum(log((a+b*(1+p)*abs((x-location))^p ) * exp(-((x-location)*(a+b*abs((x-location))^p))))) -
             sum(log((exp(-((x-location)*(a+b*abs((x-location))^p))) + 1)^2))
  if(!is.finite(z)){
    z <- 1e+20
  }
  
  return(-z)
}

