### this is a function of internal use for MLE in optim_mle_as() function

genlogis.as.loglikehood <- function(param = c(sqrt(2/pi),0.5, 2, 0, .5), x){
  
  if(length(param) != 5){
    stop('Incorrect number of parameters: param = c(a, b, p, mu, skew)')
  }
  
  a <- param[1]
  b <- param[2]
  p <- param[3]
  mu <- param[4]
  skew <- param[5]
  
  
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
  
  if(skew > 1 | skew < -1){
    stop('The skew parameter must be in the interval (-1,1).')
  }
  
  z <- -sum(log(dgenlog_as(x, a, b, p, mu, skew)))
  
  return(z)
  }

