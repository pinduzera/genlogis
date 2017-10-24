### this is a function of internal use for MLE in optim_mle() function

genlogis.likeli.gr <- function(param = c(sqrt(2/pi),0.5, 2, 0), x){

  if(length(param) == 3){
    #warning('mu parameter is set to 0')
    mu = 0
  }
  
  if(length(param) == 4){
    mu = param[4]
  }
  
  a = param[1]
  b = param[2]
  p = param[3]
  
  gr.a <- ((exp((x-mu)*(b*abs(x-mu)^p+a))*((mu-x)*exp(-(x-mu)*(b*abs(x-mu)^p+a))*(b*(p+1)*abs(x-mu)^p+a)+exp(-(x-mu)*(b*abs(x-mu)^p+a))))/(b*(p+1)*abs(x-mu)^p+a))-
          ((2*(mu-x)*exp(-(x-mu)*(b*abs(x-mu)^p+a)))/(exp(-(x-mu)*(b*abs(x-mu)^p+a))+1))
  
  gr.b <- ((exp((x-mu)*(b*abs(x-mu)^p+a))*((p+1)*exp(-(x-mu)*(b*abs(x-mu)^p+a))*abs(x-mu)^p-(x-mu)*exp(-(x-mu)*(b*abs(x-mu)^p+a))*abs(x-mu)^p*(b*(p+1)*abs(x-mu)^p+a)))/(b*(p+1)*abs(x-mu)^p+a))+
          ((2*(x-mu)*exp(-(x-mu)*(b*abs(x-mu)^p+a))*abs(x-mu)^p)/(exp(-(x-mu)*(b*abs(x-mu)^p+a))+1))
  
  gr.p <- ((exp((x-mu)*(b*abs(x-mu)^p+a))*(exp(-(x-mu)*(b*abs(x-mu)^p+a))*(b*(p+1)*abs(x-mu)^p*log(abs(x-mu))+b*abs(x-mu)^p)-b*(x-mu)*exp(-(x-mu)*(b*abs(x-mu)^p+a))*abs(x-mu)^p*(b*(p+1)*abs(x-mu)^p+a)*log(abs(x-mu))))/(b*(p+1)*abs(x-mu)^p+a))+
          (2*b*(x-mu)*exp(-(x-mu)*(b*abs(x-mu)^p+a))*abs(x-mu)^p*log(abs(x-mu)))/(exp(-(x-mu)*(b*abs(x-mu)^p+a))+1)
  
  gr.mu <- ((exp((x-mu)*(b*abs(x-mu)^p+a))*(exp(-(x-mu)*(b*abs(x-mu)^p+a))*(b*p*abs(x-mu)^p+b*abs(x-mu)^p+a)*(b*(p+1)*abs(x-mu)^p+a)-(b*p*(p+1)*exp(-(x-mu)*(b*abs(x-mu)^p+a))*abs(x-mu)^p)/(x-mu)))/(b*(p+1)*abs(x-mu)^p+a))-
            ((2*exp(-(x-mu)*(b*abs(x-mu)^p+a))*(b*p*abs(x-mu)^p+b*abs(x-mu)^p+a))/(exp(-(x-mu)*(b*abs(x-mu)^p+a))+1))
  
  grad <- colSums(cbind(gr.a, gr.b, gr.p, gr.mu), na.rm = TRUE)

  return(-grad)
}

### this is a function of internal use for MLE in optim_mle() function

genlogis.as.likeli.gr <- function(param = c(sqrt(2/pi),0.5, 2, 0, .5), x){
  
  if(length(param) != 5){
    stop('Incorrect number of parameters: param = c(a, b, p, mu, skew)')
  }
  
  
  a <- param[1]
  b <- param[2]
  p <- param[3]
  mu <- param[4]
  skew <- param[5]
  
  gr.a <- 1/(b*(p+1)*abs(x-mu)^p+a)+
          (skew*(x-mu)*exp(-skew*(x-mu)*(abs(skew)^p*b*abs(x-mu)^p+a)))/
          (exp(-skew*(x-mu)*(abs(skew)^p*b*abs(x-mu)^p+a))+1)-
          (2*(mu-x)*exp((mu-x)*(b*abs(x-mu)^p+a)))/(exp((mu-x)*(b*abs(x-mu)^p+a))+1)-x+mu

  gr.b <- ((p+1)*abs(x-mu)^p)/(b*(p+1)*abs(x-mu)^p+a)+
          (skew*abs(skew)^p*(x-mu)*exp(-skew*(x-mu)*(abs(skew)^p*b*abs(x-mu)^p+a))*abs(x-mu)^p)/
         (exp(-skew*(x-mu)*(abs(skew)^p*b*abs(x-mu)^p+a))+1)-
         (2*(mu-x)*exp((mu-x)*(b*abs(x-mu)^p+a))*abs(x-mu)^p)/
          (exp((mu-x)*(b*abs(x-mu)^p+a))+1)-(x-mu)*abs(x-mu)^p

  gr.p <- (b*(p+1)*abs(x-mu)^p*log(abs(x-mu))+b*abs(x-mu)^p)/
          (b*(p+1)*abs(x-mu)^p+a)+
          (skew*(x-mu)*exp(-skew*(x-mu)*(abs(skew)^p*b*abs(x-mu)^p+a))*(abs(skew)^p*b*abs(x-mu)^p*log(abs(x-mu))+abs(skew)^p*log(abs(skew))*b*abs(x-mu)^p))/
          (exp(-skew*(x-mu)*(abs(skew)^p*b*abs(x-mu)^p+a))+1)-
          (2*b*(mu-x)*exp((mu-x)*(b*abs(x-mu)^p+a))*abs(x-mu)^p*log(abs(x-mu)))/
          (exp((mu-x)*(b*abs(x-mu)^p+a))+1)-b*(x-mu)*abs(x-mu)^p*log(abs(x-mu))

  gr.mu <- -(b*p*(p+1)*abs(x-mu)^p)/((x-mu)*(b*(p+1)*abs(x-mu)^p+a))-
            (exp(-skew*(x-mu)*(abs(skew)^p*b*abs(x-mu)^p+a))*(skew*(abs(skew)^p*b*abs(x-mu)^p+a)+skew*abs(skew)^p*b*p*abs(x-mu)^p))/
            (exp(-skew*(x-mu)*(abs(skew)^p*b*abs(x-mu)^p+a))+1)-
            (2*exp((mu-x)*(b*abs(x-mu)^p+a))*(-(b*p*(mu-x)*abs(x-mu)^p)/(x-mu)+b*abs(x-mu)^p+a))/
            (exp((mu-x)*(b*abs(x-mu)^p+a))+1)+b*p*abs(x-mu)^p+b*abs(x-mu)^p+a

  gr.skew <- -(exp(-skew*(x-mu)*(abs(skew)^p*b*abs(x-mu)^p+a))*(-(x-mu)*(abs(skew)^p*b*abs(x-mu)^p+a)-abs(skew)^p*b*p*(x-mu)*abs(x-mu)^p))/
              (exp(-skew*(x-mu)*(abs(skew)^p*b*abs(x-mu)^p+a))+1)

  grad <- colSums(cbind(gr.a, gr.b, gr.p, gr.mu, gr.skew), na.rm = TRUE)
  
  return(-grad)
}

