
##########################

#' Simulating the Generalized logistic distribution
#'
#' Creating a simulation of the generalized logistic distribution maximum likelihood estimation of the parameters 
#' with parallelized processing code using the \code{foreach} package.
#' @param real.par the real parameters value of the distribution wich the random sample will be taken. It has to be a vector os length 4,
#' the parameters are the values of \code{c(a, b, p, location)} as listed in \code{rgenlog}, 
#' \code{location} can be omitted and will be set to 0. There are no default values.
#' @param init.par Initial values for the parameters to be optimized over in the following order \code{c(a, b, p, location)}.
#'  Can be an object returned by \code{genlog_slider}. There are no default values.
#' @param sample.size the sample size to be taken in each \code{k} simulation.
#' @param k  the number of simulations.
#' @param seed seed to be given to \code{set.seed()} function during the sampling process
#' @param threads the numbers of CPU threads to be used for parallel computing. If the threads 
#' number is higher than the available the maximum allowed will be used.
#' @keywords parallel, computing, maximum likelihood, mle, simulation, simulating, genlog_simu
#' 
#' @export
#' @examples 
#' genlog_simu(real.par = c(0.3, 0.9, 1.5, 0.0), init.par = c(0.9, 0.3, 0.2, 0.0), 
#'             sample.size = 100, k = 100, threads = 4, seed = 200) 
#' 
#' @usage 
#' genlog_simu(real.par, init.par, sample.size = 100,
#'             k = 1000, seed = 555, threads = 1)
#' 
#' @return 
#' It returns a data.frame with \code{k} rows (each simulation) and 7 columns with the following information:
#' \cr
#' \code{a, b,  p} and \code{location} are estimations using maximum likelihood estimation, for more info \code{help(genlogis_mle)}  
#' \cr
#' \code{sample.size} The sample size used for each \code{k} simulation.
#' \cr
#' \code{convergence} The estimation's convergence status.
#' \cr
#' \code{conv_message} Any aditional convergece status message.
#' \cr
#' 
#' @details 
#' 
#' The used distribution for this package is given by: \deqn{f(x) = ((a + b*(1+p)*(abs(x-location)^p))*exp(-(x-location)*(a+b*(|x-location|^p)))) / ((exp(-(x-location)*(a + b* (|x-location|^p)))+1)^2)}
#'  For more about the distribution use \code{help(dgenlog)}.
#' 
#' @references 
#' RATHIE, P. N., SWAMEE, P. K. \emph{On a new invertible generalized logistic distribution
#' approximation to normal distribution}, Technical Research Report in Statistics, 07/2006,
#' Dept. of Statistics, Univ. of Brasilia, Brasilia, Brazil. 2006.
#' 
#' WETSON, Steve, \emph{Using the foreach Package. Oct. 2015.
#' Consulted on August 2017 at https://cran.r-project.org/web/packages/foreach/vignettes/foreach.pdf} 
#' 
#' 

genlog_simu <- function(real.par, init.par, sample.size = 100,
                        k = 1000, seed = 555, threads = 1){

  if(length(real.par) < 3 | length(real.par) > 4 ){
    stop('Incorrect number of parameters: real.par = c(a,b,p,location)')
  }
  
  if(length(init.par) != 4 ){
    stop('Incorrect number of parameters: init.par = c(a,b,p,location)')
  }
  
  if(length(real.par) == 3){
    warning('Location parameter is set to 0')
    location = 0
  }
  
  
  if(length(real.par) == 4){
    location = real.par[4]
  }
  
  a = real.par[1]
  b = real.par[2]
  p = real.par[3]
  
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
  
  set.seed(seed)  
  cl <- snow::makeCluster(threads) #number of CPU cores
  doSNOW::registerDoSNOW(cl)
  results <- data.frame()
  
  results <- foreach::foreach(i = 1:k, .packages = c('genlogis'),
                     .combine = 'rbind', .inorder = FALSE) %dopar% {
                       
                       am1 <- genlogis::rgenlog(a = a, b = b,
                                      p = p, location = location, n = sample.size)
                       
                       mle1 <- genlogis::genlog_mle(init.par, data = am1)
                       
                       
                       ret <- rbind(c(mle1$par, sample.size, mle1$convergence, mle1$message))
                       
                       return(ret)          
                     }
  
  snow::stopCluster(cl)
  
  colnames(results) <- c('a', 'b', 'p', 'location', 'sample.size', 'convergence', 'conv_message')
  results <- as.data.frame(results)
  
  for(i in 1:7){
    if(i < 7){
      results[,i] <- as.numeric(levels(results[,i]))[results[,i]]
    }
    else{
      results[,i] <- as.character(levels(results[,i]))[results[,i]]
    }
  }
  
  return(results)
}