##########################

#' Simulating the Generalized logistic distribution
#'
#' Creating a simulation of the generalized logistic distribution maximum likelihood estimation of the parameters 
#' with parallelized processing code using the \code{foreach} package.
#' @param real.par the real parameters value of the distribution wich the random sample will be taken. It has to be a vector of length 4,
#' the parameters are the values of \code{c(a, b, p, mu)} as listed in \code{rgenlog}, 
#' \code{mu} can be omitted and will be set to 0. There are no default values.
#' @param init.par Initial values for the parameters to be optimized over in the following order \code{c(a, b, p, mu)}.
#'  Can be an object returned by \code{genlog_slider}. There are no default values.
#' @param sample.size the sample size to be taken in each \code{k} simulation.
#' @param k  the number of simulations.
#' @param seed seed to be given to \code{set.seed()} function during the sampling process
#' @param threads the numbers of CPU threads to be used for parallel computing. If the threads 
#' number is higher than the available the maximum allowed will be used.
#' @param progress.bar show progress bar for each thread during simulations, default value \code{TRUE}.
#' @keywords genlogis
#' 
#' @export
#' @import foreach
#' @examples 
#'  
#' genlog_simu(real.par = c(0.3, 0.9, 1.5, 0.0), init.par = c(0.9, 0.3, 0.2, 0.0), 
#'             sample.size = 100, k = 50, threads = 2, seed = 200) 
#' 
#' @usage 
#' genlog_simu(real.par, init.par, sample.size = 100,
#'             k = 1000, seed = 555, threads = 1, progress.bar = T)
#' 
#' @return 
#' It returns a data.frame with \code{k} rows (each simulation) and 7 columns with the following information:
#' \cr
#' \code{a, b,  p} and \code{mu} are estimations using maximum likelihood estimation, for more info \code{help(genlogis_mle)}  
#' \cr
#' \code{sample.size} The sample size used for each \code{k} simulation.
#' \cr
#' \code{convergence} The estimation's convergence status.
#' \cr
#' 
#' @details 
#' 
#' The used distribution for this package is given by: \deqn{f(x) = ((a + b*(1+p)*(abs(x-mu)^p))*exp(-(x-mu)*(a+b*(|x-mu|^p)))) / ((exp(-(x-mu)*(a + b* (|x-mu|^p)))+1)^2)}
#'  For more about the distribution use \code{help(dgenlog)}.
#' 
#' @references 
#' 
#' Rathie, P. N. and Swamee, P. K. (2006) \emph{On a new invertible generalized logistic distribution
#' approximation to normal distribution}, Technical Research Report in Statistics, 07/2006,
#' Dept. of Statistics, Univ. of Brasilia, Brasilia, Brazil.


genlog_simu <- function(real.par, init.par, sample.size = 100,
                         k = 1000, seed = 555, threads = 1, progress.bar = T){
  
  if(length(real.par) < 3 | length(real.par) > 4 ){
    stop('Incorrect number of parameters: real.par = c(a,b,p,mu)')
  }
  
  if(length(init.par) != 4 ){
    stop('Incorrect number of parameters: init.par = c(a,b,p,mu)')
  }
  
  if(length(real.par) == 3){
    warning('mu parameter is set to 0')
    mu = 0
  }
  
  
  if(length(real.par) == 4){
    mu = real.par[4]
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
  
  core <- function(){
    
    am1 <- rgenlog(a = a, b = b,
                             p = p, mu = mu, n = sample.size)
    
    mle1 <- genlog_mle(init.par, data = am1)
    
    
    ret <- rbind(c(mle1$par, sample.size, mle1$convergence))
    
    if(ret[,6] != 0){
      ret <- core()
    }
    return(ret)          
  }

  set.seed(seed)  
  
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  
  results <- data.frame()
  
  results <- foreach::foreach(i = 1:k, .packages = c('genlogis', 'tcltk'),
                              .combine = 'rbind', .inorder = T) %dopar% {
                                if(progress.bar == T){
                                if(!exists("pb")) pb <- tcltk::tkProgressBar("Parallel task", min=1, max=k)
                                tcltk::setTkProgressBar(pb, i)
                                }
                                core()        
                              }
  
  parallel::stopCluster(cl)
  
  colnames(results) <- c('a', 'b', 'p', 'mu', 'sample.size', 'convergence')
  results <- as.data.frame(results)
  i = 1
  for(i in 1:6){
    results[,i] <- as.numeric((results[,i]))
  }
  
  return(results)
}

##########################

#' Simulating the Generalized logistic distribution with skewness
#'
#' Creating a simulation of the generalized logistic distribution with skewness maximum likelihood estimation of the parameters 
#' with parallelized processing code using the \code{foreach} package.
#' @param real.par the real parameters value of the distribution wich the random sample will be taken. It has to be a vector of length 5,
#' the parameters are the values of \code{c(a, b, p, mu)} as listed in \code{rgenlog}, 
#' \code{mu} can be omitted and will be set to 0. There are no default values.
#' @param init.par Initial values for the parameters to be optimized over in the following order \code{c(a, b, p, mu, skew)}.
#'  Can be an object returned by \code{genlog_slider}. There are no default values.
#' @param sample.size the sample size to be taken in each \code{k} simulation.
#' @param k  the number of simulations.
#' @param seed seed to be given to \code{set.seed()} function during the sampling process
#' @param threads the numbers of CPU threads to be used for parallel computing. If the threads 
#' number is higher than the available the maximum allowed will be used.
#' @param progress.bar show progress bar for each thread during simulations, default value \code{TRUE}.
#' @keywords genlogis
#' 
#' @export
#' @examples 
#' genlog_simu_sk(real.par = c(0.3, 0.9, 1.5, 0.0, .9), init.par = c(0.9, 0.3, 0.2, 0.0, .9), 
#'             sample.size = 100, k = 50, threads = 2, seed = 200) 
#' 
#' @usage 
#' genlog_simu_sk(real.par, init.par, sample.size = 100,
#'             k = 1000, seed = 555, threads = 1, progress.bar = T)
#' 
#' @return 
#' It returns a data.frame with \code{k} rows (each simulation) and 7 columns with the following information:
#' \cr
#' \code{a, b,  p} and \code{mu} are estimations using maximum likelihood estimation, for more info \code{help(genlogis_mle)}  
#' \cr
#' \code{sample.size} The sample size used for each \code{k} simulation.
#' \cr
#' \code{convergence} The estimation's convergence status.
#' \cr
#' 
#' @details 
#' 
#' The used distribution for this package is given by: \deqn{f(x) = 2*((a + b*(1+p)*(abs(x-mu)^p))*exp(-(x-mu)*(a+b*(abs(x-mu)^p))))/ 
#'    ((exp(-(x-mu)*(a + b* (abs(x-mu)^p)))+1)^2) *
#'    ((exp(-(skew*(x-mu))*(a+b*(abs(skew*(x-mu))^p)))+1)^(-1)) }
#' 
#' @references 
#' 
#' Rathie, P. N. and Swamee, P. K (2006) \emph{On a new invertible generalized logistic distribution
#' approximation to normal distribution}, Technical Research Report in Statistics, 07/2006,
#' Dept. of Statistics, Univ. of Brasilia, Brasilia, Brazil. 
#' 
#' Azzalini, A. (1985) \emph{A class of distributions which includes the normal ones}. Scandinavian Journal of Statistics.

genlog_simu_sk <- function(real.par, init.par, sample.size = 100,
                        k = 1000, seed = 555, threads = 1, progress.bar = T){
  
  if(length(real.par) !=5 ){
    stop('Incorrect number of parameters: real.par = c(a, b, p, mu, skew)')
  }
  
  a <- real.par[1]
  b <- real.par[2]
  p <- real.par[3]
  mu <- real.par[4]
  skew <- real.par[5]
  
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
  
  core <- function(){
    
    am1 <- rgenlog_sk(a = a, b = b,
                   p = p, mu = mu, skew = skew,
                   n = sample.size)
    
    mle1 <- genlog_mle_sk(init.par, data = am1)
    
    
    ret <- rbind(c(mle1$par, sample.size, mle1$convergence))
    
    if(ret[,7] != 0){
      ret <- core()
    }
    return(ret)          
  }
  
  set.seed(seed)  
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  
  results <- data.frame()
  
  results <- foreach::foreach(i = 1:k, .packages = c('genlogis', 'tcltk'),
                              .combine = 'rbind', .inorder = T) %dopar% {
                                if(progress.bar == T){
                                  if(!exists("pb")) pb <- tcltk::tkProgressBar("Parallel task", min=1, max=k)
                                  tcltk::setTkProgressBar(pb, i)
                                }
                                core()        
                              }
  
  parallel::stopCluster(cl)
  
  colnames(results) <- c('a', 'b', 'p', 'mu','skew', 'sample.size', 'convergence')
  results <- as.data.frame(results)
  i = 1
  for(i in 1:7){
    results[,i] <- as.numeric((results[,i]))
  }
  
  return(results)
  }