#' Optimization for a generalized logistic distribution
#'
#' Maximum likehood estimation of parameters for a generalized logistic distribution.
#' 
#' @param parameters Initial values for the parameters to be optimized over in the following order \code{c(a, b, p, location)},
#'  \code{location} can be omioped and will be equaled to 0.
#' @param data This is the the data to be utilized for the estimation
#' @param alpha type 1 error parameter for confidente interval
#' @keywords d, p, q, r, genlogis
#' 
#' @export
#' @examples 
#' datas <- rgenlog(10000, 1.5,2,2, 0)
#' genlog_mle2(c(.5,1.6, 1.5, 0),datas, alpha = 0.05)
#' 
#' @usage 
#' genlog_mle2(parameters, data, alpha = 0.05)
#' 
#' @return 
#' Same as \code{optim} with extra:
#' \code{bounds} Confidence interval for maximum likelihood parameter estimatives with of 1-\code(alpha) confidence 
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


genlog_mle2 <- function(parameters, data, alpha = 0.05){
  
  
  op <- optim(par=parameters, fn = genlogis.loglikehood, x=data,
              lower = c(0.01,0.01,0.01,-Inf), upper = c(Inf,Inf,Inf,Inf),
              method = 'L-BFGS-B', hessian = T) 
  
  fisher_info <- solve(op$hessian)
  var_matrix <- diag(fisher_info)
  var_matrix <- sqrt(var_matrix)
  upper <- op$par + qnorm(1-alpha) * var_matrix
  lower <- op$par - qnorm(1-alpha) * var_matrix
  interval <- data.frame(value=op$par, lower=lower, upper=upper)
  op$bounds <-  as.matrix(interval)
  
  return(op)
  
}

