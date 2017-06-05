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
#' datas <- rgenlog(10000, 1.5,2,2, 0)
#' genlog_mle(c(.5,1.6, 1.5, 0),datas)
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


op <- optim(par=parameters, fn = genlogis.loglikehood, x=data,
      lower = c(0,0,0,-Inf), upper = c(Inf,Inf,Inf,Inf),
      method = 'L-BFGS-B') 

return(op)

}

