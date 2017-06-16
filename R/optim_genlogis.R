#' Optimization for a generalized logistic distribution
#'
#' Maximum likehood estimation of parameters for a generalized logistic distribution.
#' 
#' @param parameters Initial values for the parameters to be optimized over in the following order \code{c(a, b, p, location)},
#'  \code{location} can be omitted and will be set to 0.
#' @param data This is the the data to be utilized for the estimation.
#' @param hessian logical value that returns hessian, also returns the parameters estimation's confidence interval.
#' @param alpha Type I error given to calculate confidence intervals, used when \code{hessian = T}.
#' @keywords d, p, q, r genlogis
#' 
#' @export
#' @examples 
#' 
#' ## Using generic parameter starting values
#' datas <- rgenlog(10000, 1.5,2,2, 0)
#' genlog_mle(c(.5,1.6, 1.5, 0),datas)
#' 
#' ## Select parameters starting values with genlog_slider
#' datas <- rgenlog(10000, 1.5,2,2, 0)
#' genlog_slider(datas, return_var = 'parameters') ## choose parameters
#' genlog_mle(parameters,datas)
#' 
#' @usage 
#' genlog_mle(parameters, data, hessian = F, alpha = 0.05)
#' 
#' @details 
#' Maximum likehood estimation of parameters for the distribution proposed in this package.\cr
#' This function is an application of \code{optim} as a particular case needed for this 
#' distribution using the method "\code{L-BFGS-B}".  \cr
#' For more information of the output check \code{help(optim)}.\cr
#' 
#' The used distribution for this package is given by: 
#' \deqn{f(x) = ((a + b*(1+p)*(abs(x-location)^p))*exp(-(x-location)*(a+b*(|x-location|^p)))) / ((exp(-(x-location)*(a + b* (|x-location|^p)))+1)^2)}
#' 
#' \code{help(dgenlog)} for parameters restrictions.\cr
#' Despite de function is defined from parameter values > 0, the lower bound for optimization
#' is set to 0.01 otherwise there is a convergence problem, specially for the hessian matrix.
#' 
#' @return 
#' 
#' Return a list of components as \code{optim} \(for more information, check this function\)
#'  with some extras: \cr
#'  
#' \code{par} The best set of parameters found. 
#' \cr
#' 
#' \code{value} The value of the loglikelihood function corresponding to \code{par}.
#' \cr
#' 
#' \code{counts} A two-element integer vector giving the number of calls to the likelihood 
#' function and \code{L-BFGS-B} respectively. This excludes those calls needed to 
#' compute the Hessian, and any calls to likelihood function to compute a 
#' finite-difference approximation to the gradient.\cr
#' 
#' \code{convergence} An integer code. 0 indicates successful completion, 1 indicates that 
#'the iteration limit \code{maxit} had been reached. For more errors \code{help(optim)}.
#' \cr
#' 
#' \code{message} A character string giving any additional information returned by the optimizer,
#'  or NULL.
#'  \cr
#'  
#' \code{hessian} A symmetric matrix giving an estimate of the (negative) Hessian at the solution found. 
#' Note that this is the Hessian of the unconstrained problem even if the box 
#' constraints are active.
#' \cr
#' 
#' \code{bounds} Return the best parameters found and the upper and lower bounds for the estimation.
#' \cr
#' 
#' @references 
#' RATHIE, P. N., SWAMEE, P. K. \emph{On a new invertible generalized logistic distribution
#' approximation to normal distribution}, Technical Research Report in Statistics, 07/2006,
#' Dept. of Statistics, Univ. of Brasilia, Brasilia, Brazil. 2006.
#'
#' Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C. \emph{A limited memory algorithm for bound
#' constrained optimization}, Technical Report NAM-08, Dept. of Elecrtrical Engineering and
#' Computer Science, Northwestern University, United States of America, mar. 1994.


genlog_mle <- function(parameters, data, hessian = F, alpha = 0.05){
  
  op <- optim(par=parameters, fn = genlogis.loglikehood, x=data,
              lower = c(0.01,0.01,0.01, -Inf), upper = c(Inf,Inf,Inf,Inf),
              method = 'L-BFGS-B', hessian = hessian) 
  
  if(hessian == T){
  fisher_info <- solve(op$hessian)
  var_matrix <- diag(fisher_info)
  var_matrix <- sqrt(var_matrix)
  upper <- op$par + qnorm(1-alpha) * var_matrix
  lower <- op$par - qnorm(1-alpha) * var_matrix
  interval <- data.frame(value=op$par, lower=lower, upper=upper)
  op$bounds <-  as.matrix(interval)
  }
  
  return(op)
  
}


