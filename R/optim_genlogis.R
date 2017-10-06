#' Optimization for a generalized logistic distribution
#'
#' Maximum likehood estimation of parameters for a generalized logistic distribution.
#' 
#' @param parameters Initial values for the parameters to be optimized over in the following order \code{c(a, b, p, mu)},
#'  \code{mu} can be omitted and will be set to 0.
#' @param data This is the the data to be utilized for the estimation.
#' @param hessian logical value that returns hessian, also returns the parameters estimation's confidence interval.
#' @param alpha Type I error given to calculate confidence intervals, used when \code{hessian = T}.
#' @keywords genlogis
#' 
#' @export
#' @examples 
#' 
#' ## Using generic parameter starting values
#' datas <- rgenlog(10000, 1.5,2,2, 0)
#' genlog_mle(c(.5,1.6, 1.5, 0),datas)
#' 
#' ## Select parameters starting values with genlog_slider
#' \donttest{
#' datas <- rgenlog(10000, 1.5,2,2, 0)
#' genlog_slider(datas, return_var = 'parameters') ## choose parameters
#' genlog_mle(parameters,datas)
#' }
#' 
#' @usage 
#' genlog_mle(parameters, data, hessian = F, alpha = 0.05)
#' 
#' @details 
#' Maximum likehood estimation of parameters for the distribution proposed in this package.\cr
#' This function is an application of \code{constrOptim} as a particular case needed for this 
#' distribution using the method "\code{BFGS}".  \cr
#' For more information of the output check \code{help(constrOptim)}.\cr
#' 
#' The used distribution for this package is given by: 
#' \deqn{f(x) = ((a + b*(1+p)*(abs(x-mu)^p))*exp(-(x-mu)*(a+b*(|x-mu|^p)))) / ((exp(-(x-mu)*(a + b* (|x-mu|^p)))+1)^2)}
#' 
#' \code{help(dgenlog)} for parameters restrictions.\cr
#' 
#' @return 
#' 
#' Return a list of components as \code{constrOptim} \(for more information, check this function\)
#'  with some extras: \cr
#'  
#' \code{par} The best set of parameters found. 
#' \cr
#' 
#' \code{value} The value of the loglikelihood function corresponding to \code{par}.
#' \cr
#' 
#' \code{counts} A two-element integer vector giving the number of calls to the likelihood 
#' function and \code{BFGS} respectively. This excludes those calls needed to 
#' compute the Hessian, and any calls to likelihood function to compute a 
#' finite-difference approximation to the gradient.\cr
#' 
#' \code{convergence} An integer code. 0 indicates successful completion, 1 indicates that 
#'the iteration limit \code{maxit} had been reached. For more errors \code{help(constrOptim)}.
#' \cr
#' 
#' \code{message} A character string giving any additional information returned by the optimizer,
#'  or NULL.
#'  \cr
#'  
#' \code{outer.iterations} gives the number of outer iterations (calls to \code{optim}).
#'  \cr
#'  
#'  \code{barrier.value} giving the value of the barrier function at the optimum.
#'  \cr
#'   
#' For \code{hessian = T} add:
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
#' Azzalini, A.. "A class of distributions which includes the normal ones". Scandinavian Journal of Statistics, 1985.
#' Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C. \emph{A limited memory algorithm for bound
#' constrained optimization}, Technical Report NAM-08, Dept. of Elecrtrical Engineering and
#' Computer Science, Northwestern University, United States of America, mar. 1994.


genlog_mle <- function(parameters, data, hessian = F, alpha = 0.05){
  
  fn <- function(p) genlogis.loglikehood(param = p, x = data)
  grr <- function(k) genlogis.likeli.gr(param = k, x = data)
  ui <- rbind(c(1, 0, 0, 0), c(0,1,0,0), c(0,0,1,0))
  ci <- c(0,0,0)    
  op <- stats::constrOptim(theta = parameters, f = fn, grad = grr,
                    ui = ui, ci = ci, method = "BFGS",
                    hessian = hessian)  
  
  if(hessian == T){
    fisher_info <- solve(op$hessian)
    var_matrix <- diag(fisher_info)
    var_matrix <- sqrt(var_matrix)
    upper <- op$par + stats::qnorm(1-alpha) * var_matrix
    lower <- op$par - stats::qnorm(1-alpha) * var_matrix
    interval <- data.frame(value=op$par, lower=lower, upper=upper)
    op$bounds <-  as.matrix(interval)
  }
  
  return(op)
  
}

#' Optimization for a generalized logistic distribution with skewness
#'
#' Maximum likehood estimation of parameters for a generalized logistic distribution with skewness.
#' 
#' @param parameters Initial values for the parameters to be optimized over in the following order \code{c(a, b, p, mu, skew)}.
#' @param data This is the the data to be utilized for the estimation.
#' @param hessian logical value that returns hessian, also returns the parameters estimation's confidence interval.
#' @param alpha Type I error given to calculate confidence intervals, used when \code{hessian = T}.
#' @keywords genlogis
#' 
#' @export
#' @examples 
#' 
#' ## Using generic parameter starting values
#' datas <- rgenlog_as(10000, 0.3,0.9,1.5, 0, 0.9)
#' genlog_mle_as(c(0.3,0.9,1.5, 0, 0.9),datas)
#' 
#' ## Select parameters starting values with genlog_slider
#' \donttest{
#' datas <- rgenlog(10000, 1.5,2,2, 0)
#' genlog_slider(datas, return_var = 'parameters', skew = T) ## choose parameters
#' genlog_mle_as(parameters,datas)
#' }
#' 
#' @usage 
#' genlog_mle_as(parameters, data, hessian = F, alpha = 0.05)
#' 
#' @details 
#' Maximum likehood estimation of parameters for the distribution proposed in this package.\cr
#' This function is an application of \code{constrOptim} as a particular case needed for this 
#' distribution using the method "\code{BFGS}".  \cr
#' For more information of the output check \code{help(constrOptim)}.\cr
#' 
#' The used distribution for this package is given by: 
#' \deqn{f(x) = 2*((a + b*(1+p)*(abs(x-mu)^p))*exp(-(x-mu)*(a+b*(abs(x-mu)^p))))/ 
#'    ((exp(-(x-mu)*(a + b* (abs(x-mu)^p)))+1)^2) *
#'    ((exp(-(skew*(x-mu))*(a+b*(abs(skew*(x-mu))^p)))+1)^(-1)) }
#' 
#' \code{help(dgenlog_as)} for parameters restrictions.\cr
#' 
#' @return 
#' 
#' Return a list of components as \code{constrOptim} \(for more information, check this function\)
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
#'the iteration limit \code{maxit} had been reached. For more errors \code{help(constrOptim)}.
#' \cr
#'  
#' \code{message} A character string giving any additional information returned by the optimizer,
#'  or NULL.
#'  \cr
#'  
#' \code{outer.iterations} gives the number of outer iterations (calls to \code{optim}).
#'  \cr
#'  
#'  \code{barrier.value} giving the value of the barrier function at the optimum.
#'  \cr
#'  
#' For \code{hessian = T} add:
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
#' 
#' Azzalini, A. \emph{A class of distributions which includes the normal ones}. Scandinavian Journal of Statistics, 1985.


genlog_mle_as <- function(parameters, data, hessian = F, alpha = 0.05){

  fn <- function(p) genlogis.as.loglikehood(param = p, x = data)
  grr <- function(k) genlogis.as.likeli.gr(param = k, x = data)
  ui <- rbind(c(1, 0, 0, 0, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0), 
              c(0, 0, 0, 0, 1), c(0, 0, 0, 0, -1))
  ci <- c(0, 0, 0, -1, -1)    
  op <- stats::constrOptim(theta = parameters, f = fn, grad = grr,
                           ui = ui, ci = ci, method = "BFGS",
                           hessian = hessian)  
  
  if(hessian == T){
    fisher_info <- solve(op$hessian)
    var_matrix <- diag(fisher_info)
    var_matrix <- sqrt(var_matrix)
    upper <- op$par + stats::qnorm(1-alpha) * var_matrix
    lower <- op$par - stats::qnorm(1-alpha) * var_matrix
    interval <- data.frame(value=op$par, lower=lower, upper=upper)
    op$bounds <-  as.matrix(interval)
  }
  
  return(op)
  
}