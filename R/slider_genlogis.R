#' Slider for generalized logistic
#'
#' Make a generalized logistic distribution slider to compare histogram with theoretical distribution
#' @param data vector of data to compare.
#' @param return_var a char string to name where parameters are assigned
#' @param loc_range a number to setup the minimum and maximum range value of the mu parameter
#' @param skew logical, if \code{TRUE}, a model with skewness should be used.
#' @keywords genlogis
#' 
#' @export
#' @examples 
#' datas <- rgenlog(1000)
#' genlog_slider(datas, return_var = 'parameters')
#' 
#' @usage 
#' genlog_slider(data, return_var = NULL, loc_range = 10, skew = F)
#' 
#' @return 
#' The function plots a interactive graphic in RStudio Viewer panel. \cr
#' Also, the parameters \code{a, b, p} and \code{mu} can be returned to 
#' \code{return_var} if asked in the graphic.
#' 
#' @details 
#' 
#' There is a small gear in the top left of the graphic where you can slide the parameters \code{@param a,b,p,mu}.
#' The used distribution for this package is given by: \deqn{f(x) = ((a + b*(1+p)*(abs(x-mu)^p))*exp(-(x-mu)*(a+b*(|x-mu|^p)))) / ((exp(-(x-mu)*(a + b* (|x-mu|^p)))+1)^2)}
#' If the density function is not printed it is not defined for these parameters. \cr 
#' 
#' \code{help(dgenlog)} for parameters restrictions.\cr 
#' 
#' This function requires RStudio to run.
#' \cr 
#'
#' @references 
#' RATHIE, P. N., SWAMEE, P. K. \emph{On a new invertible generalized logistic distribution
#' approximation to normal distribution}, Technical Research Report in Statistics, 07/2006,
#' Dept. of Statistics, Univ. of Brasilia, Brasilia, Brazil. 2006.

genlog_slider <- function(data, return_var = NULL, loc_range = 10, skew = F){
  
  if (!manipulate::isAvailable())
    stop("The genlog_slider function must be run from within RStudio")

  data <- as.data.frame(data)

if(skew == F){
    manipulate::manipulate({
    manipulate::manipulatorSetState("parameters", c(par_a, par_b, par_p, par_mu))
    
    if(returnval){
      if(!is.null(return_var)){
        assign(paste(return_var) ,manipulate::manipulatorGetState("parameters"), envir = .GlobalEnv)
        message(paste0('Parameters returned to the object \'', return_var, '\''))}
      else{
        message('No object set to receive the parameters')
      }
    }
    
    ggplot2::ggplot(data = data, ggplot2::aes(x = data)) +
      ggplot2::theme_bw()+
      ggplot2::geom_histogram(binwidth = binw,
                              colour = 'black', ggplot2::aes(y = ..density.., fill =..count..))+
      ggplot2::scale_fill_gradient("Count", low="#DCDCDC", high="#7C7C7C")+
      ggplot2::stat_function(fun= dgenlog,
                             args = c(a = par_a, b = par_b, p = par_p, mu = par_mu),
                             color="red", size = 1) +
      ggplot2::labs(y = 'Density', x = 'X', 
                    title = 'Theoretical density vs Observed histogram') +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = 'bold'))
    
  },
  
  par_a = manipulate::slider(0,10, step = 0.01, initial = sqrt(2/pi), label = 'Parameter a'),
  par_b = manipulate::slider(0,10, step = 0.01, initial = 0.5, label = 'Parameter b'),
  par_p = manipulate::slider(0,10, step = 0.01, initial = 2, label = 'Parameter p'),
  par_mu = manipulate::slider(-loc_range,loc_range, step = 0.01, initial = 0, label = 'mu parameter'),
  binw = manipulate::slider(0.1,10, step = 0.1, initial = 0.1, label = 'Binwidth'),
  
  returnval = manipulate::button("Return parameters to variable")
  )
}
  if(skew == T){
    manipulate::manipulate({
      manipulate::manipulatorSetState("parameters", c(par_a, par_b, par_p, par_mu, par_skew))
      
      if(returnval){
        if(!is.null(return_var)){
          assign(paste(return_var) ,manipulate::manipulatorGetState("parameters"), envir = .GlobalEnv)
          message(paste0('Parameters returned to the object \'', return_var, '\''))}
        else{
          message('No object set to receive the parameters')
        }
      }
      
      ggplot2::ggplot(data = data, ggplot2::aes(x = data)) +
        ggplot2::theme_bw()+
        ggplot2::geom_histogram(binwidth = binw,
                                colour = 'black', ggplot2::aes(y = ..density.., fill =..count..))+
        ggplot2::scale_fill_gradient("Count", low="#DCDCDC", high="#7C7C7C")+
        ggplot2::stat_function(fun= dgenlog_as,
                               args = c(a = par_a, b = par_b, p = par_p, mu = par_mu, skew = par_skew),
                               color="red", size = 1) +
        ggplot2::labs(y = 'Density', x = 'X', 
                      title = 'Theoretical density vs Observed histogram') +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = 'bold'))
      
    },
    
    par_a = manipulate::slider(0,10, step = 0.01, initial = sqrt(2/pi), label = 'Parameter a'),
    par_b = manipulate::slider(0,10, step = 0.01, initial = 0.5, label = 'Parameter b'),
    par_p = manipulate::slider(0,10, step = 0.01, initial = 2, label = 'Parameter p'),
    par_mu = manipulate::slider(-loc_range,loc_range, step = 0.01, initial = 0, label = 'mu parameter'),
    par_skew = manipulate::slider(-1,1, step = 0.01, initial = 0, label = 'Skewness'),
    binw = manipulate::slider(0.1,10, step = 0.1, initial = 0.1, label = 'Binwidth'),
    
    returnval = manipulate::button("Return parameters to variable")
    )
  }  

}
