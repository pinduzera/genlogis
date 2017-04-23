#' Slider for generalized logistic
#'
#' Make a generalized logistic distribution slider to compare histogram with theoretical distribution
#' @param data vector of data to compare.
#' 
#' @keywords slieder, genlogis
#' @export
#' @examples 
#' genlog_slider(rgenlog(1000))
#' 
#' @usage 
#' genlog_slider(data)
#' 
#' @details 
#' 
#' There is a small gear in the top left of the graphic where you can slide the parameters \code{@param a,b,p,location}.
#' The used distribution for this package is given by: \deqn{f(x) = ((a + b*(1+p)*(abs(x-location)^p))*exp(-(x-location)*(a+b*(|x-location|^p)))) / ((exp(-(x-location)*(a + b* (|x-location|^p)))+1)^2)}
#' If the density function is not printed it is not defined for these parameters. \cr 
#' 
#' \code{help(dgenlog)} for parameters restrictions.\cr 
#' 
#' This function requires RStudio to run.\cr 

genlog_slider <- function(data){

data <- as.data.frame(data)
    
manipulate::manipulate(
  ggplot2::ggplot(data = data, ggplot2::aes(x = data)) +
    ggplot2::theme_bw()+
    ggplot2::geom_histogram(binwidth = 0.1,
                   colour = 'black', ggplot2::aes(y = ..density.., fill =..count..))+
    ggplot2::scale_fill_gradient("Count", low="#DCDCDC", high="#7C7C7C")+
    ggplot2::stat_function(fun= genlogis::dgenlog,
                  args = c(a = par_a, b = par_b, p = par_p, location = par_location),
                  color="red", size = 1) +
    ggplot2::labs(y = 'Density', x = 'X', 
         title = 'Theoretical density vs Observed histogram') +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = 'bold')),

  par_a = manipulate::slider(0,10, step = 0.01, initial = sqrt(2/pi)),
  par_b = manipulate::slider(0,10, step = 0.01, initial = 0.5),
  par_p = manipulate::slider(0,10, step = 0.01, initial = 2),
  par_location = manipulate::slider(-10,10, step = 0.01, initial = 0)
)
}
