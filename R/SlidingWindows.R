#' @title Sliding Windows.
#
#' @description This function generates sliding windows approach of a time series.
#'
#' @details This function return the matrix with time series sliding windows.
#'
#' @param y A vector containing univariate time series.
#'
#' @param w An integer value indicating the window size \eqn{w < length(y)}.
#'          If \eqn{w = length(y)}, will be computed the function will not slide.
#'
#' @return A list containing "w", "SlidingWindows".
#'
#' @examples
#' y <- rnorm(100)
#' SlidingWindows(y,w=99)
#'
#' @references
#' Guedes, E.F. Modelo computacional para análise de movimentos e co-movimentos de mercados financeiros, Ph.D. thesis, Programa de Pós-graduação em Modelagem Computacional e Tecnologia Industrial. Centro Universitário Senai Cimatec, 2019.
#'
#' @export
SlidingWindows <- function(y,w=99){
 if(!(is.null(y) || is.numeric(y) || is.logical(y))){
    stop("Time series must be numeric")
	}
  N <- length(y)
  slide_win=c()
  for(i in 1:(N-w+1)){
    slide_win <- rbind(slide_win, y[i:(i+w-1)])
  }
  return(slide_win)
}
