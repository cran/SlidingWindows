#' @title Descritive statistics with sliding windows.
#'
#' @description This function generates descriptive statistics of a univariate time series with sliding windows approach.
#'
#' @details This function include following measures: min, max, mean, median, standard deviation, skewness and kurtosis.
#'
#' @param y A vector containing univariate time series.
#'
#' @param w An integer value indicating the window size \eqn{w < length(y)}.
#'          If \eqn{w = length(y)}, will be computed the function will not slide.
#'
#' @param skewness A non-numeric value. See PerformanceAnalytics package.
#'
#' @param kurtosis A non-numeric value. See PerformanceAnalytics package.
#'
#' @return A list containing "w", "min","max","mean", "median", "standard deviation","skewness" and "kurtosis".
#'
#' @examples
#' y <- rnorm(100)
#' descritive.SlidingWindows(y, w=99, skewness="moment", kurtosis="moment")
#'
#' @importFrom stats median sd
#' @importFrom PerformanceAnalytics skewness kurtosis
#'
#' @references
#' Guedes, E.F. Modelo computacional para análise de movimentos e co-movimentos de mercados financeiros, Ph.D. thesis, Programa de Pós-graduação em Modelagem Computacional e Tecnologia Industrial. Centro Universitário Senai Cimatec, 2019.
#'
#' @export
    descritive.SlidingWindows <- function(y, w=99, skewness="moment", kurtosis="moment"){
      if(!(is.null(y) || is.numeric(y) || is.logical(y))){
    stop("Time series must be numeric")
      }
  N <- length(y)
  if(w > N){
    stop("The window needs to be smaller than the series length")
  }
  sw <- SlidingWindows(y,w)
  min_sw <- c()
  max_sw <- c()
  mean_sw <- c()
  median_sw <- c()
  sd_sw <- c()
  skewness_sw <- c()
  kurtosis_sw <- c()
  for(i in 1:nrow(sw)){
       min_sw[i] <-  min(sw[i,])
       max_sw[i] <-  max(sw[i,])
      mean_sw[i] <- mean(sw[i,])
        sd_sw[i] <- stats::sd(sw[i,])
    median_sw[i] <- stats::median(sw[i,])
  skewness_sw[i] <- PerformanceAnalytics::skewness(sw[i,], method=skewness)
  kurtosis_sw[i] <- PerformanceAnalytics::kurtosis(sw[i,], method=kurtosis)
  }
  return(list(w = w, min=min_sw, max=max_sw,mean=mean_sw, median=median_sw, sd=sd_sw, skewness= skewness_sw, kurtosis=kurtosis_sw))
}
