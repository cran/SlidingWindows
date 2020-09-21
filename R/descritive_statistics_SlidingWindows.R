#' @title Descritive statistics with sliding windows.
#'
#' @description This function generates descritive statistics of
#' a univariate time serieswith sliding windows approach.
#'
#' @details This function include following measures:
#'
#' mean. median, standard deviation, skewness and kurtosis.
#'
#' @param y A vector contaning univariate time series.
#' @param w An integer value indicating the size of the window \eqn{w < length(y)}.
#'          If \eqn{w = length(y)}, will be computed the function will not slide.
#' @param skewness A non-numeric value. See PerformanceAnalytics package.
#' @param kurtosis A non-numeric value. See PerformanceAnalytics package.
#'
#' @return A list contaning "w", "mean", "median", "standard deviation","skewness" and "kurtosis".
#'
#' @examples
#' y <- rnorm(1000)
#' descritive_statsistics_SlidingWindows(rnorm(100), 99, skewness="moment", kurtosis="moment")
#'
#' @importFrom stats median sd
#' @importFrom PerformanceAnalytics skewness kurtosis
#'
#' @references
#' Guedes, E.F. Modelo computacional para análise de movimentos e co-movimentos de mercados financeiros, Ph.D. thesis, Programa de Pós-graduação em Modelagem Computacional e Tecnologia Industrial. Centro Universitário Senai Cimatec, 2019.
#'
#' @export
    descritive_statsistics_SlidingWindows <- function(y,w, skewness=c("moment","sample","fisher"), kurtosis=c("moment","sample","fisher", "excess", "sample_excess")){
      if(!(is.null(y) || is.numeric(y) || is.logical(y))){
    stop("Time series must be numeric")
   }
  N <- length(y)
  if(w > N){
    stop("The window needs to be smaller than the series length")
  }
  sw <- SlidingWindows(y,w)
  mean_sw <- c()
  median_sw <- c()
  sd_sw <- c()
  skewness_sw <- c()
  kurtosis_sw <- c()
  for(i in 1:nrow(sw)){
      mean_sw[i] <- mean(sw[i,])
    median_sw[i] <- stats::median(sw[i,])
        sd_sw[i] <- stats::sd(sw[i,])
  skewness_sw[i] <- PerformanceAnalytics::skewness(sw[i,], method=skewness)
  kurtosis_sw[i] <- PerformanceAnalytics::kurtosis(sw[i,], method=kurtosis)
  }
  return(list(w = w, mean_SlidingWindows=mean_sw, median_SlidingWindows=median_sw, sd_SlidingWindows=sd_sw, skewness_SlidingWindows= skewness_sw, kurtosis_SlidingWindows=kurtosis_sw))
}
path.expand("~/teste/SlidingWindows")
