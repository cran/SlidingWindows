#' @title Detrended Cross-Correlation Coefficient with sliding windows.
#
#' @description This function generates Detrended Cross-Correlation Coefficient of
#' two time series with sliding windows approach.
#'
#' @details This function include following measures:
#'
#'     w, timescale, rhodcca
#'
#' @param x A vector containing univariate time series.
#'
#' @param y A vector containing univariate time series.
#'
#' @param w An integer value indicating the window size \eqn{w < length(y)}.
#'          If \eqn{w = length(y)}, will be computed the function will not slide.
#'
#' @param k An integer value indicating the boundary of the division \eqn{(N/k)}.
#'          The smallest value of \eqn{k} is \eqn{4}.
#'
#' @param nu An integer value. See DCCA package.
#'
#' @return A list containing "w", "timescale", "rhodcca".
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' rhodcca_SlidingWindows(x,y,w=99,k=10,nu=0)
#'
#' @references
#' GUEDES, E.F.; ZEBENDE, G.F. DCCA cross-correlation coefficient with sliding windows approach. PHYSICA A, v.527, p.121286, 2019.
#'
#' ZEBENDE, G.F. DCCA cross-correlation coefficient: Quantifying level of cross-correlation, Physica A, v. 390, n. 4, p. 614-618, 2011.
#'
#' @importFrom DCCA rhodcca
#'
#' @export
rhodcca_SlidingWindows <- function(x,y,w,k,nu){
 if(!(is.null(y) || is.numeric(y) || is.logical(y))){
    stop("Time series must be numeric")
	}
 if(!(is.null(x) || is.numeric(x) || is.logical(x))){
    stop("Time series must be numeric")
	}

  Ny <- length(y)
  Nx <- length(x)
  m <- 4:round(w/k,0)

  if(Nx != Ny){
    stop("Time series have different lengths")
  }

    if(w > Nx){
    stop("The window needs to be smaller than the series length")
  }

  if(w == Nx){
        yx <- DCCA::rhodcca(y, x, m=m, nu=nu)$rhodcca
      return(list(w = w, timescale=m, rhodcca=yx))
  }

  if(w < Nx){
  x_sw <- SlidingWindows(x,w)
  y_sw <- SlidingWindows(y,w)
  rho <- matrix(data = NA, nrow = nrow(x_sw), ncol = length(m), byrow = TRUE)
  for(i in 1:nrow(x_sw)){
      for(j in 1:length(m)){
    rho[i,j] <- DCCA::rhodcca(x_sw[i,], y_sw[i,], m=m[j], nu=nu)$rhodcca
      }
	}
  return(list(w = w, timescale=m, rhodcca=rho))
  }
}
