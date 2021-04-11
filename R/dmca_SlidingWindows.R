#' @title DMCA coefficient with sliding windows.
#
#' @description This function generates Detrending moving-average cross-correlation coefficient of two time series with sliding windows approach.
#'
#' @details This function include following measures: w, timescale, dmca
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
#' @return A list containing "w", "timescale", "dmca".
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' dmca.SlidingWindows(x,y,w=99,k=10)
#'
#' @references
#' KRISTOUFEK, L. Detrending moving-average cross-correlation coefficient: Measuring cross-correlations between non-stationary series. PHYSICA A, v.406, p.169-175, 2014. doi="doi.org/10.1016/j.physa.2014.03.015".
#'
#' @importFrom stats filter
#'
#' @export
dmca.SlidingWindows <- function(x,y,w=98,k=10){
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

  dmca <- function(x,y,m){
    xx <- cumsum(x)
    yy <- cumsum(y)

    mm <- c(rep(1,m))/m
    mm_x <- stats::filter(xx,mm)
    mm_y <- stats::filter(yy,mm)

    F2_xy <- mean((xx-mm_x)[(1+floor(m/2)):(length(xx)-floor(m/2))]*(yy-mm_y)[(1+floor(m/2)):(length(yy)-floor(m/2))])
    F2_xx <- mean((xx-mm_x)[(1+floor(m/2)):(length(xx)-floor(m/2))]*(xx-mm_x)[(1+floor(m/2)):(length(xx)-floor(m/2))])
    F2_yy <- mean((yy-mm_y)[(1+floor(m/2)):(length(yy)-floor(m/2))]*(yy-mm_y)[(1+floor(m/2)):(length(yy)-floor(m/2))])

    rho <- F2_xy/sqrt(F2_xx*F2_yy)
    return(rho)
  }


  if(w == Nx){
    yx <- c()
        for(i in 1:length(m)){
              yx[i] <- dmca(x,y,m[i])
        }

       return(list(w = w, timescale=m, dmca=yx))
  }

  if(w < Nx){

  x_sw <- SlidingWindows(x,w)
  y_sw <- SlidingWindows(y,w)

  yx <- matrix(data = NA, nrow = nrow(x_sw), ncol = length(m), byrow = TRUE)

  for(i in 1:nrow(x_sw)){
      for(j in 1:length(m)){
        yx[i,j] <- dmca(y_sw[i,], x_sw[i,], m[j])
      }
	}
  return(list(w = w, timescale=m, dmca=yx))
  }
}
