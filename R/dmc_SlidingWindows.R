#' @title Multiple detrended cross-correlation coefficient with sliding windows.
#
#' @description This function generates DMC Coefficient of three time series with sliding windows approach.
#'
#' @details This function include following measures: w, timescale, dmc and cross-correlation between: yx1, yx2, x1x2
#'
#' @param x1 A vector containing univariate time series.
#'
#' @param x2 A vector containing univariate time series.
#'
#' @param y A vector containing univariate time series.
#'
#' @param w An integer value indicating the window size \eqn{w < length(y)}.
#'          If \eqn{w = length(y)}, will be computed the function will not slide.
#'
#' @param k An integer value indicating the boundary of the division \eqn{(N/k)}.
#'          The smallest value of \eqn{k} is \eqn{4}.
#'
#'@param method A character string indicating which correlation coefficient is to be used. If method = "rhodcca" the dmc coefficient is generated from the DCCA coefficient. If method = "dmca", the dmc coefficient is generated from the DMCA coefficient.
#'
#' @param nu An integer value. See the DCCA package.
#'
#' @return A list containing "w", "dmc", "yx1", "yx2", "x1x2".
#'
#' @examples
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <- rnorm(100)
#' dmc_SlidingWindows(x1,x2,y,w=99,k=10,nu=0, method="rhodcca")
#' dmc_SlidingWindows(x1,x2,y,w=99,k=10,nu=0, method="dmca")
#'
#' @references
#' ZEBENDE, G.; SILVA-FILHO, A.M. Detrended multiple cross-correlation coefficient, Physica A 510, 91-97, 2018.
#'
#' GUEDES, E.F.; ZEBENDE, G.F. DCCA cross-correlation coefficient with sliding windows approach. PHYSICA A, v.527, p.121286, 2019.
#'
#' ZEBENDE, G.F. DCCA cross-correlation coefficient: Quantifying level of cross-correlation, Physica A, v. 390, n. 4, p. 614-618, 2011.
#'
#' KRISTOUFEK, L. Detrending moving-average cross-correlation coefficient: Measuring cross-correlations between non-stationary series. PHYSICA A, v.406, p.169-175, 2014.
#'
#' @importFrom DCCA rhodcca
#' @importFrom stats filter
#'
#' @export
dmc_SlidingWindows <- function(x1,x2,y,w,k,method,nu){

   N1 <- length(x1)
   N2 <- length(x2)
   N3 <- length(y)
    n <- 4:round(w/k,0)

   dmca <- function(x,y,n){
     xx <- cumsum(x)
     yy <- cumsum(y)

     mm <- c(rep(1,n))/n
     mm_x <- stats::filter(xx,mm)
     mm_y <- stats::filter(yy,mm)

     F2_xy <- mean((xx-mm_x)[(1+floor(n/2)):(length(xx)-floor(n/2))]*(yy-mm_y)[(1+floor(n/2)):(length(yy)-floor(n/2))])
     F2_xx <- mean((xx-mm_x)[(1+floor(n/2)):(length(xx)-floor(n/2))]*(xx-mm_x)[(1+floor(n/2)):(length(xx)-floor(n/2))])
     F2_yy <- mean((yy-mm_y)[(1+floor(n/2)):(length(yy)-floor(n/2))]*(yy-mm_y)[(1+floor(n/2)):(length(yy)-floor(n/2))])

     rho <- F2_xy/sqrt(F2_xx*F2_yy)
     return(rho)
   }

 if(!(is.null(y) || is.numeric(y) || is.logical(y))){
    stop("Time series must be numeric")
	}
 if(!(is.null(x1) || is.numeric(x1) || is.logical(x1))){
    stop("Time series must be numeric")
	}
 if(!(is.null(x2) || is.numeric(x2) || is.logical(x2))){
    stop("Time series must be numeric")
	}
  if(N1 != N2){
    stop("Time series have different lengths")
  }
  if(N1 != N3){
    stop("Time series have different lengths")
  }
  if(N2 != N3){
    stop("Time series have different lengths")
  }
  if(w > N1){
    stop("The window needs to be smaller than the series length")
  }

  if(w == N1){
    if(method =='rhodcca'){
        yx1 <- DCCA::rhodcca(y, x1, m=n, nu=nu)$rhodcca
        yx2 <- DCCA::rhodcca(y, x2, m=n, nu=nu)$rhodcca
       x1x2 <- DCCA::rhodcca(x1, x2, m=n, nu=nu)$rhodcca
        dmc <- (yx1^2 + yx2^2 - (2*yx1*yx2*x1x2))/(1-x1x2^2)
      return(list(w = w, timescale=n, dmc=dmc, yx1=yx1, yx2=yx2, x1x2=x1x2))
      }
    if(method =='dmca'){
        yx1 <- dmca(y, x1, n)
        yx2 <- dmca(y, x2, n)
       x1x2 <- dmca(x1, x2, n)
        dmc <- (yx1^2 + yx2^2 - (2*yx1*yx2*x1x2))/(1-x1x2^2)
      return(list(w = w, timescale=n, dmc=dmc, yx1=yx1, yx2=yx2, x1x2=x1x2))
      }
  }

   if(w < N1){

     x1_sw <- SlidingWindows(x1,w)
     x2_sw <- SlidingWindows(x2,w)
     y_sw <- SlidingWindows(y,w)
     yx1 <- matrix(data = NA, nrow = nrow(x1_sw), ncol = length(n), byrow = TRUE)
     yx2 <- matrix(data = NA, nrow = nrow(x1_sw), ncol = length(n), byrow = TRUE)
     x1x2 <- matrix(data = NA, nrow = nrow(x1_sw), ncol = length(n), byrow = TRUE)
     dmc <- matrix(data = NA, nrow = nrow(x1_sw), ncol = length(n), byrow = TRUE)

     if(method =='rhodcca'){
        for(i in 1:nrow(x1_sw)){
      		   for(j in 1:length(n)){
           yx1[i,j]  <- DCCA::rhodcca(y_sw[i,],  x1_sw[i,], m=n[j], nu=nu)$rhodcca
           yx2[i,j]  <- DCCA::rhodcca(y_sw[i,],  x2_sw[i,], m=n[j], nu=nu)$rhodcca
           x1x2[i,j] <- DCCA::rhodcca(x1_sw[i,], x2_sw[i,], m=n[j], nu=nu)$rhodcca
           dmc[i,j]  <- (yx1[i,j]^2 + yx2[i,j]^2-(2*yx1[i,j]*yx2[i,j]*x1x2[i,j]))/(1-x1x2[i,j]^2)
		   }
		  }
	  return(list(w = w, timescale=n, dmc=dmc, yx1 = yx1, yx2 = yx2, x1x2 = x1x2))
 }

 if(method =='dmca'){
        for(i in 1:nrow(x1_sw)){
      		   for(j in 1:length(n)){
           yx1[i,j]  <- dmca(y_sw[i,],  x1_sw[i,], n[j])
           yx2[i,j]  <- dmca(y_sw[i,],  x2_sw[i,], n[j])
           x1x2[i,j] <- dmca(x1_sw[i,], x2_sw[i,], n[j])
           dmc[i,j]  <- (yx1[i,j]^2 + yx2[i,j]^2-(2*yx1[i,j]*yx2[i,j]*x1x2[i,j]))/(1-x1x2[i,j]^2)
		   }
		  }
	  return(list(w = w, timescale=n, dmc=dmc, yx1 = yx1, yx2 = yx2, x1x2 = x1x2))
   }
 }
}
