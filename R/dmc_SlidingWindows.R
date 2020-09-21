#' @title Detrended Multiple Cross-Correlation Coefficient with sliding windows.
#
#' @description This function generates Detrended Cross-Correlation Coefficient of
#' three time series with sliding windows approach.
#'
#' @details This function include following measures:
#'
#'     w, timescale, dmc, rhodcca_yx1, rhodcca_yx2, rhodcca_x1x2
#'
#' @param x1 A vector contaning univariate time series.
#' @param x2 A vector contaning univariate time series.
#' @param y A vector contaning univariate time series.
#' @param w An integer value indicating the size of the window \eqn{w < length(y)}.
#'          If \eqn{w = length(y)}, will be computed the function will not slide.
#' @param k An integer value indicating the boundary of the division \eqn{(N/k)}.
#'          The smallest value of \eqn{k} is \eqn{4}.
#' @param nu An integer value. See the DCCA package.
#'
#' @return A list contaning "w", "dmc", "rhodcca_yx1", "rhodcca_yx2", "rhodcca_x1x2".
#'
#' @examples
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' y <- rnorm(100)
#' dmc_SlidingWindows(x1,x2,y,w=99,k=10,nu=0)
#'
#' @references
#' ZEBENDE, G.F.; SILVA-FILHO, A.M. Detrended Multiple Cross-Correlation Coefficient. PHYSICA A, v.510, p.91-97, 2018. doi = "https://doi.org/10.1016/j.physa.2018.06.119".
#'
#' @importFrom DCCA rhodcca
#'
#' @export
dmc_SlidingWindows <- function(x1,x2,y,w,k,nu){
  N1 <- length(x1)
  N2 <- length(x2)
  N3 <- length(y)
   m <- 4:round(w/k,0)
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
        yx1 <- DCCA::rhodcca(y, x1, m=m, nu=nu)$rhodcca
        yx2 <- DCCA::rhodcca(y, x2, m=m, nu=nu)$rhodcca
       x1x2 <- DCCA::rhodcca(x1, x2, m=m, nu=nu)$rhodcca
        dmc <- (yx1^2 + yx2^2 - (2*yx1*yx2*x1x2))/(1-x1x2^2)
      return(list(w = w, timescale=m, dmc=dmc, rhodcca_yx1=yx1, rhodcca_yx2=yx2, rhodcca_x1x2=x1x2))
  }
  if(w < N1){
  x1_sw <- SlidingWindows(x1,w)
  x2_sw <- SlidingWindows(x2,w)
   y_sw <- SlidingWindows(y,w)
    yx1 <- matrix(data = NA, nrow = nrow(x1_sw), ncol = length(m), byrow = TRUE)
    yx2 <- matrix(data = NA, nrow = nrow(x1_sw), ncol = length(m), byrow = TRUE)
   x1x2 <- matrix(data = NA, nrow = nrow(x1_sw), ncol = length(m), byrow = TRUE)
    dmc <- matrix(data = NA, nrow = nrow(x1_sw), ncol = length(m), byrow = TRUE)
        for(i in 1:nrow(x1_sw)){
      		   for(j in 1:length(m)){
           yx1[i,j]  <- DCCA::rhodcca(y_sw[i,],  x1_sw[i,], m=m[j], nu=nu)$rhodcca
           yx2[i,j]  <- DCCA::rhodcca(y_sw[i,],  x2_sw[i,], m=m[j], nu=nu)$rhodcca
           x1x2[i,j] <- DCCA::rhodcca(x1_sw[i,], x2_sw[i,], m=m[j], nu=nu)$rhodcca
           dmc[i,j]  <- (yx1[i,j]^2 + yx2[i,j]^2-(2*yx1[i,j]*yx2[i,j]*x1x2[i,j]))/(1-x1x2[i,j]^2)
		   }
		  }
	  return(list(w = w, timescale=m, dmc=dmc, rhodcca_yx1 = yx1, rhodcca_yx2 = yx2, rhodcca_x1x2 = x1x2))
 }
}
path.expand("~/teste/SlidingWindows")

