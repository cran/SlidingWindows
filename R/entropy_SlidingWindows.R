#' @title Approximate entropy with sliding windows.
#'
#' @description This function computes approximate entropy of a univariate time series  with sliding windows approach.
#'
#' @details This function return the list with time series sliding windows.
#'
#' @param y A vector containing univariate time series.
#'
#' @param w An integer value indicating the window size \eqn{w < length(y)}.
#'          If \eqn{w = length(y)}, will be computed the function will not slide.
#'
#' @param k An integer value indicating the boundary of the division \eqn{(N/k)}.
#'          The smallest value of \eqn{k} is \eqn{4}.
#'
#' @param dim The dimension of given time series. See TSEntropies package.
#'
#' @param r The radius of searched areas. See TSEntropies package.
#'
#' @param lag The downsampling. See TSEntropies package.
#'
#' @return A list contaning "w", "ApEn", "FastApEn".
#'
#' @examples
#' y <- rnorm(100)
#' entropy.SlidingWindows(y, w=99, k=4, dim=2, r=.2,lag=1)
#'
#' @references
#' Pincus, S.M. (1991). Approximate entropy as a measure of system complexity. Proc. Natl. Acad. Sci. USA, Vol. 88, pp. 2297–2301. doi="doi.org/10.1073/pnas.88.6.2297".
#'
#' @importFrom TSEntropies ApEn FastApEn
#' @importFrom stats sd
#'
#' @export
entropy.SlidingWindows <- function(y,w=99,k=4,dim=2,r=0.5,lag=1){
  if(!(is.null(y) || is.numeric(y) || is.logical(y))){
   stop("Time series must be numeric")
   }

  Ny <- length(y)

  if(w > Ny){
    stop("The window needs to be smaller than the series length")
  }

  if(w == Ny){

    m <- 4:round(w/k,0)

    r3 <- c()
    r4 <- c()
    ApEn <- c()
    FastApEn <- c()

    for(i in 1:length(m)){
      ## Divide in to overlapping boxes of size m
      ny <- SlidingWindows(y, m[i])

      for(j in 1:nrow(ny)){

        r3[j] <-     TSEntropies::ApEn(ny[j,], dim = dim, lag = lag, r = r*stats::sd(ny[j,]))
        r4[j] <- TSEntropies::FastApEn(ny[j,], dim = dim, lag = lag, r = r*stats::sd(ny[j,]))
      }

      ApEn[i] <- mean(r3)
      FastApEn[i] <- mean(r4)

      ny <- NULL
    }

    return(list(window = w,
                timescale = m,
                ApEn=ApEn,
                FastApEn=FastApEn))
  }


  if(w < Ny){
     sw <- SlidingWindows(y,w)

     m <- 4:round(w/k,0)

     r3 <- c()
     r4 <- c()

     ApEn <- matrix(data = NA, nrow = nrow(sw), ncol = length(m), byrow = TRUE)
 FastApEn <- matrix(data = NA, nrow = nrow(sw), ncol = length(m), byrow = TRUE)

for(h in 1:nrow(sw)){
  for(i in 1:length(m)){
    ## Divide in to overlapping boxes of size n
    ny <- SlidingWindows(sw[h,], m[i])

    for(j in 1:nrow(ny)){

      r3[j] <-       TSEntropies::ApEn(ny[j,], dim = dim, lag = lag, r = r*sd(ny[j,]))
      r4[j] <-   TSEntropies::FastApEn(ny[j,], dim = dim, lag = lag, r = r*sd(ny[j,]))
    }

        ApEn[h, i] <- mean(r3)
    FastApEn[h, i] <- mean(r4)

    ny <- NULL

   }
  }
  return(list(window = w,
              timescale = m,
              ApEn=ApEn,
              FastApEn=FastApEn))
  }
}

