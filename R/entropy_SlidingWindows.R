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
#' @param dim The dimension of given time series. See TSEntropies r package.
#'
#' @param r The radius of searched areas. See TSEntropies package.
#'
#' @param lag The downsampling. See TSEntropies r package.
#'
#' @return A list contaning "w", "ApEn", "FastApEn", "SampEn".
#'
#' @examples
#' y <- rnorm(100)
#' entropy_SlidingWindows(y, w=100, k=10, dim=2,r=.5,lag=1)
#'
#' @references
#' Pincus, S.M. (1991). Approximate entropy as a measure of system complexity. Proc. Natl. Acad. Sci. USA, Vol. 88, pp. 2297–2301.
#'
#' @importFrom TSEntropies ApEn FastApEn SampEn
#' @importFrom stats sd
#'
#' @export
entropy_SlidingWindows <- function(y,w,k,dim,r,lag){
  if(!(is.null(y) || is.numeric(y) || is.logical(y))){
   stop("Time series must be numeric")
   }

  Ny <- length(y)

  if(w > Ny){
    stop("The window needs to be smaller than the series length")
  }

  if(w == Ny){

    n <- 4:round(w/k,0)

    r3 <- c()
    r4 <- c()
    r5 <- c()
    ApEn <- c()
    FastApEn <- c()
    SampEn <- c()

    for(i in 1:length(n)){
      ## Divide in to overlapping boxes of size n
      ny <- SlidingWindows(y, n[i])

      for(j in 1:nrow(ny)){

        r3[j] <-     TSEntropies::ApEn(ny[j,], dim = dim, lag = lag, r = r*stats::sd(ny[j,]))
        r4[j] <- TSEntropies::FastApEn(ny[j,], dim = dim, lag = lag, r = r*stats::sd(ny[j,]))
        r5[j] <-   TSEntropies::SampEn(ny[j,], dim = dim, lag = lag, r = r*stats::sd(ny[j,]))

      }

          ApEn[i] <- mean(r3)
      FastApEn[i] <- mean(r4)
        SampEn[i] <- mean(r5)

      ny <- NULL
    }

    return(list(window = w,
                timescale = n,
                ApEn=ApEn,
                FastApEn=FastApEn,
                SampEn = SampEn))
  }


  if(w < Ny){
     sw <- SlidingWindows(y,w)

     n <- 4:round(w/k,0)

     r3 <- c()
     r4 <- c()
     r5 <- c()
  ApEn <- matrix(data = NA, nrow = nrow(sw), ncol = length(n), byrow = TRUE)
 FastApEn <- matrix(data = NA, nrow = nrow(sw), ncol = length(n), byrow = TRUE)
  SampEn <- matrix(data = NA, nrow = nrow(sw), ncol = length(n), byrow = TRUE)

for(h in 1:nrow(sw)){
  for(i in 1:length(n)){
    ## Divide in to overlapping boxes of size n
    ny <- SlidingWindows(sw[h,], n[i])

    for(j in 1:nrow(ny)){

      r3[j] <-       TSEntropies::ApEn(ny[j,], dim = dim, lag = lag, r = r*sd(ny[j,]))
      r4[j] <-   TSEntropies::FastApEn(ny[j,], dim = dim, lag = lag, r = r*sd(ny[j,]))
      r5[j] <-     TSEntropies::SampEn(ny[j,], dim = dim, lag = lag, r = r*sd(ny[j,]))

    }

        ApEn[h, i] <- mean(r3)
    FastApEn[h, i] <- mean(r4)
      SampEn[h, i] <- mean(r5)

    ny <- NULL

   }
  }
  return(list(window = w,
              timescale = n,
              ApEn=ApEn,
              FastApEn=FastApEn,
              SampEn=SampEn))
  }
}

