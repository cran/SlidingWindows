#' @title Detrended Fluctuation Analysis with sliding windows.
#'
#' @description This function generates scaling exponents (long-range correlations) of
#' a univariate time series with sliding windows approach.
#'
#' @details This function include following measures:
#'
#'     alpha_dfa, se_alpha_dfa, r2_alpha_dfa.
#'
#' @param y A vector contaning univariate time series.
#' @param w an integer value indicating the size of the window \eqn{w < length(y)}.
#'          If \eqn{w = length(y)}, will be computed the function will not slide.
#' @param k an integer value indicating the boundary of the division \eqn{(N/k)}.
#'          The smallest value of \eqn{k} is \eqn{4}.
#' @param npoints The number of different window sizes that will be used to estimate the Fluctuation function in each zone. See nonlinearTseries package.
#'
#' @return A list contaning "w", "alpha_dfa", "se_alpha_dfa", "r2_alpha_dfa".
#'
#' @examples
#' y <- rnorm(100)
#' dfa_SlidingWindows(y,w=99,k=10,npoints=15)
#'
#' @references
#' GUEDES, E.F.;FERREIRA, P.;DIONISIO, A.; ZEBENDE,G.F. An econophysics approach to study the effect of BREXIT referendum on European Union stock markets. PHYSICA A, v.523, p.1175-1182, 2019. doi = "https://doi.org/10.1016/j.physa.2019.04.132".
#'
#' FERREIRA, P.; DIONISIO, A.;GUEDES, E.F.; ZEBENDE, G.F. A sliding windows approach to analyse the evolution of bank shares in the European Union. PHYSICA A, v.490, p.1355-1367, 2018. doi = "https://doi.org/10.1016/j.physa.2017.08.095".
#'
#' @importFrom nonlinearTseries dfa
#' @importFrom stats coef lm
#'
#' @export
dfa_SlidingWindows <- function(y,w,k,npoints){
  if(!(is.null(y) || is.numeric(y) || is.logical(y))){
   stop("Time series must be numeric")
   }

  Ny <- length(y)

  if(w > Ny){
    stop("The window needs to be smaller than the series length")
  }

  if(w == Ny){
      dfa <- nonlinearTseries::dfa(y,
                               window.size.range=c(4,round(w/k,0)),
                               npoints=npoints,
                               do.plot=FALSE)
             model <- lm(log10(dfa$fluctuation.function)~log10(dfa$window.sizes))
         alpha_dfa <- coef(summary(model))[2, "Estimate"]
      se_alpha_dfa <- coef(summary(model))[2, "Std. Error"]
      r2_alpha_dfa <- summary(model)$r.squared
  return(list(alpha_dfa=alpha_dfa, se_alpha_dfa=se_alpha_dfa, r2_alpha_dfa=r2_alpha_dfa))
  }

  if(w < Ny){
  sw <- SlidingWindows(y,w)
  alpha_dfa <- c()
  se_alpha_dfa <- c()
  r2_alpha_dfa <- c()
  for(i in 1:nrow(sw)){
     dfa <- nonlinearTseries::dfa(sw[i,],
                               window.size.range=c(4,round(w/k,0)),
                               npoints=npoints,
                               do.plot=FALSE)
               model <- lm(log10(dfa$fluctuation.function)~log10(dfa$window.sizes))
        alpha_dfa[i] <- coef(summary(model))[2, "Estimate"]
     se_alpha_dfa[i] <- coef(summary(model))[2, "Std. Error"]
     r2_alpha_dfa[i] <- summary(model)$r.squared
  }
  return(list(window = w, alpha_dfa=alpha_dfa, se_alpha_dfa=se_alpha_dfa, r2_alpha_dfa=r2_alpha_dfa))
 }
}
