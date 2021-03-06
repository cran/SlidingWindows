% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmca_SlidingWindows.R
\name{dmca.SlidingWindows}
\alias{dmca.SlidingWindows}
\title{DMCA coefficient with sliding windows.}
\usage{
dmca.SlidingWindows(x, y, w = 98, k = 10)
}
\arguments{
\item{x}{A vector containing univariate time series.}

\item{y}{A vector containing univariate time series.}

\item{w}{An integer value indicating the window size \eqn{w < length(y)}.
If \eqn{w = length(y)}, will be computed the function will not slide.}

\item{k}{An integer value indicating the boundary of the division \eqn{(N/k)}.
The smallest value of \eqn{k} is \eqn{4}.}
}
\value{
A list containing "w", "timescale", "dmca".
}
\description{
This function generates Detrending moving-average cross-correlation coefficient of two time series with sliding windows approach.
}
\details{
This function include following measures: w, timescale, dmca
}
\examples{
x <- rnorm(100)
y <- rnorm(100)
dmca.SlidingWindows(x,y,w=99,k=10)

}
\references{
KRISTOUFEK, L. Detrending moving-average cross-correlation coefficient: Measuring cross-correlations between non-stationary series. PHYSICA A, v.406, p.169-175, 2014. doi="doi.org/10.1016/j.physa.2014.03.015".
}
