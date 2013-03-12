\name{qda_ews}
\alias{qda_ews}
\title{Description: Quick Detection Analysis for Generic Early Warning Signals}
\usage{
  qda_ews(timeseries, param = NULL, winsize = 50,
    detrending = c("no", "gaussian", "linear", "first-diff"),
    bandwidth = NULL, boots = 100, s_level = 0.05,
    cutoff = 0.05, detection.threshold = 0.002,
    grid.size = 50, logtransform = FALSE,
    interpolate = FALSE)
}
\arguments{
  \item{timeseries}{a numeric vector of the observed
  univariate timeseries values or a numeric matrix where
  the first column represents the time index and the second
  the observed timeseries values. Use vectors/matrices with
  headings.}

  \item{param}{values corresponding to observations in
  timeseries}

  \item{winsize}{is the size of the rolling window
  expressed as percentage of the timeseries length (must be
  numeric between 0 and 100). Default is 50\%.}

  \item{detrending}{the timeseries can be
  detrended/filtered prior to analysis. There are four
  options: \code{gaussian} filtering, \code{linear}
  detrending and \code{first-differencing}. Default is
  \code{no} detrending.}

  \item{bandwidth}{is the bandwidth used for the Gaussian
  kernel when gaussian filtering is applied. It is
  expressed as percentage of the timeseries length (must be
  numeric between 0 and 100). Alternatively it can be given
  by the bandwidth selector \code{\link{bw.nrd0}}
  (Default).}

  \item{boots}{the number of surrogate data to generate
  from fitting an ARMA(p,q) model. Default is 100.}

  \item{s_level}{significance level. Default is 0.05.}

  \item{cutoff}{the cutoff value to visualize the potential
  landscape}

  \item{detection.threshold}{detection threshold for
  potential minima}

  \item{grid.size}{grid size (for potential analysis)}

  \item{logtransform}{logical. If TRUE data are
  logtransformed prior to analysis as log(X+1). Default is
  FALSE.}

  \item{interpolate}{logical. If TRUE linear interpolation
  is applied to produce a timeseries of equal length as the
  original. Default is FALSE (assumes there are no gaps in
  the timeseries).

  \code{qda_ews} returns three plots. The first plot
  contains the original data, the detrending/filtering
  applied and the residuals (if selected), autocorrelation
  and variance. For each statistic trends are estimated by
  the nonparametric Kendall tau correlation.  The second
  plot, returns a histogram of the distributions of the
  Kendall trend statistic for autocorrelation and variance
  estimated on the surrogated data. Vertical lines
  represent the level of significance, whereas the black
  dots the actual trend found in the time series. The third
  plot is the reconstructed potential landscape in 2D.}
}
\description{
  \code{qda_ews} is used to estimate autocorrelation,
  variance within rolling windows along a timeseries, test
  the significance of their trends, and reconstruct the
  potential landscape of the timeseries
}
\details{
  see ref below
}
\examples{
data(foldbif)
out <- qda_ews(foldbif, param = NULL, winsize = 50, detrending="gaussian", bandwidth=NULL, boots = 50, s_level = 0.05, cutoff=0.05, detection.threshold = 0.002, grid.size = 50, logtransform=FALSE, interpolate=FALSE)
}
\author{
  Vasilis Dakos, Leo Lahti, March 1, 2013
  \email{vasilis.dakos@gmail.com}
}
\references{
  Dakos, V., et al (2012)."Methods for Detecting Early
  Warnings of Critical Transitions in Time Series
  Illustrated Using Simulated Ecological Data." \emph{PLoS
  ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010
}
\seealso{
  \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}};
  \code{\link{bdstest_ews}}; \code{\link{sensitivity_ews}};
  \code{\link{surrogates_ews}}; \code{\link{ch_ews}};
  \code{\link{movpotential_ews}};
  \code{\link{livpotential_ews}};
}
\keyword{early-warning}

