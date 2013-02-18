\name{sensitivity_ews}
\alias{sensitivity_ews}
\title{Description: Sensitivity Early Warning Signals}
\usage{
  sensitivity_ews(timeseries,
    indicator = c("ar1", "sd", "acf1", "sk", "kurt", "cv", "returnrate", "densratio"),
    winsizerange = c(25, 75), incrwinsize = 25,
    detrending = c("no", "gaussian", "linear", "first-diff"),
    bandwidthrange = c(5, 100), incrbandwidth = 20,
    logtransform = FALSE, interpolate = FALSE)
}
\arguments{
  \item{timeseries}{a numeric vector of the observed
  univariate timeseries values or a numeric matrix where
  the first column represents the time index and the second
  the observed timeseries values. Use vectors/matrices with
  headings.}

  \item{indicator}{is the statistic (leading indicator)
  selected for which the sensitivity analysis is perfomed.
  Currently, the indicators supported are: \code{ar1}
  autoregressive coefficient of a first order AR model,
  \code{sd} standard deviation, \code{acf1} autocorrelation
  at first lag, \code{sk} skewness, \code{kurt} kurtosis,
  \code{cv} coeffcient of variation, \code{returnrate}, and
  \code{densratio} density ratio of the power spectrum at
  low frequencies over high frequencies.}

  \item{winsizerange}{is the range of the rolling window
  sizes expressed as percentage of the timeseries length
  (must be numeric between 0 and 100). Default is 25\% -
  75\%.}

  \item{incrwinsize}{increments the rolling window size
  (must be numeric between 0 and 100). Default is 25.}

  \item{detrending}{the timeseries can be
  detrended/filtered. There are three options:
  \code{gaussian} filtering, \code{linear} detrending and
  \code{first-differencing}. Default is \code{no}
  detrending.}

  \item{bandwidthrange}{is the range of the bandwidth used
  for the Gaussian kernel when gaussian filtering is
  selected. It is expressed as percentage of the timeseries
  length (must be numeric between 0 and 100). Default is
  5\% - 100\%.}

  \item{incrbandwidth}{is the size to increment the
  bandwidth used for the Gaussian kernel when gaussian
  filtering is applied. It is expressed as percentage of
  the timeseries length (must be numeric between 0 and
  100). Default is 20.}

  \item{logtransform}{logical. If TRUE data are
  logtransformed prior to analysis as log(X+1). Default is
  FALSE.}

  \item{interpolate}{logical. If TRUE linear interpolation
  is applied to produce a timeseries of equal length as the
  original. Default is FALSE (assumes there are no gaps in
  the timeseries).}
}
\value{
  \code{sensitivity_ews} returns a matrix that contains the
  Kendall tau rank correlation estimates for the rolling
  window sizes (rows) and bandwidths (columns), if
  \code{gaussian filtering} is selected.

  In addition, \code{sensitivity_ews} returns a plot with
  the Kendall tau estimates and their p-values for the
  range of rolling window sizes used, together with a
  histogram of the distributions of the statistic and its
  significance. When \code{gaussian filtering} is chosen, a
  contour plot is produced for the Kendall tau estimates
  and their p-values for the range of both rolling window
  sizes and bandwidth used. A reverse triangle indicates
  the combination of the two parameters for which the
  Kendall tau was the highest
}
\description{
  \code{sensitivity_ews} is used to estimate trends in
  statistical moments for different sizes of rolling
  windows along a timeseries. The trends are estimated by
  the nonparametric Kendall tau correlation coefficient.
}
\details{
  see ref below
}

\author{
  Vasilis Dakos \email{vasilis.dakos@gmail.com}
}
\references{
  Dakos, V., et al (2008). "Slowing down as an early
  warning signal for abrupt climate change."
  \emph{Proceedings of the National Academy of Sciences}
  105(38): 14308-14312

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
  \code{\link{livpotential_ews}}
}
\keyword{early-warning}

