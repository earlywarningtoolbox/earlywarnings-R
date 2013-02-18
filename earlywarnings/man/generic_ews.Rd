\name{generic_ews}
\alias{generic_ews}
\title{Description: Generic Early Warning Signals}
\usage{
  generic_ews(timeseries, winsize = 50,
    detrending = c("no", "gaussian", "linear", "first-diff"),
    bandwidth = NULL, logtransform = FALSE,
    interpolate = FALSE, AR_n = FALSE,
    powerspectrum = FALSE)
}
\arguments{
  \item{timeseries}{a numeric vector of the observed
  univariate timeseries values or a numeric matrix where
  the first column represents the time index and the second
  the observed timeseries values. Use vectors/matrices with
  headings. If the powerspectrum is to be plotted as well,
  the timeseries lenght should be even number.}

  \item{winsize}{is the size of the rolling window
  expressed as percentage of the timeseries length (must be
  numeric between 0 and 100). Default is 50\%.}

  \item{bandwidth}{is the bandwidth used for the Gaussian
  kernel when gaussian filtering is applied. It is
  expressed as percentage of the timeseries length (must be
  numeric between 0 and 100). Alternatively it can be given
  by the bandwidth selector \code{\link{bw.nrd0}}
  (Default).}

  \item{detrending}{the timeseries can be
  detrended/filtered prior to analysis. There are four
  options: \code{gaussian} filtering, \code{linear}
  detrending and \code{first-differencing}. Default is
  \code{no} detrending.}

  \item{logtransform}{logical. If TRUE data are
  logtransformed prior to analysis as log(X+1). Default is
  FALSE.}

  \item{interpolate}{logical. If TRUE linear interpolation
  is applied to produce a timeseries of equal length as the
  original. Default is FALSE (assumes there are no gaps in
  the timeseries).}

  \item{AR_n}{logical. If TRUE the best fitted AR(n) model
  is fitted to the data. Default is FALSE.}

  \item{powerspectrum}{logical. If TRUE the power spectrum
  within each rolling window is plotted. Default is FALSE.}
}
\value{
  \code{generic_ews} returns a matrix that contains:

  \item{tim}{the time index.}

  \item{ar1}{the \code{autoregressive coefficient ar(1)} of
  a first order AR model fitted on the data within the
  rolling window.}

  \item{sd}{the \code{standard deviation} of the data
  estimated within each rolling window.}

  \item{sk}{the \code{skewness} of the data estimated
  within each rolling window.}

  \item{kurt}{the \code{kurtosis} of the data estimated
  within each rolling window.}

  \item{cv}{the \code{coefficient of variation} of the data
  estimated within each rolling window.}

  \item{returnrate}{the return rate of the data estimated
  as \code{1-ar(1)} cofficient within each rolling window.}

  \item{densratio}{the \code{density ratio} of the power
  spectrum of the data estimated as the ratio of low
  frequencies over high frequencies within each rolling
  window.}

  \item{acf1}{the \code{autocorrelation at first lag} of
  the data estimated within each rolling window.}

  In addition, \code{generic_ews} returns three plots. The
  first plot contains the original data, the
  detrending/filtering applied and the residuals (if
  selected), and all the moment statistics. For each
  statistic trends are estimated by the nonparametric
  Kendall tau correlation.  The second plot, if asked,
  quantifies resilience indicators fitting AR(n) selected
  by the Akaike Information Criterion. The third plot, if
  asked, is the power spectrum estimated by
  \code{\link{spec.ar}} for all frequencies within each
  rolling window.
}
\description{
  \code{generic_ews} is used to estimate statistical
  moments within rolling windows along a timeserie
}
\details{
  see ref below
}
\examples{
data(foldbif)
 out=generic_ews(foldbif,winsize=50,detrending="gaussian",
 bandwidth=5,logtransform=FALSE,interpolate=FALSE)
}
\author{
  Vasilis Dakos \email{vasilis.dakos@gmail.com}
}
\references{
  Ives, A. R. (1995). "Measuring resilience in stochastic
  systems." \emph{Ecological Monographs} 65: 217-233

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

