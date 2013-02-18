\name{bdstest_ews}
\alias{bdstest_ews}
\title{Description: BDS test Early Warning Signals}
\usage{
  bdstest_ews(timeseries, ARMAoptim = TRUE,
    ARMAorder = c(1, 0), GARCHorder = c(0, 1), embdim = 3,
    epsilon = c(0.5, 0.75, 1), boots = 1000,
    logtransform = FALSE, interpolate = FALSE)
}
\arguments{
  \item{timeseries}{a numeric vector of the observed
  univariate timeseries values or a numeric matrix where
  the first column represents the time index and the second
  the observed timeseries values. Use vectors/matrices with
  headings.}

  \item{ARMAoptim}{is the order of the \code{ARMA(p,q)}
  model to be fitted on the original timeseries. If TRUE
  the best ARMA model based on AIC is applied. If FALSE the
  \code{ARMAorder} is used.}

  \item{ARMAorder}{is the order of the \code{AR(p)} and
  \code{MA(q)} process to be fitted on the original
  timeseries. Default is \code{p=1} \code{q=0}.}

  \item{GARCHorder}{fits a GARCH model on the original
  timeseries where \code{GARCHorder[1]} is the GARCH part
  and \code{GARCHorder[2]} is the ARCH part.}

  \item{embdim}{is the embedding dimension (2, 3,...
  \code{embdim}) up to which the BDS test will be estimated
  (must be numeric). Default value is 3.}

  \item{epsilon}{is a numeric vector that is used to scale
  the standard deviation of the timeseries. The BDS test is
  computed for each element of epsilon. Default is 0.5,
  0.75 and 1.}

  \item{boots}{is the number of bootstraps performed to
  estimate significance p values for the BDS test. Default
  is 1000.}

  \item{logtransform}{logical. If TRUE data are
  logtransformed prior to analysis as log(X+1). Default is
  FALSE.}

  \item{interpolate}{logical. If TRUE linear interpolation
  is applied to produce a timeseries of equal length as the
  original. Default is FALSE (assumes there are no gaps in
  the timeseries).}
}
\value{
  \code{bdstest_ews} returns output on the R console that
  summarizes the BDS test statistic for all embedding
  dimensions and \code{epsilon} values used, and for
  first-differenced data, ARMA(p.q) residuals, and
  GARCH(0,1) residuals). Also the significance p values are
  returned estimated both by comparing to a standard normal
  distribution and by bootstrapping.

  In addition, \code{bdstest_ews} returns a plot with the
  original timeseries, the residuals after
  first-differencing, and fitting the ARMA(p,q) and
  GARCH(0,1) models. Also the autocorrelation
  \code{\link{acf}} and partial autocorrelation
  \code{\link{pacf}} functions are estimated serving as
  guides for the choice of lags of the linear models fitted
  to the data.
}
\description{
  \code{bdstest_ews} is used to estimate the BDS statistic
  to detect nonlinearity in the residuals of a timeseries
  after first-difference detrending, fitting an ARMA(p,q)
  model, and fitting a GARCH(0,1) model. The function is
  making use of \code{bds.test}.
}
\details{
  See also \code{bds.test{tseries}} for more details. The
  function requires the installation of packages
  \code{tseries} and \code{quadprog} that are not available
  under Linux and need to be manually installed under
  Windows.

  Example to run after installing the mentioned packages:

  data(foldbif)
  bdstest_ews(foldbif,ARMAoptim=FALSE,ARMAorder=c(1,0),embdim=3,epsilon=0.5,
  boots=200,logtransform=FALSE,interpolate=FALSE)
}
\author{
  S. R. Carpenter, modified by V. Dakos
}
\references{
  J. B. Cromwell, W. C. Labys and M. Terraza (1994):
  Univariate Tests for Time Series Models, Sage, Thousand
  Oaks, CA, pages 32-36.

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

