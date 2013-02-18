\name{ch_ews}
\alias{ch_ews}
\title{Description: Conditional Heteroskedasticity}
\usage{
  ch_ews(timeseries, winsize = 10, alpha = 0.1,
    optim = TRUE, lags = 4, logtransform = FALSE,
    interpolate = FALSE)
}
\arguments{
  \item{timeseries}{a numeric vector of the observed
  timeseries values or a numeric matrix where the first
  column represents the time index and the second the
  observed timeseries values. Use vectors/matrices with
  headings.}

  \item{winsize}{is length of the rolling window expressed
  as percentage of the timeseries length (must be numeric
  between 0 and 100). Default is 10\%.}

  \item{alpha}{is the significance threshold (must be
  numeric). Default is 0.1.}

  \item{optim}{logical. If TRUE an autoregressive model is
  fit to the data within the rolling window using AIC
  optimization. Otherwise an autoregressive model of
  specific order \code{lags} is selected.}

  \item{lags}{is a parameter that determines the specific
  order of an autoregressive model to fit the data. Default
  is 4.}

  \item{logtransform}{logical. If TRUE data are
  logtransformed prior to analysis as log(X+1). Default is
  FALSE.}

  \item{interpolate}{logical. If TRUE linear interpolation
  is applied to produce a timeseries of equal length as the
  original. Default is FALSE (assumes there are no gaps in
  the timeseries).}
}
\value{
  \code{ch_ews} returns a matrix that contains:

  \item{time}{the time index.}

  \item{r.squared}{the R2 values of the regressed
  residuals.}

  \item{critical.value}{the chi-square critical value based
  on the desired \code{alpha} level for 1 degree of freedom
  divided by the number of residuals used in the
  regression.}

  \item{test.result}{logical. It indicates whether
  conditional heteroskedasticity was significant.}

  \item{ar.fit.order}{the order of the specified
  autoregressive model- only informative if \code{optim}
  FALSE was selected.}

  In addition, \code{ch_ews} plots the original timeseries
  and the R2 where the level of significance is also
  indicated.
}
\description{
  \code{ch_ews} is used to estimate changes in conditional
  heteroskedasticity within rolling windows along a
  timeseries
}
\details{
  see ref below
}

\author{
  T. Cline, modified by V. Dakos
}
\references{
  Seekell, D. A., et al (2011). "Conditional
  heteroscedasticity as a leading indicator of ecological
  regime shifts." \emph{American Naturalist} 178(4):
  442-451

  Dakos, V., et al (2012)."Methods for Detecting Early
  Warnings of Critical Transitions in Time Series
  Illustrated Using Simulated Ecological Data." \emph{PLoS
  ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010
}
\seealso{
  \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}};
  \code{\link{bdstest_ews}}; \code{\link{sensitivity_ews}};
  \code{\link{surrogates_ews}}; \code{\link{ch_ews}};
  \code{movpotential_ews}; \code{livpotential_ews}
}
\keyword{early-warning}

