\name{ddjnonparam_ews}
\alias{ddjnonparam_ews}
\title{Description: Drift Diffusion Jump Nonparametrics Early Warning Signals}
\usage{
  ddjnonparam_ews(timeseries, bandwidth = 0.6, na = 500,
    logtransform = TRUE, interpolate = FALSE)
}
\arguments{
  \item{timeseries}{a numeric vector of the observed
  univariate timeseries values or a numeric matrix where
  the first column represents the time index and the second
  the observed timeseries values. Use vectors/matrices with
  headings.}

  \item{bandwidth}{is the bandwidht of the kernel regressor
  (must be numeric). Default is 0.6.}

  \item{na}{is the number of points for computing the
  kernel (must be numeric). Default is 500.}

  \item{logtransform}{logical. If TRUE data are
  logtransformed prior to analysis as log(X+1). Default is
  FALSE.}

  \item{interpolate}{logical. If TRUE linear interpolation
  is applied to produce a timeseries of equal length as the
  original. Default is FALSE (assumes there are no gaps in
  the timeseries).}
}
\value{
  \code{ddjnonparam_ews} returns an object with elements:

  \item{avec}{is the mesh for which values of the
  nonparametric statistics are estimated.}

  \item{S2.vec}{is the conditional variance of the
  timeseries \code{x} over \code{avec}.}

  \item{TotVar.dx.vec}{is the total variance of \code{dx}
  over \code{avec}.}

  \item{Diff2.vec}{is the diffusion estimated as
  \code{total variance - jumping intensity} vs
  \code{avec}.}

  \item{LamdaZ.vec}{is the jump intensity over
  \code{avec}.}

  \item{Tvec1}{is the timeindex.}

  \item{S2.t}{is the conditional variance of the timeseries
  \code{x} data over \code{Tvec1}.}

  \item{TotVar.t}{is the total variance of \code{dx} over
  \code{Tvec1}.}

  \item{Diff2.t}{is the diffusion over \code{Tvec1}.}

  \item{Lamda.t}{is the jump intensity over \code{Tvec1}.}

  In addition, \code{ddjnonparam_ews} returns a first plot
  with the original timeseries and the residuals after
  first-differencing. A second plot shows the nonparametric
  conditional variance, total variance, diffusion and jump
  intensity over the data, and a third plot the same
  nonparametric statistics over time.
}
\description{
  \code{ddjnonparam_ews} is used to compute
  nonparametrically conditional variance, drift, diffusion
  and jump intensity in a timeseries. It also interpolates
  to obtain the evolution of the nonparametric statistics
  in time.
}
\details{
  The approach is based on estimating terms of a
  drift-diffusion-jump model as a surrogate for the unknown
  true data generating process: [1] \eqn{dx = f(x,\theta)dt
  + g(x,\theta)dW + dJ} Here x is the state variable, f()
  and g() are nonlinear functions, dW is a Wiener process
  and dJ is a jump process. Jumps are large, one-step,
  positive or negative shocks that are uncorrelated in
  time.
}
\examples{
data(foldbif)
output<-ddjnonparam_ews(foldbif,bandwidth=0.6,na=500,
logtransform=TRUE,interpolate=FALSE)
}
\author{
  S. R. Carpenter, modified by V. Dakos
}
\references{
  Carpenter, S. R. and W. A. Brock (2011). "Early warnings
  of unknown nonlinear shifts: a nonparametric approach."
  \emph{Ecology} 92(12): 2196-2201

  Dakos, V., et al (2012)."Methods for Detecting Early
  Warnings of Critical Transitions in Time Series
  Illustrated Using Simulated Ecological Data." \emph{PLoS
  ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010
}
\seealso{
  \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}};
  \code{\link{bdstest_ews}};
  \code{\link{sensitivity_ews}};\code{\link{surrogates_ews}};
  \code{\link{ch_ews}}; \code{\link{movpotential_ews}};
  \code{\link{livpotential_ews}}
}
\keyword{early-warning}

