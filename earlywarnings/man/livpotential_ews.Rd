\name{livpotential_ews}
\alias{livpotential_ews}
\title{Description: Potential Analysis}
\usage{
  livpotential_ews(x, std = 1, bw = -1, xi = NULL,
    weights = c(), grid.size = 200)
}
\arguments{
  \item{x}{data vector}

  \item{std}{the standard deviation of the noise (defaults
  to 1, so then you use scaled potentials}

  \item{bw}{bandwidth for kernel estimation}

  \item{xi}{x values at which the potential is estimated}

  \item{weights}{optional weights in ksdensity (used by
  movpotentials).}

  \item{grid.size}{grid size}
}
\value{
  \code{livpotential} returns a list with the following
  elements:

  \item{xi}{the grid of points on which the potential is
  estimated}

  \item{pot}{the actual value of the potential}

  \item{minima}{the grid points at which the potential has
  minimum values}

  \item{maxima}{the grid points at which the potential has
  maximum values}

  \item{bw}{bandwidth of kernel used}
}
\description{
  \code{livpotential_ews} performs one-dimensional
  potential estimation derived from a uni-variate
  timeseries
}
\details{
  see ref below
}
\examples{
data(foldbif)
res <- livpotential_ews(foldbif)
plot(res$xi, res$pot)
}
\author{
  Based on Matlab code from Egbert van Nes modified by Leo
  Lahti. Implemented in early warnings package by V. Dakos.
}
\references{
  Livina, VN, F Kwasniok, and TM Lenton, 2010. Potential
  analysis reveals changing number of climate states during
  the last 60 kyr . \emph{Climate of the Past}, 6, 77-82.

  Dakos, V., et al (2012)."Methods for Detecting Early
  Warnings of Critical Transitions in Time Series
  Illustrated Using Simulated Ecological Data." \emph{PLoS
  ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010
}
\seealso{
  \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}};
  \code{\link{bdstest_ews}};
  \code{\link{sensitivity_ews}};\code{\link{surrogates_ews}};
  \code{\link{ch_ews}};\code{\link{movpotential_ews}}
}
\keyword{early-warning}

