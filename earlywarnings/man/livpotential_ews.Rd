\name{livpotential_ews}
\alias{livpotential_ews}
\title{Description: Potential Analysis for univariate data}
\usage{
  livpotential_ews(x, grid.points = NULL, std = 1, bw = -1,
    weights = c(), grid.size = 200,
    detection.threshold = 0.002)
}
\arguments{
  \item{x}{Univariate data (vector) for which the
  potentials shall be estimated}

  \item{grid.points}{Data points at which the potential is
  estimated. If not provided, equally sized grid with
  grid.size is used.}

  \item{std}{Standard deviation of the noise (defaults to
  1; this will set scaled potentials)}

  \item{bw}{bandwidth for kernel estimation}

  \item{weights}{optional weights in ksdensity (used by
  movpotentials).}

  \item{grid.size}{Grid size for potential estimation.}

  \item{detection.threshold}{Threshold for local minima to
  be discarded.}
}
\value{
  \code{livpotential} returns a list with the following
  elements:

  \item{xi}{the grid of points on which the potential is
  estimated}

  \item{pot}{the actual value of the potential}

  \item{min.inds}{indices of the grid points at which the
  density has minimum values; (-potentials; neglecting
  local optima)}

  \item{max.inds}{indices the grid points at which the
  density has maximum values; (-potentials; neglecting
  local optima)}

  \item{bw}{bandwidth of kernel used}

  \item{min.points}{grid point values at which the density
  has minimum values; (-potentials; neglecting local
  optima)}

  \item{max.points}{grid point values at which the density
  has maximum values; (-potentials; neglecting local
  optima)}
}
\description{
  \code{livpotential_ews} performs one-dimensional
  potential estimation derived from a uni-variate
  timeseries
}
\examples{
data(foldbif)
res <- livpotential_ews(foldbif)
plot(res$grid.points, res$pot)
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

