\name{movpotential_ews}
\alias{movpotential_ews}
\title{Description: Moving Average Potential}
\usage{
  movpotential_ews(X, param = NULL, bw = -1,
    detection.threshold = 0.002, std = 1, grid.size = 50,
    plot.cutoff = 0.5)
}
\arguments{
  \item{X}{a vector of the X observations of the state
  variable of interest}

  \item{param}{parameter values corresponding to the
  observations in X}

  \item{bw}{Bandwidth for smoothing kernels. Automatically
  determined by default.}

  \item{detection.threshold}{Threshold for local minima to
  be discarded}

  \item{std}{std}

  \item{grid.size}{number of evaluation points; number of
  steps between min and max potential; also used as kernel
  window size}

  \item{plot.cutoff}{cuttoff for potential minima and
  maxima in visualization

  Returns:}
}
\value{
  A list with the following elements:

  \item{pars}{values of the covariate parameter as matrix}

  \item{xis}{values of the x as matrix}

  \item{pots}{smoothed potentials}

  \item{mins}{minima in the densities (-potentials;
  neglecting local optima)}

  \item{maxs}{maxima in densities (-potentials; neglecting
  local optima)}

  \item{plot}{an object that displays the potential
  estimated in 2D}
}
\description{
  \code{movpotential_ews} reconstructs a potential derived
  from data along a gradient of a given parameter the
  \code{movpotential_ews} calculates the potential for
  values that correspond to a particular parameter. see ref
  below
}
\examples{
X = c(rnorm(1000, mean = 0), rnorm(1000, mean = -2), rnorm(1000, mean = 2))
 param = seq(0,5,length=3000)
 res <- movpotential_ews(X, param)
}
\author{
  Based on Matlab code from Egbert van Nes modified by Leo
  Lahti. Implemented in early warnings package by V. Dakos.
}
\references{
  Hirota, M., Holmgren, M., van Nes, E.H. & Scheffer, M.
  (2011). Global resilience of tropical forest and savanna
  to critical transitions. \emph{Science}, 334, 232-235.
}
\seealso{
  \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}};
  \code{\link{bdstest_ews}};
  \code{\link{sensitivity_ews}};\code{\link{surrogates_ews}};
  \code{\link{ch_ews}}; \code{livpotential_ews}
}
\keyword{early-warning}

