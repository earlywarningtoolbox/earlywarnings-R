\name{movpotential_ews}
\alias{movpotential_ews}
\title{Description: Moving Average Potential}
\usage{
  movpotential_ews(X, param, sdwindow = NULL, bw = -1,
    minparam = NULL, maxparam = NULL, npoints = 50,
    thres = 0.002, std = 1, grid.size = 200, cutoff = 0.5)
}
\arguments{
  \item{X}{a vector of the X observations of the state
  variable of interest}

  \item{param}{parameter values that correspond to the X
  observations}

  \item{sdwindow}{window for smoothing kernels (over the
  \code{param} axis)}

  \item{bw}{bandwidth used for smoothing kernels}

  \item{minparam}{minimum value of parameter on which to
  estimate potential}

  \item{maxparam}{maximum value of parameter on which to
  estimate potential}

  \item{npoints}{number of potentials}

  \item{thres}{threshold for local minima to be discarded}

  \item{std}{std}

  \item{grid.size}{number of evaluation points}

  \item{cutoff}{the cuttof value to estimate minima and
  maxima in the potential

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
 res <- movpotential_ews(X, param, npoints = 100, thres = 0.003)
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

