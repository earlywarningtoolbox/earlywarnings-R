\name{find.optima}
\alias{find.optima}
\title{Description: find.optima}
\usage{
  find.optima(fpot, detection.limit = 0, bw, x)
}
\arguments{
  \item{fpot}{potential}

  \item{detection.limit}{detection threshold will be
  determined by multiplying this scalar with kernel
  height.}

  \item{bw}{bandwidth}

  \item{x}{original data

  Returns:}
}
\value{
  A list with the following elements: min potential minima
  max potential maxima
}
\description{
  Detect optima from the potential, excluding very local
  optima below detection.threshold
}
\details{
  Arguments:
}
\examples{
#
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("TBA")
}
\keyword{utilities}

