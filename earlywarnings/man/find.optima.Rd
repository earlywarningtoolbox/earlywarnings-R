\name{find.optima}
\alias{find.optima}
\title{Description: find.optima}
\usage{
  find.optima(f, detection.multiplier = 0, bw, x,
    detection.limit = 0)
}
\arguments{
  \item{f}{density}

  \item{detection.multiplier}{multiple of kernel height to
  be used as the detection threshold}

  \item{bw}{bandwidth}

  \item{x}{original data}

  \item{detection.limit}{ignore maxima that are below
  detection.limit given by multiple of kernel height

  Returns:}
}
\value{
  A list with the following elements: min minima max maxima
}
\description{
  Detect optima, excluding very local optima below
  detection.threshold
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

