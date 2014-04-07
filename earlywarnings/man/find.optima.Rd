\name{find.optima}
\alias{find.optima}
\title{Description: find.optima}
\usage{
find.optima(f, detection.threshold = 0, bw, x, detection.limit = 0)
}
\arguments{
  \item{f}{density}

  \item{detection.threshold}{detection threshold for peaks}

  \item{bw}{bandwidth}

  \item{x}{original data}

  \item{detection.limit}{ignore maxima that are below this
  value

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
See citation('TBA')
}
\keyword{utilities}

