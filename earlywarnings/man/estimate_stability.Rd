\name{estimate_stability}
\alias{estimate_stability}
\title{estimate_stability}
\usage{
estimate_stability(df, reference.point = NULL)
}
\arguments{
  \item{df}{Input data frame (samples x variables). This
  should contain for each sample the following
  self-explanatory fields: subjectID, time, data.}

  \item{reference.point}{Optional. Calculate stability of
  the data w.r.t. this point.  By default the intermediate
  range is used (min + (max - min)/2)}
}
\value{
Estimated intermediate stability for the given data.
}
\description{
Description: Quantify intermediate stability with respect
to a given reference point.
}
\details{
The stability is quantified by calculating correlation
between (i) deviation from the reference point at the first
time point (ii) deviation between the first and second time
points.  The first and last time point for each subject is
used.  Samples with missing data, and subjects with less
than two time point are excluded.
}
\examples{
#df <- data.frame(list(
  #	  subjectID = rep(paste("subject", 1:50, sep = "-"), each = 2),
  #	  time = rep(1:2, 50),
  #	  data = rnorm(100)))
  # s <- estimate_stability(df, reference.point = NULL)
}
\author{
Leo Lahti \email{leo.lahti@iki.fi}
}
\keyword{internal}

