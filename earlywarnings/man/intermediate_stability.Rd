\name{intermediate_stability}
\alias{intermediate_stability}
\title{Description: Quantify intermediate stability}
\usage{
intermediate_stability(dat, meta, reference.point = NULL)
}
\arguments{
  \item{dat}{Input data matrix (variables x samples)}

  \item{meta}{Metadata (samples x factors). This should
  contain for each sample the following self-explanatory
  fields: subjectID, time, data.}

  \item{reference.point}{Calculate stability of the data
  w.r.t. this point.  By default the intermediate range is
  used (min + (max - min)/2)}
}
\value{
A vector listing the intermediate stability for each
variable
}
\description{
Description: Quantify intermediate stability
}
\examples{
# Create simulated example data
  dat <- matrix(rnorm(1000), nrow = 10)
  rownames(dat) <- paste("Variable", 1:nrow(dat), sep = "")
  colnames(dat) <- paste("Sample", 1:ncol(dat), sep = "")
  meta <- data.frame(list(
    	  sampleID = colnames(dat),
  	  subjectID = rep(paste("subject", 1:50, sep = "-"), each = 2),
	  time = rep(1:2, 50)))
  # Intentionally make point very unstable around 0
  dat[1, meta$time == 2] <- 1/abs(dat[1, meta$time == 1])
  s <- intermediate_stability(dat, meta, reference.point = 0)
}
\author{
Leo Lahti \email{leo.lahti@iki.fi}
}
\keyword{early-warning}

