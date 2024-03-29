% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/surrogates_ews.R
\name{surrogates_ews}
\alias{surrogates_ews}
\title{Surrogates Early Warning Signals}
\usage{
surrogates_ews(
  timeseries,
  indicator = c("ar1", "sd", "acf1", "sk", "kurt", "cv", "returnrate", "densratio"),
  winsize = 50,
  detrending = c("no", "gaussian", "loess", "linear", "first-diff"),
  bandwidth = NULL,
  span = NULL,
  degree = NULL,
  boots = 100,
  logtransform = FALSE,
  interpolate = FALSE
)
}
\arguments{
\item{timeseries}{a numeric vector of the observed univariate timeseries values or a numeric matrix where the first column represents the time index and the second the observed timeseries values. Use vectors/matrices with headings.}

\item{indicator}{is the statistic (leading indicator) selected for which the surrogate timeseries are produced. Currently, the indicators supported are: \code{ar1} autoregressive coefficient of a first order AR model, \code{sd} standard deviation, \code{acf1} autocorrelation at first lag, \code{sk} skewness, \code{kurt} kurtosis, \code{cv} coeffcient of variation, \code{returnrate}, and \code{densratio} density ratio of the power spectrum at low frequencies over high frequencies.}

\item{winsize}{is the size of the rolling window expressed as percentage of the timeseries length (must be numeric between 0 and 100). Default valuise 50\%.}

\item{detrending}{the timeseries can be detrended/filtered prior to analysis. There are three options: \code{gaussian} filtering, \code{loess} fitting, \code{linear} detrending and \code{first-diff}erencing. Default is \code{no} detrending.}

\item{bandwidth}{is the bandwidth used for the Gaussian kernel when gaussian filtering is selected. It is expressed as percentage of the timeseries length (must be numeric between 0 and 100). Alternatively it can be given by the bandwidth selector \code{\link{bw.nrd0}} (Default).}

\item{span}{parameter that controls the degree of smoothing (numeric between 0 and 100, Default 25). see more on loess{stats}}

\item{degree}{the degree of polynomial to be used for when loess fitting is applied, normally 1 or 2 (Default). see more on loess{stats}}

\item{boots}{the number of surrogate data. Default is 100.}

\item{logtransform}{logical. If TRUE data are logtransformed prior to analysis as log(X+1). Default is FALSE.}

\item{interpolate}{logical. If TRUE linear interpolation is applied to produce a timeseries of equal length as the original. Default is FALSE (assumes there are no gaps in the timeseries).}
}
\value{
\code{surrogates_ews} returns a matrix that contains:
  Kendall tau estimate original the trends of the original timeseries;
  Kendall tau p-value original the p-values of the trends of the original timeseries;
  Kendall tau estimate surrogates the trends of the surrogate timeseries;
  Kendall tau p-value surrogates the associated p-values of the trends of the surrogate timeseries;
  significance p the p-value for the original Kendall tau rank correlation estimate compared to the surrogates;
}
\description{
\code{\link{surrogates_ews}} is used to estimate distributions of trends in statistical moments from different surrogate timeseries generated after fitting an ARMA(p,q) model on the data. The trends are estimated by the nonparametric Kendall tau correlation coefficient and can be compared to the trends estimated in the original timeseries to produce probabilities of false positives.
}
\details{
In addition, \code{surrogates_ews} returns a plot with the distribution of the surrogate Kendall tau estimates and the Kendall tau estimate of the original series. Vertical lines indicate the 5\% and 95\% significance levels.
}
\examples{
data(foldbif) 
output <- surrogates_ews(foldbif,indicator='sd',winsize=50,detrending='gaussian', bandwidth=10,
             boots=200, logtransform=FALSE,interpolate=FALSE)
}
\references{
Dakos, V., et al (2008). 'Slowing down as an early warning signal for abrupt climate change.' \emph{Proceedings of the National Academy of Sciences} 105(38): 14308-14312 

Dakos, V., et al (2012).'Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data.' \emph{PLoS ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010
}
\author{
Vasilis Dakos \email{vasilis.dakos@gmail.com}
}
\keyword{early-warning}
