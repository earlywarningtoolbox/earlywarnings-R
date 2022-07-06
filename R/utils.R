#' @import Kendall
#' @import KernSmooth
#' @import lmtest
#' @import nortest
#' @import som
#' @import spam
#' @import stats
NULL

#' Get group assigment indices for univariate data points, given cluster break points
#'
#' @param x Univariate data vector
#' @param breakpoints Cluster breakpoints
#' @return A vector of cluster indices
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords early-warning
UnivariateGrouping <- function(x, breakpoints) {
    g <- rep.int(NA, length(x))
    mps <- c(breakpoints, Inf)
    for (i in 1:length(mps)) {
        g[x <= mps[[i]] & is.na(g)] <- i
    }
    g
} 
