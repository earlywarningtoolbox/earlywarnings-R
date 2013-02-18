# Description: Plot Potential
#
# Visualization of the potential function from the movpotential function
#
#  Arguments:
#    @param res output from movpotential function
#    @param title title text
#    @param xlab.text xlab text
#    @param ylab.text ylab text
#    @param cutoff parameter determining the upper limit of potential for visualizations
# 
# Returns:
#   @return \item{ggplot2}{potential plotted}
#
# @export
#
# @references Dakos, V., et al (2012)."Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data." \emph{PLoS ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010
# @author L. Lahti
# @examples #
#
# @keywords early-warning

PlotPotential <- function (res, title = "", xlab.text, ylab.text, cutoff = 0.5) {

  library(akima)
  library(ggplot2)
  cut.potential <- max(apply(res$pots, 1, min)) + cutoff*abs(max(apply(res$pots, 1, min))) # Ensure all minima are visualized
  pots <- res$pots
  pots[pots > cut.potential] <- cut.potential

  # Static contour
  # Interpolate potential grid
  intp <- interp(as.vector(res$pars), as.vector(res$xis), as.vector(pots)) 
  xy <- expand.grid(intp$x, intp$y)
  z <- as.vector(intp$z)   	
  z[is.na(z)] <- max(na.omit(z))
  bg.var=NULL
  phylotype=NULL

  df <- data.frame(list(bg.var = xy[,1], phylotype = xy[,2], potential = z))
  p<-ggplot(df, aes(bg.var, phylotype, z = potential)) + geom_tile(aes(fill = potential)) + stat_contour(binwidth = 0.2) 
  p <- p + xlab(xlab.text) + ylab(ylab.text) + labs(title = title)
  p

}



#' Description: Potential Analysis
#'
#' \code{livpotential_ews} performs one-dimensional potential estimation derived from a uni-variate timeseries
#'
# Details:
#' see ref below
#' 
#  Arguments:
#'    @param x data vector
#'    @param std the standard deviation of the noise (defaults to 1, so then you use scaled potentials
#'    @param bw bandwidth for kernel estimation
#'    @param xi x values at which the potential is estimated
#'    @param weights optional weights in ksdensity (used by movpotentials).
#'    @param grid.size grid size
#'
# Returns:
#'   @return \code{livpotential} returns a list with the following elements:
#'   @return \item{xi}{the grid of points on which the potential is estimated}
#'   @return \item{pot}{the actual value of the potential}
#'   @return \item{minima}{the grid points at which the potential has minimum values}
#'   @return \item{maxima}{the grid points at which the potential has maximum values}
#'   @return \item{bw}{bandwidth of kernel used}
#'
#' @export
#'
#' @references Livina, VN, F Kwasniok, and TM Lenton, 2010. Potential analysis reveals changing number of climate states during the last 60 kyr . Climate of the Past, 6, 77-82.
#' 
#' Dakos, V., et al (2012)."Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data." \emph{PLoS ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010
#' @author Based on Matlab code from Egbert van Nes modified by Leo Lahti. Implemented in early warnings package by V. Dakos.
#' @seealso 
#' \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}}; \code{\link{bdstest_ews}}; \code{\link{sensitivity_ews}};\code{\link{surrogates_ews}}; \code{\link{ch_ews}};\code{\link{movpotential_ews}}
# ; \code{\link{timeVAR_ews}}; \code{\link{thresholdAR_ews}}
#' @examples 
#' data(foldbif)
#' res <- livpotential_ews(foldbif)
#' plot(res$xi, res$pot) 
#' @keywords early-warning

livpotential_ews <- function (x, std = 1, bw = -1, xi = NULL, weights = c(), grid.size = 200) {

  library(stats)
  x=data.frame(x)
  
  if (is.null(xi)) {
    xi <- seq(min(x), max(x), length = grid.size)
  } 

  if (bw < 0) {
   # following Silverman, B. W.: Density estimation for statistics and data analysis,Chapman and Hall, 1986.
    
   bw <- 1.06 * sapply(x,sd) / nrow(x)^(1/5);
  }

  # Density estimation
  # Matlab version: res <- ksdensity(x, xi, 'width', bw, 'npoints', 500, 'weights', weights);
  de <- density(ts(x), bw = bw, adjust = 1, kernel = "gaussian",
        weights = weights, window = kernel, n = length(xi), 
	from = min(x), to = max(x), cut = 3, na.rm = FALSE)

  # Estimated density
  f <- de$y

  # Final grid points and bandwidth
  xi <- de$x
  bw <- de$bw

  # Detect minima and maxima of the density (see Livina et al.)
  minima <- which(diff(sign(diff(c(0, f, 0)))) > 0)
  maxima <- which(diff(sign(diff(c(0, f, 0)))) < 0)

  # Compute potential
  U <- -log(f)*std^2/2

  list(xi = xi, pot = U, minima = minima, maxima = maxima, bw = bw)

}


#' Description: Moving Average Potential
#'
#' \code{movpotential} reconstructs a potential derived from data along a gradient of a given parameter
#Detail
#' the \code{movpotential} calculates the potential for values that correspond to a particular parameter. see ref below
#'
# Arguments:
#'  @param X a vector of the X observations of the state variable of interest
#'  @param param parameter values that correspond to the X observations
#'  @param sdwindow window for smoothing kernels (over the \code{param} axis)
#'  @param bw bandwidth used for smoothing kernels
#'  @param minparam minimum value of parameter on which to estimate potential
#'  @param maxparam maximum value of parameter on which to estimate potential
#'  @param npoints number of potentials
#'  @param thres threshold for local minima to be discarded
#'  @param std std
#'  @param grid.size number of evaluation points 
#'  @param cutoff the cuttof value to estimate minima and maxima in the potential
#'
#' Returns:
#'   @return A list with the following elements:
#'   @return \item{pars}{values of the covariate parameter as matrix}
#'   @return \item{xis}{values of the x as matrix}
#'   @return \item{pots}{smoothed potentials}
#'   @return \item{mins}{minima in the densities (-potentials; neglecting local optima)}
#'   @return \item{maxs}{maxima in densities (-potentials; neglecting local optima)}
#'   @return \item{plot}{an object that displays the potential estimated in 2D}
#' 
#' @export
#'
#' @references  Hirota, M., Holmgren, M., van Nes, E.H. & Scheffer, M. (2011). Global resilience of tropical forest and savanna to critical transitions. Science, 334, 232–235
#' @author Based on Matlab code from Egbert van Nes modified by Leo Lahti. Implemented in early warnings package by V. Dakos.
#' @seealso \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}}; \code{\link{bdstest_ews}}; \code{\link{sensitivity_ews}};\code{\link{surrogates_ews}}; \code{\link{ch_ews}}; \code{livpotential_ews}
# ; \code{\link{timeVAR_ews}}; \code{\link{thresholdAR_ews}}
#' @examples 
#'  X = c(rnorm(1000, mean = 0), rnorm(1000, mean = -2), rnorm(1000, mean = 2))
#'  param = seq(0,5,length=3000)
#'  res <- movpotential_ews(X, param, npoints = 100, thres = 0.003)
#' @keywords early-warning

movpotential_ews <- function (X, param, sdwindow = NULL, bw = -1, minparam = NULL, maxparam = NULL, npoints = 50, thres = 0.002, std = 1, grid.size = 200, cutoff=0.5) {

  # To debug:
  # set.seed(1324); X <- dat[, pt]; param <- annot[, varname]; npoints <- 100; thres <- 0.003; sdwindow <- NULL; bw = -1; minparam = NULL; maxparam = NULL; res <- movpotential(dat[, pt], annot[, varname], npoints = 100, thres = 0.003); std <- 1; cut.potential <- NULL

  if (is.null(minparam)) { 
    minparam <- min(param)
  }

  if (is.null(maxparam)) { 
    maxparam <- max(param)
  }

  if (is.null(sdwindow)) {
    sdwindow <- (maxparam - minparam) * 0.05
  }

  # Place evaluation points evenly across data range 
  xi <- seq(min(X), max(X), length = 200)

  # Determine step size
  step <- (maxparam - minparam) / npoints

  # Initialize
  xis <- pars <- pots <- matrix(0, nrow = npoints, ncol = length(xi))
  maxs <- mins <- matrix(0, nrow = npoints, ncol = length(xi))
  
  for (i in 1:npoints) {

    # print(i)
  
    # Increase the parameter at each step
    par <- minparam + (i - 0.5) * step

    # Check which elements in evaluation range (param) are within 2*sd of par
    param2 <- abs(par - param) / sdwindow
    #index <- which(param2 < 2)
    index <- 1:length(param2) # Use all. LL converted into this after getting into problems with occasionally sparse data.    

    weights <- exp(-param2[index]^2 / 2)
    weights <- weights/sum(weights) # LL needed to add normalization in the R implementation 16.5.2012

    # Calculate the potential
    tmp <- livpotential_ews(X[index], std, bw, xi, weights, grid.size)

    # Store variables
    pots[i, ] <- tmp$pot
     xis[i, ] <- tmp$xi 
    pars[i, ] <- par + rep(0, length(tmp$xi))

    # Mark the final minima and maxima for this evaluation point
    # (minima and maxima are for the density)    
    fpot <- exp(-2*tmp$pot/std^2) # backtransform to density distribution
    ops  <- find.optima(tmp$minima, tmp$maxima, fpot, thres)
    mins[i, ops$min] <- 1
    maxs[i, ops$max] <- 1

  }  

  res=list(pars = pars, xis = xis, pots = pots, mins = mins, maxs = maxs, std = std)
  p=PlotPotential(res, title = "Moving Average Potential", 'parameter', 'state variable', cutoff = cutoff)
  list(res=res,plot=p)
}



# Description: find.optima
#
# Detect optima from the potential
#
#  Arguments:
#    @param minima minima
#    @param maxima maxima
#    @param fpot potential
#    @param thres threshold
#
# Returns:
#   @return A list with the following elements:
#     min potential minima
#     max potential maxima
#
# @export
#
# @references See citation("TBA") 
# @author Leo Lahti \email{leo.lahti@@iki.fi}
# @examples #
#
# @keywords utilities

find.optima <- function (minima, maxima, fpot, thres) {

    # Remove minima and maxima that are too shallow
    delmini <- logical(length(minima))
    delmaxi <- logical(length(maxima))

    for (j in 1:length(maxima)) {

      # Calculate distance of this maximum to all minima
      s <- minima - maxima[[j]]

      # Set distances to deleted minima to zero
      s[delmini] <- 0

      # identify the closest minima
      i1 <- i2 <- NULL
      if (length(s)>0) {

        minima.spos <- minima[s > 0]
	minima.sneg <- minima[s < 0]

        if (length(minima.spos) > 0) {i1 <- min(minima.spos)}
        if (length(minima.sneg) > 0) {i2 <- max(minima.sneg)}

      } 
        
      # if no positive differences available, set it to same value with i2
      if (is.null(i1) && !is.null(i2)) {
         i1 <- i2
      }

      # if no negative differences available, set it to same value with i1
      if (is.null(i2) && !is.null(i1)) {
         i2 <- i1
      }

      # If a closest minimum exists, check differences and remove if difference is under threshold
      if (!is.null(i1)) {
  
        # Smallest difference between this maximum and the closest minima
        diff <- min( abs(fpot[i1] - fpot[maxima[[j]]]), 
	     	     abs(fpot[i2] - fpot[maxima[[j]]]))

        if (diff < thres) {

	  # If difference is below threshold, delete this maximum 
          delmaxi[[j]] <- 1

	  # Delete the larger of the two neighboring minima 
          if (fpot[[i1]] > fpot[[i2]]) {
            delmini[minima == i1] <- TRUE
	  } else {
            delmini[minima == i2] <- TRUE
	  }
        }   
      } else {
        # if both i1 and i2 are NULL, do nothing	 
      }

    }

    # Delete the shallow minima and maxima 
    if (length(minima) > 0 && sum(delmini)>0) {
      minima <- minima[-delmini]
    }
    if (length(maxima) > 0 && sum(delmaxi)>0) {
      maxima <- maxima[-delmaxi]
    }


    list(min = minima, max = maxima)
  
}
