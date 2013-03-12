#' Description: Quick Detection Analysis for Generic Early Warning Signals
#'
#' \code{qda_ews} is used to estimate autocorrelation, variance within rolling windows along a timeseries, test the significance of their trends, and reconstruct the potential landscape of the timeseries
#'
# Details:
#' see ref below
#'
# Arguments:
#'    @param timeseries a numeric vector of the observed univariate timeseries values or a numeric matrix where the first column represents the time index and the second the observed timeseries values. Use vectors/matrices with headings.
#'    @param param values corresponding to observations in timeseries
#'    @param winsize is the size of the rolling window expressed as percentage of the timeseries length (must be numeric between 0 and 100). Default is 50\%.
#'    @param detrending the timeseries can be detrended/filtered prior to analysis. There are four options: \code{gaussian} filtering, \code{linear} detrending and \code{first-differencing}. Default is \code{no} detrending.
#'    @param bandwidth is the bandwidth used for the Gaussian kernel when gaussian filtering is applied. It is expressed as percentage of the timeseries length (must be numeric between 0 and 100). Alternatively it can be given by the bandwidth selector \code{\link{bw.nrd0}} (Default).
#    @param incrwinsize increments the rolling window size (must be numeric between 0 and 100). Default is 25.
#'    @param boots the number of surrogate data to generate from fitting an ARMA(p,q) model. Default is 100.
#'    @param s_level significance level. Default is 0.05.
#    @param incrbandwidth is the size to increment the bandwidth used for the Gaussian kernel when gaussian filtering is applied. It is expressed as percentage of the timeseries length (must be numeric between 0 and 100). Default is 20.
#'    @param cutoff the cutoff value to visualize the potential landscape
#'    @param detection.threshold detection threshold for potential minima
#'    @param grid.size grid size (for potential analysis)
#'    @param logtransform logical. If TRUE data are logtransformed prior to analysis as log(X+1). Default is FALSE.
#'    @param interpolate logical. If TRUE linear interpolation is applied to produce a timeseries of equal length as the original. Default is FALSE (assumes there are no gaps in the timeseries).
#' 
# Returns:
#' \code{qda_ews} returns three plots. The first plot contains the original data, the detrending/filtering applied and the residuals (if selected), autocorrelation and variance. For each statistic trends are estimated by the nonparametric Kendall tau correlation.  The second plot, returns a histogram of the distributions of the Kendall trend statistic for autocorrelation and variance estimated on the surrogated data. Vertical lines represent the level of significance, whereas the black dots the actual trend found in the time series. The third plot is the reconstructed potential landscape in 2D.
#'  
#' @export
#' 
#' @author Vasilis Dakos, Leo Lahti, March 1, 2013 \email{vasilis.dakos@@gmail.com}
#' @references
#' Dakos, V., et al (2012)."Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data." \emph{PLoS ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010 
#' @seealso \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}}; \code{\link{bdstest_ews}}; \code{\link{sensitivity_ews}}; \code{\link{surrogates_ews}}; \code{\link{ch_ews}}; \code{\link{movpotential_ews}}; \code{\link{livpotential_ews}}; 
#'
#' @examples 
#' data(foldbif)
#' out <- qda_ews(foldbif, param = NULL, winsize = 50, detrending="gaussian", bandwidth=NULL, boots = 50, s_level = 0.05, cutoff=0.05, detection.threshold = 0.002, grid.size = 50, logtransform=FALSE, interpolate=FALSE)
#' @keywords early-warning

qda_ews <- function(timeseries, param = NULL, winsize = 50, detrending=c("no","gaussian","linear","first-diff"), bandwidth=NULL, boots = 100, s_level = 0.05, cutoff=0.05, detection.threshold = 0.002, grid.size = 50, logtransform=FALSE, interpolate=FALSE){
  
  timeseries <- data.matrix(timeseries)
    message("Indicator trend analysis")
    generic_RShiny(timeseries,winsize,detrending,bandwidth,logtransform,interpolate,AR_n=FALSE,powerspectrum=FALSE)
      
    message("Trend significance analysis")
    s <- surrogates_RShiny(timeseries,winsize,detrending,bandwidth,boots,s_level,logtransform,interpolate)
    print(s)

    message("Potential analysis")
    p <- movpotential_RShiny(timeseries, param, detection.threshold = detection.threshold, grid.size = grid.size, plot.cutoff = cutoff)
    print(p)
    
    
    #   message("Sensitivity of trends")
    #   sens <- sensitivity_RShiny(timeseries,winsizerange=c(25,75),incrwinsize,detrending=detrending, bandwidthrange=c(5,100),incrbandwidth,logtransform=FALSE,interpolate=FALSE)

}
  
# generic_Rshiny for estimating only AR1 and Variance in moving windows with various options for pretreating the data
# 26 Feb 2013

generic_RShiny<-function(timeseries, winsize = 50, detrending = c("no", "gaussian", "linear", "first-diff"), bandwidth, logtransform, interpolate, AR_n = FALSE, powerspectrum = FALSE){
	
# 	require(lmtest)
# 	require(nortest)
# 	require(stats)
# 	require(som)
# 	require(Kendall)
# 	require(KernSmooth)
# 	require(moments)
    
    #timeseries<-ts(timeseries)
    timeseries<-data.matrix(timeseries) #strict data-types the input data as tseries object for use in later steps
    if (ncol(timeseries) == 1){
		Y=timeseries
		timeindex=1:dim(timeseries)[1]
    }else if(dim(timeseries)[2]==2){
		Y<-timeseries[,2]
		timeindex<-timeseries[,1]
    }else{
		warning("not right format of timeseries input")
    }
    #return(timeindex)
	# Interpolation
	if (interpolate){
		YY<-approx(timeindex,Y,n=length(Y),method="linear")
		Y<-YY$y
    }else{
		Y<-Y}
    
	# Log-transformation
	if (logtransform){
		Y<-log(Y+1)}
	
	# Detrending
	detrending<-match.arg(detrending)
	if (detrending=="gaussian"){
		if (is.null(bandwidth)){
			bw<-round(bw.nrd0(timeindex))}else{
                bw<-round(length(Y)*bandwidth/100)}
        smYY<-ksmooth(timeindex,Y,kernel="normal",bandwidth=bw,range.x=range(timeindex),x.points=timeindex)
        if (timeindex[1]>timeindex[length(timeindex)]){
            nsmY<-Y-rev(smYY$y)
            smY<-rev(smYY$y)
        }else{  	nsmY<-Y-smYY$y
            smY<-smYY$y}
	}else if(detrending=="linear"){
		nsmY<-resid(lm(Y~timeindex))
		smY<-fitted(lm(Y~timeindex))
	}else if(detrending=="first-diff"){
		nsmY<-diff(Y)
		timeindexdiff<-timeindex[1:(length(timeindex)-1)]
	}else if(detrending=="no"){
		smY<-Y
		nsmY<-Y
	}
    
    
	# Rearrange data for indicator calculation
	mw<-round(length(Y)*winsize/100)
	omw<-length(nsmY)-mw+1 ##number of moving windows
	low<-6
	high<-omw
	nMR<-matrix(data=NA,nrow=mw,ncol=omw)
	x1<-1:mw
    for (i in 1:omw){
		Ytw<-nsmY[i:(i+mw-1)]
		nMR[,i]<-Ytw}
    
	# Calculate indicators
	nARR<-numeric()
	nSD<-numeric()
    
	nSD<-apply(nMR, 2, sd, na.rm = TRUE)
	for (i in 1:ncol(nMR)){
		nYR<-ar.ols(nMR[,i],aic= FALSE, order.max=1, dmean=FALSE, 		intercept=FALSE)
		nARR[i]<-nYR$ar
	}
	
    nVAR=sqrt(nSD)
    
	# Estimate Kendall trend statistic for indicators
	timevec<-seq(1,length(nARR))
	KtAR<-cor.test(timevec,nARR,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
	KtVAR<-cor.test(timevec,nVAR,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
    
	# Plotting
	# Generic Early-Warnings
	#dev.new()
	par(mar=(c(1,2,0.5,2)+0),oma=c(2,2,2,2),mfrow=c(4,1))
	plot(timeindex,Y,type="l",ylab="",xlab="",xaxt="n",lwd=2,las=1,xlim=c(timeindex[1],timeindex[length(timeindex)]))
	legend("bottomleft","data",,bty = "n")
	if(detrending=="gaussian"){
		lines(timeindex,smY,type="l",ylab="",xlab="",xaxt="n",lwd=2,col=2,las=1,xlim=c(timeindex[1],timeindex[length(timeindex)]))
	}
	if(detrending=="no"){
		plot(c(0,1),c(0,1),ylab="",xlab="",yaxt="n",xaxt="n",type="n",las=1)
		text(0.5,0.5,"no detrending - no residuals")
    }else if (detrending=="first-diff"){
		limit<-max(c(max(abs(nsmY))))
		plot(timeindexdiff,nsmY,ylab="",xlab="",type="l",xaxt="n",lwd=2,las=1,ylim=c(-	limit,limit),xlim=c(timeindexdiff[1],timeindexdiff[length(timeindexdiff)]))
		legend("bottomleft","first-differenced",bty = "n")		}else{
            limit<-max(c(max(abs(nsmY))))
            plot(timeindex,nsmY,ylab="",xlab="",type="h",xaxt="n",las=1,lwd=2,ylim=c(-	limit,limit),xlim=c(timeindex[1],timeindex[length(timeindex)]))
            legend("bottomleft","residuals",bty = "n")}
	plot(timeindex[mw:length(nsmY)],nARR,ylab="",xlab="",type="l",xaxt="n",col="green",lwd=2,las=1,xlim=c(timeindex[1],timeindex[length(timeindex)])) #3
	legend("bottomright",paste("trend ",round(KtAR$estimate,digits=3)),bty = "n")
	legend("bottomleft","autocorrelation",bty = "n")
	plot(timeindex[mw:length(nsmY)],nVAR,ylab="",xlab="",type="l",col="blue", lwd=2, las=1,xlim=c(timeindex[1],timeindex[length(timeindex)]))
	legend("bottomright",paste("trend ",round(KtVAR$estimate,digits=3)),bty = "n")
	legend("bottomleft","variance",bty = "n")
	mtext("time",side=1,line=2,cex=0.8)
	mtext("Generic Early-Warnings: Autocorrelation - Variance",side=3,line=0.2, outer=TRUE)#outer=TRUE print on the outer margin
	
    # Output
    out<-data.frame(timeindex[mw:length(nsmY)],nARR,nSD)
    colnames(out)<-c("timeindex","ar1","sd")
    #return(out)
    
}

# surrogates_Rshiny for estimating significance of trends for variance and autocorrelation
# 6 March 2013

surrogates_RShiny<-function(timeseries,winsize=50,detrending=c("no","gaussian","linear","first-diff"),bandwidth=NULL,boots=100,s_level=0.05,logtransform=FALSE,interpolate=FALSE){
    
    timeseries<-data.matrix(timeseries)
    if (dim(timeseries)[2]==1){
		Y=timeseries
		timeindex=1:dim(timeseries)[1]
    }else if(dim(timeseries)[2]==2){
		Y<-timeseries[,2]
		timeindex<-timeseries[,1]
    }else{
		warning("not right format of timeseries input")
    }
    
	# Interpolation
	if (interpolate){
		YY<-approx(timeindex,Y,n=length(Y),method="linear")
		Y<-YY$y
    }else{
		Y<-Y}
    
	# Log-transformation
	if (logtransform){
		Y<-log(Y+1)}
	
	# Detrending
	detrending<-match.arg(detrending)
	if (detrending=="gaussian"){
		if (is.null(bandwidth)){
			bw<-round(bw.nrd0(timeindex))}else{
                bw<-round(length(Y)*bandwidth)/100}
		smYY<-ksmooth(timeindex,Y,kernel=c("normal"), bandwidth=bw, range.x=range(timeindex),n.points=length(timeindex))
		nsmY<-Y-smYY$y
		smY<-smYY$y
	}else if(detrending=="linear"){
		nsmY<-resid(lm(Y~timeindex))
		smY<-fitted(lm(Y~timeindex))
	}else if(detrending=="first-diff"){
		nsmY<-diff(Y)
		timeindexdiff<-timeindex[1:(length(timeindex)-1)]
	}else if(detrending=="no"){
		smY<-Y
		nsmY<-Y
	}
	
	
	# Rearrange data for indicator calculation
    mw<-round(length(Y)*winsize)/100
    omw<-length(nsmY)-mw+1
    low<-6
    high<-omw
    nMR<-matrix(data=NA,nrow=mw,ncol=omw)
    for (i in 1:omw){
        Ytw<-nsmY[i:(i+mw-1)]
        nMR[,i]<-Ytw}
    # Estimate indicator
    
    indic_ar1<-apply(nMR,2,function(x){nAR1<-ar.ols(x,aic= FALSE, order.max=1,dmean=FALSE,intercept=FALSE)
        nAR1$ar})
    
    indic_var<-apply(nMR,2,var)
    
    # Calculate trend statistics
    timevec<-seq(1,length(indic_ar1))
    Kt_ar1<-cor.test(timevec,indic_ar1,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
    Ktauestind_ar1orig<-Kt_ar1$estimate
    
    Kt_var<-cor.test(timevec,indic_var,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
    Ktauestind_varorig<-Kt_var$estimate
    
	# Fit ARMA model based on AIC
	arma=matrix(,4,5)
	for (ij in 1:4){
        for (jj in 0:4){
            ARMA<-arima(nsmY, order = c(ij,0,jj),include.mean = FALSE)
            arma[ij,jj+1]=ARMA$aic
            print(paste("AR","MA", "AIC"),quote=FALSE)
            print(paste(ij,jj,ARMA$aic),zero.print=".",quote=FALSE)
		}
    }
    
	# Simulate ARMA(p,q) model fitted on residuals
	ind=which(arma==min(arma),arr.ind=TRUE)
	ARMA<-arima(nsmY, order = c(ind[1],0,ind[2]-1),include.mean = FALSE)
    
    Ktauestind_ar1 <- numeric()
    Ktauestind_var <- numeric()
    
	for (jjj in 1:boots){
        x=arima.sim(n = length(nsmY), list(ar = c(ARMA$coef[1:ind[1]]), ma = c(ARMA$coef[(1+ind[1]):(ind[1]+ind[2]-1)])),sd=sqrt(ARMA$sigma2))
        
        ## Rearrange data for indicator calculation
        nMR1<-matrix(data=NA,nrow=mw,ncol=omw)
        for (i in 1:omw){
            Ytw<-x[i:(i+mw-1)]
            nMR1[,i]<-Ytw}
        
        # Estimate indicator
        
        indic_ar1<-apply(nMR1,2,function(x){nAR1<-ar.ols(x,aic= FALSE, order.max=1,dmean=FALSE,intercept=FALSE)
            nAR1$ar})
        
        indic_var<-apply(nMR1,2,var)
        
        # Calculate trend statistics
        timevec<-seq(1,length(indic_ar1))
        Kt_ar1<-cor.test(timevec,indic_ar1,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
        Ktauestind_ar1[jjj]<-Kt_ar1$estimate
        
        Kt_var<-cor.test(timevec,indic_var,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
        Ktauestind_var[jjj]<-Kt_var$estimate
        
    }
    
	# Estimate probability of false positive
	q_ar1<-sort(Ktauestind_ar1,na.last=NA)
	Kpos_ar1<-max(which(Ktauestind_ar1orig>q_ar1),na.rm=TRUE)
	p<-(boots+1-Kpos_ar1)/boots
	print(paste('significance autocorrelation p = ',p,' estimated from ',boots,' surrogate ARMA timeseries'))
    
    q_var<-sort(Ktauestind_var,na.last=NA)
    Kpos_var<-max(which(Ktauestind_varorig>q_var),na.rm=TRUE)
    p<-(boots+1-Kpos_var)/boots
    print(paste('significance variance p = ',p,' estimated from ',boots,' surrogate ARMA timeseries'))
    
	# Plotting
    layout(matrix(1:2,1,2))
    par(font.main=10,mar=(c(4.6,3.5,0.5,2)+0.2),mgp=c(2,1,0),oma=c(0.5,0.5,2,0),cex.axis=0.8,cex.lab=0.8,cex.main=0.8)
    hist(Ktauestind_ar1,freq=TRUE,nclass=20,xlim=c(-1,1),col="green",main=NULL,xlab="Surrogate trend estimates",ylab="occurrence")#,ylim=c(0,boots))
    abline(v=q_ar1[s_level*boots],col="red",lwd=2)
    abline(v=q_ar1[(1-s_level)*boots],col="red",lwd=2)
    points(Ktauestind_ar1orig,0,pch=21, bg="black", col = "black", cex=4)
    title("Autocorrelation",cex.main=1.3)
    
    hist(Ktauestind_var,freq=TRUE,nclass=20,xlim=c(-1,1),col="blue",main=NULL,xlab="Surrogate trend estimates",ylab="occurrence")#,ylim=c(0,boots))
    abline(v=q_var[s_level*boots],col="red",lwd=2)
    abline(v=q_var[(1-s_level)*boots],col="red",lwd=2)
    points(Ktauestind_varorig,0,pch=21, bg="black", col = "black", cex=4)
    title("Variance" ,cex.main=1.3)
}

# Description: Plot Potential

PlotPotential <- function (res, title = "", xlab.text, ylab.text, cutoff = 0.5) {
  
  cut.potential <- max(apply(res$pots, 1, min)) + cutoff*abs(max(apply(res$pots, 1, min))) # Ensure all minima are visualized
  pots <- res$pots
  pots[pots > cut.potential] <- cut.potential
  
  # Static contour
  # Interpolate potential grid
  intp <- tgp::interp.loess(as.vector(res$pars), as.vector(res$xis), as.vector(pots))
  xy <- expand.grid(intp$x, intp$y)
  z <- as.vector(intp$z)
  z[is.na(z)] <- max(na.omit(z))
  
  potential <- NULL
  df <- data.frame(list(bg.var = xy[,1], phylotype = xy[,2], potential = z))
  bg.var <- NULL
  phylotype <- NULL
  
  p <- ggplot2::ggplot(df, aes(bg.var, phylotype, z = potential)) + geom_tile(aes(fill = potential)) + stat_contour(binwidth = 0.2)
  
  p <- p + xlab(xlab.text) + ylab(ylab.text) + labs(title = title)
  
  p
  
}



# Description: Potential Analysis
 
livpotential_ews <- function (x, std = 1, bw = -1, xi = NULL, weights = c(), grid.size = 200) {
  
  x <- data.frame(x)
  
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



movpotential_RShiny <- function (X, param = NULL, bw = -1, detection.threshold = 0.002, std = 1, grid.size = 50, plot.cutoff = 0.5) {
  
  if (is.null(param)) {
    param <- seq(1, nrow(X), 1)
  }
  minparam <- min(param)    
  maxparam <- max(param)
  
  # Determine step size
  sdwindow <- step <- (maxparam - minparam) / grid.size
  
  # Place evaluation points evenly across data range 
  xi <- seq(min(X), max(X), length = grid.size)
  
  # Initialize
  xis <- pars <- pots <- matrix(0, nrow = grid.size, ncol = length(xi))
  maxs <- mins <- matrix(0, nrow = grid.size, ncol = length(xi))
  
  for (i in 1:grid.size) {
    
    # Increase the parameter at each step
    par <- minparam + (i - 0.5) * step
    
    # Check which elements in evaluation range (param) are within 2*sd of par
    weights <- exp(-.5*(abs(par - param) / sdwindow)^2)
    weights <- weights/sum(weights) # LL added normalization in the R implementation 16.5.2012
    
    # Calculate the potential
    tmp <- livpotential_ews(X, std, bw, xi, weights, grid.size)
    
    # Store variables
    pots[i, ] <- tmp$pot
    xis[i, ] <- tmp$xi 
    pars[i, ] <- par + rep(0, length(tmp$xi))
    
    # Mark the final minima and maxima for this evaluation point
    # (minima and maxima are for the density)    
    fpot <- exp(-2*tmp$pot/std^2) # backtransform to density distribution
    ops  <- find.optima(tmp$minima, tmp$maxima, fpot, detection.threshold)
    mins[i, ops$min] <- 1
    maxs[i, ops$max] <- 1
    
  }  
  
  res <- list(pars = pars, xis = xis, pots = pots, mins = mins, maxs = maxs, std = std)
  p <- PlotPotential(res, title = "Moving Average Potential", 'parameter/time', 'state variable', cutoff = plot.cutoff)
  
  list(res = res, plot = p)
  
}



# Description: find.optima

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
