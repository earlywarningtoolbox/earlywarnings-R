# BDS test Early Warning Signals
# Author: Stephen R Carpenter, 22 Oct 2011
# Modified by: Vasilis Dakos, January 1, 2012

#install.packages(c("tseries","quadprog"), repos = c("http://R-Forge.R-project.org", "http://cran-mirror.cs.uu.nl/"), dep = TRUE)

BDSboot <- function(X,varname,nboot,epsvec,emb) { # begin function

require(tseries)
	
StdEpsAll <- X  # name of variable for BDS
neps <- length(epsvec)
# Compute and print BDS test
print('***********************************************',quote=FALSE)
print(c('BDS test for ',varname),quote=FALSE)
print(c('Embedding dimension = ',emb),quote=FALSE)
BDS.data <- bds.test(StdEpsAll,m=emb,epsvec)
print('BDS statistics for Nominal Data at each Epsilon',quote=FALSE)
print(round(BDS.data$statistic,3))
print('P value based on standard normal',quote=FALSE)
print(round(BDS.data$p.value,3))
# Bootstrap the BDS test
nobs <- length(StdEpsAll)
bootmat <- matrix(0,nrow=emb-1,ncol=neps)  # matrix to count extreme BDS values
for(i in 1:nboot) { # start bootstrap loop
 epsboot <- sample(StdEpsAll,nobs,replace=TRUE)
 BDS.boot <- bds.test(epsboot,m=emb,epsvec)
 for(im in 1:(emb-1)) {  # loop over embedding dimensions
   bootvec <- BDS.boot$statistic[im,]
   N.above <- ifelse(bootvec>BDS.data$statistic[im,],1,0)
   bootmat[im,] <- bootmat[im,]+N.above
   }
 # Report progress: if hash is removed from the next two lines, the program
 # will report each time an iteration is completed
 # cat('iteration = ',i,' of ',nboot,'\n') # 
 # flush.console()
 } # end bootstrap loop
print(' ',quote=FALSE)
print(c('Bootstrap P estimates for ',varname),quote=FALSE)
print(c('Bootstrap iterations = ',nboot),quote=FALSE)
Pboot <- bootmat/nboot
for(im in 1:(emb-1)) {
  print(c('For embedding dimension =',im+1),quote=FALSE)
  print(c('For epsilon = ',round(epsvec,3),'bootstrap P = '),quote=FALSE)
  print(Pboot[im,])
  }
print('**********************************************************',quote=FALSE)
} # end function

bdstest_ews<-function(timeseries,ARMAoptim=TRUE,ARMAorder=c(1,0),GARCHorder=c(0,1),embdim=3,epsilon=c(0.5,0.75,1),boots=1000,logtransform=FALSE,interpolate=FALSE){
	
	require(quadprog)
	
	timeseries<-ts(timeseries) #strict data-types the input data as tseries object for use in later steps
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
		
	# Detrend the data
	Eps1 <- diff(Y)
		
	# Define BDS parameters
	nboot <- boots
	emb <- embdim # embedding dimension
	eps.sd <- sd(Eps1)
	epsvec <- round(eps.sd* epsilon,6)

	# Run BDS with bootstrapping
	BDSboot(Eps1,c('Detrended data'),nboot,epsvec,emb)
	
	# Fit ARMA model based on AIC
	if (ARMAoptim==TRUE){
		arma=matrix(,4,5)
	for (ij in 1:4){
	for (jj in 0:4){
		ARMA<-arima(Y, order = c(ij,0,jj),method="ML",include.mean = TRUE)
		arma[ij,jj+1]=ARMA$aic	
		ind=which(arma==min(arma),arr.ind=TRUE)
		armafit<-arima(Y, order = c(ind[1],0,ind[2]-1),include.mean = TRUE)	
		print('ARMA model',quote=FALSE)
		print(armafit,digits=4)
		Eps2 <- armafit$residuals#[2:armafit$n.used]
		}
		}
	}else{	
	armafit <- arima(Y, order = c(ARMAorder[1],0,ARMAorder[2]),include.mean = TRUE)
	print('ARMA model',quote=FALSE)
	print(armafit,digits=4)
	Eps2 <- armafit$residuals#[2:armafit$n.used]
	}

	# Define BDS parameters
	nboot <- boots
	emb <- embdim # embedding dimension
	eps.sd <- sd(Eps2)
	epsvec <- round(eps.sd*epsilon,6)

	# Run BDS with bootstrapping
	BDSboot(Eps2,c('ARMA model residuals'),nboot,epsvec,emb)

	# Fit GARCH(0,1) model to detrended data
	Gfit <- garch(Y,order=c(GARCHorder[1],GARCHorder[2]))
	print('GARCH(0,1) model fit to detrended data',quote=FALSE)
	print(Gfit,digits=4)

	Eps3 <- Gfit$residuals[2:length(Y)]

	# Define BDS parameters
	nboot <- boots
	emb <- embdim # embedding dimension
	eps.sd <- sd(Eps3)
	epsvec <- round(eps.sd*epsilon,6)

	# Run BDS with bootstrapping
	BDSboot(Eps3,c('GARCH(0,1) model residuals'),nboot,epsvec,emb)
	
	# Plot the data
	dev.new()
	par(fig=c(0,0.5,0.5,1),mar=c(3, 4, 3, 2),cex.axis=0.8,cex.lab=1,mfrow=c(2,2),mgp=c(1.5,0.5,0),oma=c(1,1,2,1))
	plot(timeindex,Y,type='l',col='black',lwd=1.7,xlab='time',ylab='data')
	par(fig=c(0.5,1,0.8,0.95),mar=c(0, 4, 0, 2),new=TRUE)
	plot(timeindex[1:(length(timeindex)-1)],Eps1,type="l",col="red",lwd=1,xlab="",ylab="",xaxt="n")
	legend("topright","first-diff",bty="n",cex=0.8)
	par(fig=c(0.5,1,0.65,0.8), new=TRUE)
	plot(timeindex[1:(length(timeindex))],Eps2,type="l",xlab="",ylab="residuals",xaxt="n",col="green",lwd=1)
	legend("topright","AR",bty="n",cex=0.8)
	par(fig=c(0.5,1,0.5,0.65), new=TRUE)
	plot(timeindex[1:(length(timeindex)-1)],Eps3,type="l",col="blue",xlab="time",ylab="",lwd=1)
	legend("topright","GARCH",bty="n",cex=0.8)
	mtext("time",side=1,outer=FALSE,line=1.4,cex=0.8)
	# Time series diagnostics
	par(fig=c(0,0.5,0,0.4),mar=c(3, 4, 0, 2), new=TRUE)
	acf(Y,lag.max=25,main="")
	par(fig=c(0.5,1,0,0.4), new=TRUE)
	pacf(Y,lag.max=25,main="")
	mtext("BDS_test Diagnostics",side=3,line=0.2, outer=TRUE)
	}