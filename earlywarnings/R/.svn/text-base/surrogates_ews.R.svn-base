# Surrogates Early Warning Signals
# Author: Vasilis Dakos, January 3, 2012
	
# Load required packages
  #install.packages(c("lmtest","nortest","stats","som","Kendall","KernSmooth","e1071"), repos = c("http://R-Forge.R-project.org", "http://cran-mirror.cs.uu.nl/"), dep = TRUE)
	
surrogates_ews<-function(timeseries,indicator=c("ar1","sd","acf1","sk","kurt","cv","returnrate","densratio"),winsize=50,detrending=c("no","gaussian","linear","first-diff"),bandwidth=NULL,boots=100,logtransform=FALSE,interpolate=FALSE){
	
	require(lmtest)
	require(nortest)
	require(stats)
	require(som)
	require(Kendall)
	require(KernSmooth)
	require(moments)
	
#timeseries<-ts(timeseries) #strict data-types the input data as tseries object for use in later steps
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
			      			indicator = match.arg(indicator)
							if(indicator == "ar1"){
							indic<-apply(nMR,2,function(x){nAR1<-ar.ols(x,aic= FALSE, order.max=1,dmean=FALSE,intercept=FALSE)
							nAR1$ar})}
							else if(indicator == "sd"){
							indic<-apply(nMR,2,sd)}
							else if(indicator == "sk"){
							indic<-apply(nMR,2,skewness)}
							else if(indicator == "kurt"){
							indic<-apply(nMR,2,kurtosis)}
							else if(indicator == "acf1"){
							indic<-apply(nMR,2,function(x){nACF<-acf(x,lag.max = 1, type = c("correlation"),plot=FALSE)
							nACF$acf[2]})}
							else if(indicator == "returnrate"){
							indic<-apply(nMR,2,function(x){nACF<-acf(x,lag.max = 1, type = c("correlation"),plot=FALSE)
							1-nACF$acf[2]})}
							else if(indicator == "cv"){
							indic<-apply(nMR,2,function(x){sd(x)/mean(x)})}
							else if(indicator == "densratio"){
							indic<- apply(nMR,2,function(x){
							spectfft<-spec.ar(x,n.freq=omw,plot=FALSE,order=1)
							spectfft$spec
							spectfft$spec[low]/spectfft$spec[high]})}
				# Calculate trend statistics
				timevec<-seq(1,length(indic))
				Kt<-cor.test(timevec,indic,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
				Ktauestindorig<-Kt$estimate
				Ktaupindorig<-Kt$p.value

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

	Ktauestind<-numeric()
	Ktaupind<-numeric()

	for (jjj in 1:boots){
	x=arima.sim(n = length(nsmY), list(ar = c(ARMA$coef[1:ind[1]]), ma = c(ARMA$coef[(1+ind[1]):(ind[1]+ind[2]-1)])),sd=sqrt(ARMA$sigma2))

	## Rearrange data for indicator calculation
	nMR1<-matrix(data=NA,nrow=mw,ncol=omw)
	for (i in 1:omw){   
	Ytw<-x[i:(i+mw-1)]
	nMR1[,i]<-Ytw}
	
	# Estimate indicator
	indicator = match.arg(indicator)
	if(indicator == "ar1"){
	indic<-apply(nMR1,2,function(x){nAR1<-ar.ols(x,aic= FALSE, order.max=1,dmean=FALSE,intercept=FALSE)
	nAR1$ar})}
	else if(indicator == "sd"){
	indic<-apply(nMR1,2,sd)}
	else if(indicator == "sk"){
	indic<-apply(nMR1,2,skewness)}
	else if(indicator == "kurt"){
	indic<-apply(nMR1,2,kurtosis)}
	else if(indicator == "acf1"){
	indic<-apply(nMR1,2,function(x){nACF<-acf(x,lag.max = 1, type = c("correlation"),plot=FALSE)
	nACF$acf[2]})}
	else if(indicator == "returnrate"){
	indic<-apply(nMR1,2,function(x){nACF<-acf(x,lag.max = 1, type = c("correlation"),plot=FALSE)
	1-nACF$acf[2]})}
	else if(indicator == "cv"){
	indic<-apply(nMR1,2,function(x){sd(x)/mean(x)})}
	else if(indicator == "densratio"){
	indic<- apply(nMR1,2,function(x){
	spectfft<-spec.ar(x,n.freq=omw,plot=FALSE,order=1)
	spectfft$spec
	spectfft$spec[low]/spectfft$spec[high]})
	}
	
	# Calculate trend statistics
	timevec<-seq(1,length(indic))
	Kt<-cor.test(timevec,indic,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
	Ktauestind[jjj]<-Kt$estimate
	Ktaupind[jjj]<-Kt$p.value
}

	# Estimate probability of false positive
	q<-sort(Ktauestind,na.last=NA)
	Kpos<-max(which(Ktauestindorig>q),na.rm=TRUE)
	p<-(boots+1-Kpos)/boots
	print(paste('significance p = ',p,' estimated from ',boots,' surrogate ARMA timeseries'))

	# Plotting
	dev.new()
	par(font.main=10,mar=(c(4.6,4,0.5,2)+0.2),oma=c(0.5,1,2,1))
	hist(Ktauestind,freq=TRUE,nclass=20,xlim=c(-1,1),col="blue",main=NULL,xlab="Surrogate Kendall tau estimates",ylab="occurrence",ylim=c(0,boots))
	abline(v=q[0.05*boots],col="black",lwd=1)
	abline(v=q[0.95*boots],col="black",lwd=1)
	points(Ktauestindorig,0,pch=21, bg="black", col = "black", cex=1)
	mtext(paste("Indicator ",toupper(indicator)),side=3,line=0.2, outer=TRUE)

	# Output
	out<-data.frame(Ktauestindorig,Ktaupindorig,Ktauestind,Ktaupind,p)
	colnames(out)<-c("Kendall tau estimate original","Kendall tau p-value original","Kendall tau estimate surrogates","Kendall tau p-value surrogates","significance p")
	return(out)
}