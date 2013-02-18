# Sensitivity Early Warning Signals
# Author: Vasilis Dakos, January 4, 2012
	
# Load required packages
  install.packages("fields", repos = c("http://R-Forge.R-project.org", "http://cran-mirror.cs.uu.nl/"), dep = TRUE)

	
sensitivity_ews<-function(timeseries,indicator=c("ar1","sd","acf1","sk","kurt","cv","returnrate","densratio"),winsizerange=c(25,75),incrwinsize=25,detrending=c("no","gaussian","linear","first-diff"),bandwidthrange=c(5,100),incrbandwidth=20,logtransform=FALSE,interpolate=FALSE){
	
	require(lmtest)
	require(nortest)
	require(stats)
	require(som)
	require(Kendall)
	require(KernSmooth)
	require(moments)
	require(fields)
	
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

	# Determine the step increases in rolling windowsize
	incrtw<-incrwinsize
	tw<-seq(floor(winsizerange[1]*length(Y)/100),floor(winsizerange[2]*length(Y)/100),by=incrtw)
	twcol<-length(tw)
	low<-6
		
		# Detrending	
	detrending<-match.arg(detrending)	
	if (detrending=="gaussian"){
		incrbw<-incrbandwidth
		width<-seq(floor(bandwidthrange[1]*length(Y)/100),floor(bandwidthrange[2]*length(Y)/100),by=incrbw)
		bwrow<-length(width)
		# Create matrix to store Kendall trend statistics
		Ktauestind<-matrix(,bwrow,twcol)
		Ktaupind<-matrix(,bwrow,twcol)
			# Estimation
			for (wi in 1:(length(width))){
					width1<-width[wi]
					smYY<-ksmooth(timeindex,Y,kernel=c("normal"), bandwidth=width1, range.x=range(timeindex),n.points=length(timeindex))
					nsmY<-Y-smYY$y
				for (ti in 1:length(tw)){	
					tw1<-tw[ti]
					# Rearrange data for indicator calculation
					omw1<-length(nsmY)-tw1+1 ##number of overlapping moving windows
					high<-omw1
					nMR1<-matrix(data=NA,nrow=tw1,ncol=omw1)
						for (i in 1:omw1){
				   			Ytw<-nsmY[i:(i+tw1-1)]
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
							spectfft<-spec.ar(x,n.freq=omw1,plot=FALSE,order=1)
							spectfft$spec
							spectfft$spec[low]/spectfft$spec[high]})
							}
				# Calculate trend statistics
				timevec<-seq(1,length(indic))
				Kt<-cor.test(timevec,indic,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
				Ktauestind[wi,ti]<-Kt$estimate
				Ktaupind[wi,ti]<-Kt$p.value
				}
				}				
		# Plot
		layout(matrix(1:4,2,2))
		par(font.main=10,mar=(c(4.6,4,0.5,2)+0.2),mgp=c(2,1,0),oma=c(0.5,1,2,1))
		image.plot(width,tw,Ktauestind,zlim=c(-1,1),xlab="filtering bandwidth",ylab="rolling window size",main="Kendall tau estimate",cex.main=0.8,log="y",nlevel=20,col=rainbow(20))
		ind=which(Ktauestind==max(Ktauestind),arr.ind=TRUE)
		lines(width[ind[1]],tw[ind[2]],type="p",cex=1.2,pch=17,col=1)
		hist(Ktauestind,breaks=12,col="lightblue",main=NULL, xlab="Kendall tau estimate", ylab="occurence",border="black",xlim=c(-1,1))
		image.plot(width,tw,Ktaupind,zlim=c(0,max(Ktaupind)),xlab="filtering bandwidth",ylab="rolling window size",main="Kendall tau p-value",log="y",cex.main=0.8,nlevel=20,col=rainbow(20))
		lines(width[ind[1]],tw[ind[2]],type="p",cex=1.2,pch=17,col=1)
		hist(Ktaupind,breaks=12,col="yellow",main=NULL, xlab="Kendall tau p-value", ylab="occurence", border="black",xlim=c(0,max(Ktaupind)))
		mtext(paste("Indicator ",toupper(indicator)),side=3,line=0.2, outer=TRUE)
		
		# Output
		out<-data.frame(Ktauestind)
		colnames(out)<-tw
    rownames(out)<-width
		return(out)		

	}else if(detrending=="linear"){
		nsmY<-resid(lm(Y~timeindex))
	}else if(detrending=="first-diff"){
		nsmY<-diff(Y)
	}else if(detrending=="no"){
		nsmY<-Y
	}
	
	# Create matrix to store Kendall trend statistics
	Ktauestind<-matrix(,twcol,1)
	Ktaupind<-matrix(,twcol,1)
		
for (ti in 1:length(tw)){	
					tw1<-tw[ti]
					# Rearrange data for indicator calculation
					omw1<-length(nsmY)-tw1+1 ##number of overlapping moving windows
					high=omw1
					nMR1<-matrix(data=NA,nrow=tw1,ncol=omw1)
						for (i in 1:omw1){
				   			Ytw<-nsmY[i:(i+tw1-1)]
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
							spectfft<-spec.ar(x,n.freq=omw1,plot=FALSE,order=1)
							spectfft$spec
							spectfft$spec[low]/spectfft$spec[high]})}
				# Calculate trend statistics
				timevec<-seq(1,length(indic))
				Kt<-cor.test(timevec,indic,alternative=c("two.sided"),method=c("kendall"),conf.level=0.95)
				Ktauestind[ti]<-Kt$estimate
				Ktaupind[ti]<-Kt$p.value
				}
		# Plotting
		dev.new()
		layout(matrix(1:4,2,2))
		par(font.main=10,mar=(c(4.6,4,0.5,2)+0.2),mgp=c(2,1,0),oma=c(0.5,1,2,1))
		plot(tw,Ktauestind,type="l",ylim=c(-1,1),log="x",xlab="rolling window size",ylab="Kendall tau estimate")
		hist(Ktauestind,breaks=12,col="lightblue", xlab="Kendall tau estimate", ylab="occurence",border="black",xlim=c(-1,1),main=NULL)
		plot(tw,Ktaupind,type="l",xlab="rolling window size",log="x",ylab="Kendall tau p-value",ylim=c(0,max(Ktaupind)))
		hist(Ktaupind,breaks=12,col="yellow", xlab="Kendall tau p-value",main=NULL, ylab="occurence",border="black",xlim=c(0,max(Ktaupind)))
		mtext(paste("Indicator ",toupper(indicator)),side=3,line=0.2, outer=TRUE)
		
		# Output
# 		out<-data.frame(tw,Ktauestind,Ktaupind)
# 		colnames(out)<-c("rolling window","Kendall tau estimate","Kendall tau p-value")
	  out<-data.frame(Ktauestind)
    rownames(out)<-tw
		return(out)	
}