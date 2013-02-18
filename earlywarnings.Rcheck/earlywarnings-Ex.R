pkgname <- "earlywarnings"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('earlywarnings')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ch_ews")
### * ch_ews

flush(stderr()); flush(stdout())

### Name: ch_ews
### Title: Description: Conditional Heteroskedasticity
### Aliases: ch_ews
### Keywords: early-warning

### ** Examples

data(foldbif)
out=ch_ews(foldbif, winsize=50, alpha=0.05, optim=TRUE, lags)



cleanEx()
nameEx("ddjnonparam_ews")
### * ddjnonparam_ews

flush(stderr()); flush(stdout())

### Name: ddjnonparam_ews
### Title: Description: Drift Diffusion Jump Nonparametrics Early Warning
###   Signals
### Aliases: ddjnonparam_ews
### Keywords: early-warning

### ** Examples

data(foldbif)
output<-ddjnonparam_ews(foldbif,bandwidth=0.6,na=500,
logtransform=TRUE,interpolate=FALSE)



cleanEx()
nameEx("generic_ews")
### * generic_ews

flush(stderr()); flush(stdout())

### Name: generic_ews
### Title: Description: Generic Early Warning Signals
### Aliases: generic_ews
### Keywords: early-warning

### ** Examples

data(foldbif)
 out=generic_ews(foldbif,winsize=50,detrending="gaussian",
 bandwidth=5,logtransform=FALSE,interpolate=FALSE)



cleanEx()
nameEx("livpotential_ews")
### * livpotential_ews

flush(stderr()); flush(stdout())

### Name: livpotential_ews
### Title: Description: Potential Analysis
### Aliases: livpotential_ews
### Keywords: early-warning

### ** Examples

data(foldbif)
res <- livpotential_ews(foldbif)
plot(res$xi, res$pot)



cleanEx()
nameEx("movpotential_ews")
### * movpotential_ews

flush(stderr()); flush(stdout())

### Name: movpotential_ews
### Title: Description: Moving Average Potential
### Aliases: movpotential_ews
### Keywords: early-warning

### ** Examples

X = c(rnorm(1000, mean = 0), rnorm(1000, mean = -2), rnorm(1000, mean = 2))
 param = seq(0,5,length=3000)
 res <- movpotential_ews(X, param, npoints = 100, thres = 0.003)



cleanEx()
nameEx("sensitivity_ews")
### * sensitivity_ews

flush(stderr()); flush(stdout())

### Name: sensitivity_ews
### Title: Description: Sensitivity Early Warning Signals
### Aliases: sensitivity_ews
### Keywords: early-warning

### ** Examples

data(foldbif)
output=sensitivity_ews(foldbif,indicator="sd",detrending="gaussian",
incrwinsize=25,incrbandwidth=20)



cleanEx()
nameEx("surrogates_ews")
### * surrogates_ews

flush(stderr()); flush(stdout())

### Name: surrogates_ews
### Title: Description: Surrogates Early Warning Signals
### Aliases: surrogates_ews
### Keywords: early-warning

### ** Examples

data(foldbif);
output=surrogates_ews(foldbif,indicator="sd",winsize=50,detrending="gaussian",
bandwidth=10,boots=200,logtransform=FALSE,interpolate=FALSE)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
