# earlywarnings-R

This repository contains the R package toolbox for analysing Early Warning Signals for Critical Transitions that is also submitted to CRAN.

From this repository it is also possible to install the toolbox directly to R by install_git(devtools).

The project is based on:
Dakos V, Carpenter SR, Brock WA, Ellison AM, Guttal V, et al. (2012) Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data. PLoS ONE 7(7): e41010. doi:10.1371/journal.pone.0041010

More can be found on http://www.early-warning-signals.org/

------------------------------------------------------------

## Installing the package

### Install the package from Github in R

Note: if dependencies are missing you may wish to run the installation.script first

```{r}
library(devtools); 
install_github(repo = "earlywarnings-R", username = "earlywarningtoolbox", subdir = "earlywarnings", ref = "master")
```

### Clone the repository & install locally:

Run on command line:
<pre><code>git clone git@github.com:earlywarningtoolbox/earlywarnings-R.git
./build.sh
</pre></code>

