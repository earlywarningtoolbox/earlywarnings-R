---
title: "earlywarnings vignette"
author: "Vasilis Dakos and Leo Lahti"
date: "2015-07-06"
output:
  html_document:
    toc: true
    number_sections: true
    theme: united
    highlight: pygments
---

<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial}
  %\usepackage[utf8]{inputenc}
-->



earlywarnings R package
===========

## Installation

### Installing and loading the release version

Note: if dependencies are missing you may wish to run the [../installation.script](installation.script) first


```r
install.packages("earlywarnings")
```

### Installing and loading the experimental development version


```r
library(devtools)
install_github("earlywarningtoolbox/earlywarnings-R/earlywarnings")
```

### Clone the repository & install locally:

Run on command line:
<pre><code>git clone git@github.com:earlywarningtoolbox/earlywarnings-R.git
./build.sh
</pre></code>

### Loading the package


```r
library(earlywarnings)  
```


## Potential analysis

Potential analysis, used for instance in [Hirota et al. Science, 334, 232-235.](http://www.sciencemag.org/content/334/6053/232.long)


```r
# Create simulated example data
X <- c(rnorm(1000, mean = 0), rnorm(1000, mean = -2), 
 	           rnorm(1000, mean = 2))
param <- seq(0,5,length=3000) 

# Run potential analysis
res <- movpotential_ews(X, param)

# Visualize
p <- PlotPotential(res$res, title = '', 
	       	   xlab.text = '', ylab.text = '', 
		   cutoff = 0.5, plot.contours = TRUE, binwidth = 0.2)
print(p)
```

![plot of chunk movpotential](figure/movpotential-1.png) 


### Licensing and Citations

This work can be freely used, modified and distributed under the 
[Two-clause (Free)BSD license](http://en.wikipedia.org/wiki/BSD\_licenses).

Kindly cite the work as 'Vasilis Dakos and Leo Lahti (2014). earlywarnings R package. URL: https://github.com/earlywarningtoolbox/earlywarnings-R/tree/master/earlywarnings'.


```r
citation("earlywarnings")
```

```
## Error in tools:::.parse_CITATION_file(file, meta$Encoding): non-ASCII input in a CITATION file without a declared encoding
```

### Session info

This vignette was created with


```r
sessionInfo()
```

```
## R version 3.2.1 (2015-06-18)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## Running under: Ubuntu 15.04
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] earlywarnings_1.1.21 tseries_0.10-34      tgp_2.4-11          
## [4] moments_0.14         ggplot2_1.0.1        knitr_1.10.5        
## [7] scimapClient_0.2.1  
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.11.6        cluster_2.0.2      magrittr_1.5      
##  [4] maps_2.3-9         MASS_7.3-41        nortest_1.0-3     
##  [7] Kendall_2.2        munsell_0.4.2      som_0.3-5         
## [10] colorspace_1.2-6   lattice_0.20-31    quadprog_1.5-5    
## [13] stringr_1.0.0      plyr_1.8.3         fields_8.2-1      
## [16] tools_3.2.1        grid_3.2.1         spam_1.0-1        
## [19] gtable_0.1.2       KernSmooth_2.23-14 lmtest_0.9-34     
## [22] digest_0.6.8       RJSONIO_1.3-0      reshape2_1.4.1    
## [25] formatR_1.2        rpart_4.1-9        evaluate_0.7      
## [28] maptree_1.4-7      labeling_0.3       stringi_0.5-5     
## [31] scales_0.2.5       boot_1.3-16        proto_0.3-10      
## [34] zoo_1.7-12
```




