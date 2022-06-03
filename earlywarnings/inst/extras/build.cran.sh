#!/bin/sh
#~/bin/R-4.2.0/bin/R CMD BATCH document.R
~/bin/R-4.2.0/bin/R CMD build ../../
~/bin/R-4.2.0/bin/R CMD check --as-cran earlywarnings_1.1.24.tar.gz
~/bin/R-4.2.0/bin/R CMD INSTALL earlywarnings_1.1.24.tar.gz


