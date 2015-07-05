/usr/local/bin/R CMD BATCH document.R
/usr/local/bin/R CMD build ../../
/usr/local/bin/R CMD check --as-cran earlywarnings_1.1.21.tar.gz
/usr/local/bin/R CMD INSTALL earlywarnings_1.1.21.tar.gz


