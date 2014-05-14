/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../
/usr/bin/R CMD check --as-cran earlywarnings_1.0.64.tar.gz
/usr/bin/R CMD INSTALL earlywarnings_1.0.64.tar.gz

