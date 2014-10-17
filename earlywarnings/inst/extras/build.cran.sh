/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../
/usr/bin/R CMD check --as-cran earlywarnings_1.1.15.tar.gz
/usr/bin/R CMD INSTALL earlywarnings_1.1.15.tar.gz

