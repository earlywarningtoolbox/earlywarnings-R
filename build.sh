R CMD BATCH document.R
R CMD check earlywarnings --as-cran
R CMD build earlywarnings
R CMD INSTALL earlywarnings_1.0.53.tar.gz
