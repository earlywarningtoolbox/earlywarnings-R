R CMD BATCH document.R
R CMD check earlywarnings
R CMD build earlywarnings
R CMD INSTALL earlywarnings_1.0.34.tar.gz
