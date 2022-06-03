# R Installation script for earlywarnings package.
# Install dependencies, then the package.

# Needed for earlywarnings
install.packages("moments")
install.packages("quadprog")
install.packages("nortest")
install.packages("Kendall")
install.packages("som")
install.packages("akima")
install.packages("tseries")
install.packages("splus2R")
install.packages("ifultools")
install.packages("sapa")
install.packages("wmtsa")
install.packages("scatterplot3d")
#system("wget http://cran.r-project.org/src/contrib/Archive/fractal/fractal_1.1-1.tar.gz")
#install.packages("fractal_1.1-1.tar.gz", repos = NULL)

library(devtools)
install_github(repo = "earlywarnings-R", username = "earlywarningtoolbox", subdir = "earlywarnings", ref ="master")
#install_github(repo = "earlywarnings-R", username = "earlywarningtoolbox", type = "source", ref = "develop")

