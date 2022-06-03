

find.minima <- function (f) {
  find.maxima(-f)
}

find.maxima <- function (f) {

    f2 <- c(Inf, -f, Inf)
    cnt <- 1
    ops <- c()
    opcnt <- 0
    while (cnt < length(f2)) {
      if (f2[[cnt + 1]] - f2[[cnt]] <= 0) {
	while (f2[[cnt + 1]] - f2[[cnt]] <= 0) {
	  cnt <- cnt + 1
        }
	ind1 <- cnt - 1
	while (f2[[cnt + 1]] - f2[[cnt]] == 0) {
	  cnt <- cnt + 1
        }
	if (f2[[cnt + 1]] - f2[[cnt]] > 0) {
	  ind2 <- cnt - 1
    	  opcnt <- opcnt + 1
	  ops[[opcnt]] <- round(mean(c(ind1, ind2)))
	} else if (f2[[cnt + 1]] - f2[[cnt]] < 0) {
	  ind2 <- NULL
	}
      }
      cnt <- cnt + 1
    }
    ops 
}

