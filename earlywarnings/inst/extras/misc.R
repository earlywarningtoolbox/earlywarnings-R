## Multimodality detection

Estimate multimodality score and number of modes based on bootstrapped potential analysis.

```{r bimodality, message=FALSE, warning=FALSE, fig.width=5, fig.height=5}
# Potential analysis with bootstrap
library(earlywarnings)

# Example Data
X <- as.matrix(rbind(c(rnorm(100, mean = 0), rnorm(100, mean = 5)), 
           c(rnorm(200, mean = 0))))
m <- multimodality_score(X, detection.threshold = 1, bs.iterations = 20, detection.limit = 3)

# Plot the original data for feature i 
# together with the estimated density maxima and minima
i <- 1
plot(density(X[i,])); 
abline(v = m$results[[i]]$maxima)
abline(v = m$results[[i]]$minima, lty = 2)
```

