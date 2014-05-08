#' intermediate_stability
#'
#' Quantify stability with respect to a given reference point for all variables.
#'
#' @param dat Input data matrix (variables x samples)
#' @param meta Metadata (samples x factors). This should contain for each sample 
#'           the following self-explanatory fields: subjectID and time. 
#' @param reference.point Calculate stability of the data w.r.t. this point. 
#'                    By default the intermediate range is used (min + (max - min)/2)
#' @param method "lm" (linear model) or "correlation"; the linear model takes time into account 
#' 	         as a covariate 
#'
#' @return A list with following elements: 
#' 	     stability: vector listing the intermediate stability for each variable
#'	     data: processed data sets used for calculations	    
#'
#'
#' @details This method decomposes the data set into differences between
#' consecutive time points. For each variable and time point we calculate
#' for the data values: (i) the distance from reference point; (ii)
#' distance from the data value at the consecutive time point. The
#' "correlation" method calculates correlation between these two
#' variables. Negative correlations indicate that values closer to
#' reference point tend to have larger shifts in the consecutive time
#' point. The "lm" method takes the time lag between the consecutive time
#' points into account as this may affect the comparison and is not taken
#' into account by the straightforward correlation. Here the coefficients
#' of the following linear model are used to assess stability:
#' abs(change) ~ time + abs(start.reference.distance)
#'
#' @export
#'
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples 
#'   # Create simulated example data
#'   dat <- matrix(rnorm(1000), nrow = 10)
#'   rownames(dat) <- paste("Variable", 1:nrow(dat), sep = "")
#'   colnames(dat) <- paste("Sample", 1:ncol(dat), sep = "")
#'   meta <- data.frame(list(
#'     	  sampleID = colnames(dat), 
#'   	  subjectID = rep(paste("subject", 1:50, sep = "-"), each = 2), 
#'	  time = rep(1:2, 50)))
#'   # Intentionally make point very unstable around 0
#'   dat[1, meta$time == 2] <- 1/abs(dat[1, meta$time == 1])
#'   s <- intermediate_stability(dat, meta, reference.point = 0)
#'
#' @keywords early-warning

intermediate_stability <- function (dat, meta, reference.point = NULL, method = "lm") {

  df <- meta
  stabilities <- c()	      
  stabilities.right <- c()	      
  stabilities.left <- c()	      
  data <- list()
  for (i in 1:nrow(dat)) {	      
    df$data <- dat[i,]
    stab <- estimate_stability(df = df, reference.point = reference.point, method = method)
    stabilities[[i]] <- stab$stability
    stabilities.right[[i]] <- stab$stability.right
    stabilities.left[[i]] <- stab$stability.left
    data[[i]] <- stab$data
  }

  if (!is.null(rownames(dat))){
    names(stabilities) <- rownames(dat)
    names(stabilities.right) <- rownames(dat)
    names(stabilities.left) <- rownames(dat)
    names(data) <- rownames(dat)
  }

  list(stability = stabilities, stability.right = stabilities.right, stability.left = stabilities.left, data = data)

}



#' estimate_stability
#'
#' Description: Quantify intermediate stability with respect to a given reference point. 
#'
#' @param df Input data frame (samples x variables). This should contain for each sample 
#'           the following self-explanatory fields: subjectID, time, data. 
#' @param reference.point Optional. Calculate stability of the data w.r.t. this point. 
#'                        By default the intermediate range is used (min + (max - min)/2)
#' @param method "lm" (linear model) or "correlation"; the linear model takes time into account 
#' 	         as a covariate 
#' 
#' @return A list with following elements: 
#' 	     stability: estimated stability
#'	     data: processed data set used in calculations	    
#'
#' @details This method decomposes the data set into differences between
#' consecutive time points. For each variable and time point we calculate
#' for the data values: (i) the distance from reference point; (ii)
#' distance from the data value at the consecutive time point. The
#' "correlation" method calculates correlation between these two
#' variables. Negative correlations indicate that values closer to
#' reference point tend to have larger shifts in the consecutive time
#' point. The "lm" method takes the time lag between the consecutive time
#' points into account as this may affect the comparison and is not taken
#' into account by the straightforward correlation. Here the coefficients
#' of the following linear model are used to assess stability:
#' abs(change) ~ time + abs(start.reference.distance). 
#' Samples with missing data, and subjects with less than two time point are excluded.	   
#'
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @export
#' @examples 
#'   #df <- data.frame(list(
#'   #	  subjectID = rep(paste("subject", 1:50, sep = "-"), each = 2), 
#'   #	  time = rep(1:2, 50), 
#'   #	  data = rnorm(100)))
#'   # s <- estimate_stability(df, reference.point = NULL, method = "lm")
#'
#' @keywords internal


estimate_stability <- function (df, reference.point = NULL, method = "lm") {

  # Remove NAs
  df <- df[!is.na(df$data),]

  # Detect intermediate value in the overall data if reference point not given
  if (is.null(reference.point)) {
    reference.point <- mean(range(df$data))
  }
  
  # Remove subjects with only one measurement
  df <- df[df$subjectID %in% names(which(table(df$subjectID) > 1)),]

  # Split data by subject
  spl <- split(df, df$subjectID)    

  dfis <- NULL
  for (spli in spl) {

    # Ensure the measurements are ordered in time
    spli <- spli[order(spli$time),]

    # Calculate differences between consecutive time points; 
    # and the initial values; and their distance from reference
    data.difs <- diff(spli$data)
    time.difs <- diff(spli$time)
    start.points <- spli$data[-nrow(spli)]
    start.reference.distance <- start.points - reference.point

    # Organize into data frame
    dfi <- data.frame(change = data.difs, time = time.difs, start = start.points, start.reference.distance = start.reference.distance)

    # Add to the collection
    dfis <- rbind(dfis, dfi)

  }

  # Simplified stability calculation (do not consider time effect)
  stability <- NULL
  if (method == "correlation") {
    # For each subject, check distance from the stability point
    # at the baseline time point
    baseline.distance <- abs(dfis$start.reference.distance)
    # For each subject, calculate deviation between the first and second time point
    followup.distance <- abs(dfis$change)
    stability <- cor(baseline.distance, followup.distance)  
  } else if (method == "lm") {
    # Advanced calculation, take time into account with linear model (also possible to check p-values later if needed)
    stability <- coef(summary(lm(abs(change) ~ time + abs(start.reference.distance), data = dfis)))["abs(start.reference.distance)", "Estimate"]
  }

  stability.left <- stability.right <- NA

  dfis.left <- subset(dfis, start.reference.distance < 0)
  if (nrow(dfis.left)>10) {
    # Negative values for low stability
    stability.left <- -coef(summary(lm(change ~ time + abs(start.reference.distance), data = dfis.left)))["abs(start.reference.distance)", "Estimate"]
  }

  dfis.right <- subset(dfis, start.reference.distance > 0)
  if (nrow(dfis.right)>10) {
    # Negative values for low stability
    stability.right <- -coef(summary(lm(change ~ time + abs(start.reference.distance), data = dfis.right)))["abs(start.reference.distance)", "Estimate"]
  }

  list(stability = stability, stability.right = stability.right, stability.left = stability.left, data = dfis)

}
