#' Description: Quantify intermediate stability
#'
#' @param dat Input data matrix (variables x samples)
#' @param meta Metadata (samples x factors). This should contain for each sample 
#'           the following self-explanatory fields: subjectID, time, data. 
#' @param reference.point Calculate stability of the data w.r.t. this point. 
#'                    By default the intermediate range is used (min + (max - min)/2)
#' 
#' @return A vector listing the intermediate stability for each variable
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

intermediate_stability <- function (dat, meta, reference.point = NULL) {

  df <- meta
  stabilities <- c()	      
  for (i in 1:nrow(dat)) {	      
    df$data <- dat[i,]
    stabilities[[i]] <- estimate_stability(df, reference.point)
  }

  if (!is.null(rownames(dat))){
    names(stabilities) <- rownames(dat)
  }

  stabilities

}


#' estimate_stability
#'
#' Description: Quantify intermediate stability with respect to a given reference point. 
#'
#' @param df Input data frame (samples x variables). This should contain for each sample 
#'           the following self-explanatory fields: subjectID, time, data. 
#' @param reference.point Optional. Calculate stability of the data w.r.t. this point. 
#'                        By default the intermediate range is used (min + (max - min)/2)
#' 
#' @return Estimated intermediate stability for the given data.
#'
#' @details The stability is quantified by calculating correlation between 
#'          (i) deviation from the reference point at the first time point
#' 	    (ii) deviation between the first and second time points.
#'          The first and last time point for each subject is used.  
#' 	    Samples with missing data, and subjects with less than two time point are excluded.	   
#'
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples 
#'   #df <- data.frame(list(
#'   #	  subjectID = rep(paste("subject", 1:50, sep = "-"), each = 2), 
#'   #	  time = rep(1:2, 50), 
#'   #	  data = rnorm(100)))
#'   # s <- estimate_stability(df, reference.point = NULL)
#'
#' @keywords internal

estimate_stability <- function (df, reference.point = NULL) {

  # Remove NAs
  df <- df[!is.na(df$data),]

  # Detect intermediate value in the overall data
  if (is.null(reference.point)) {
    reference.point <- mean(range(df$data))
  }
  
  # Remove subjects with only one measurement
  df <- df[df$subjectID %in% names(which(table(df$subjectID) > 1)),]

  # Split data by subject
  spl <- split(df, df$subjectID)    

  # If there are multiple time points, keep the first and last one
  spl <- lapply(spl, function (tab) {tab[c(which.min(tab$time), which.max(tab$time)),]})

  # For each subject, calculate distance from the stability point
  # at the baseline time point
  baseline.distance <- abs(sapply(spl, function (x) {x[which.min(x$time), "data"] - reference.point}))

  # For each subject, calculate deviation between the first and second time point
  followup.distance <- abs(sapply(spl, function (x) {x[which.max(x$time), "data"] - x[which.min(x$time), "data"]}))

  stability <- cor(baseline.distance, followup.distance)  

  stability

}



