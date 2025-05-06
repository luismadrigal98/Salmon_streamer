mark_outliers_IQR <- function(x, method = 'heuristic', 
                            multiplier = 1.5, na.rm = TRUE) {
  #' Remove outliers using the IQR method (heuristic) or a cluster based on (slower)
  #' 
  #' @param x A numeric vector or data frame column.
  #' @param multiplier A numeric value to adjust the sensitivity of the outlier detection.
  #' @param na.rm Logical. Should missing values be removed?
  #' 
  #' ___________________________________________________________________________ 
  
  if (!is.numeric(x)) {
    warning("Input to is_not_outlier_iqr is not numeric.")
    break
  }
  
  q1 <- stats::quantile(x, 0.25, na.rm = na.rm)
  q3 <- stats::quantile(x, 0.75, na.rm = na.rm)
  iqr_val <- stats::IQR(x, na.rm = na.rm)
  
  lower_bound <- q1 - multiplier * iqr_val
  upper_bound <- q3 + multiplier * iqr_val
  
  # Handle cases where IQR is 0 (e.g., all values are the same in a small group)
  # In such cases, no value should be considered an outlier by this method.
  if (is.na(iqr_val) || iqr_val == 0) {
    return(rep(TRUE, length(x)))
  }
  
  return(x >= lower_bound & x <= upper_bound)
}