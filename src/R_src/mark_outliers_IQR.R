remove_outliers <- function(x, group, method = 'heuristic', 
                            multiplier = 1.5, na.rm = TRUE) {
  #' Remove outliers using the IQR method (heuristic) or a cluster based on (slower)
  #' 
  #' @param x A numeric vector or data frame column.
  #' @param group A grouping variable (e.g., a factor or character vector) to define clusters.
  #' @param method A character string specifying the method for outlier detection.
  #'               Options are: `heuristic` (IQR based) and `k-means`.
  #' @param multiplier A numeric value to adjust the sensitivity of the outlier detection.
  #' @param na.rm Logical. Should missing values be removed?
  #' 
  #' ___________________________________________________________________________ 
  
  if (!is.numeric(x)) {
    warning("Input to is_not_outlier_iqr is not numeric. Returning TRUE for all non-numeric elements.")
    return(rep(TRUE, length(x))) # Or handle as an error, or return FALSE
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

PCA_data_filtered_iqr <- PCA_data |>
  dplyr::group_by(Group_merged) |>
  dplyr::filter(
    dplyr::if_all(
      .cols = dplyr::starts_with("PC"), # Selects all PCs
      .fns = ~ is_not_outlier_iqr(.)    # Apply the IQR non-outlier check
    )
  ) |>
  dplyr::ungroup()

return(PCA_data_filtered_iqr)