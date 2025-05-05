process_PCA_results <- function(PCA_res, PC_to_retain)
{
  #' This function will process the PCA results obtained trhough prcomp fucntion.
  #' it will return the PC based dataframe with additional information related
  #' to the amount of variance explained by the reatined components.
  #' 
  #' @param PCA_res A list with the PCA results obtained through prcomp function.
  #' 
  #' @param PC_to_retain A numeric value with the number of PC to retain.
  #' 
  #' @return A data frame
  #' ___________________________________________________________________________
  
  # Check if PCA_res is a prcomp object
  assertthat::see_if(class(PCA_res) == 'prcomp')
  
  PCA_data <- PCA_res$x[, 1:PC_to_retain]
  
  # Retrieve the variance explained
  
  PCA_summary <- summary(PCA_res)$importance['Cumulative Proportion', 
                                                 PC_to_retain] 
  
  message(paste0("With ", PC_to_retain, "PCs, you capture ", PCA_summary * 100, 
                 "% of the variance."))
  
  return(PCA_data)
}