assign_groups_from_labels <- function(labels, rules, default_group = "Unknown") 
{
  #' Assign Groups Based on Label Patterns
  #'
  #' Assigns a group name to each label based on a list of patterns.
  #' The group name assigned corresponds to the first pattern in the rules list
  #' that the label matches.
  #'
  #' @param labels A character vector of sample labels.
  #' @param rules A named list where names are the desired group names
  #'              and values are the regular expression patterns to identify 
  #'              those groups.
  #'              The order of rules matters, as the first match determines 
  #'              the group.
  #' @param default_group A character string for labels that don't match any rule.
  #'                      Defaults to "Unknown".
  #'
  #' @return A character vector of the same length as labels, containing the
  #'         assigned group name for each label.
  #'
  #' @examples
  #' sample_labels <- c("SFcross_s6_J.P.032", "SWBcross_s2_D.069",
  #'                    "SWBcross_s1_767.P9", "SFcross_s5_J.102", 
  #'                    "Some_Other_Label")
  #' grouping_rules <- list(
  #'   "SF_Parental" = "SFcross.*_.*P",
  #'   "SWB_Parental" = "SWBcross.*_.*P",
  #'   "SF_Cross" = "^SFcross", # Use ^ to anchor at the beginning
  #'   "SWB_Cross" = "^SWBcross" # Use ^ to anchor at the beginning
  #' )
  #' assigned_groups <- assign_groups_from_labels(sample_labels, grouping_rules)
  #' print(assigned_groups)
  #' # Output: [1] "SF_Parental" "SWB_Cross" "SWB_Parental" "SF_Cross" "Unknown"
  #' ___________________________________________________________________________
  
  # Initialize the result vector with the default group
  assigned_groups <- rep(default_group, length(labels))
  
  # Iterate through each label
  for (i in seq_along(labels)) {
    label <- labels[i]
    
    # Iterate through the rules in the specified order
    for (group_name in names(rules)) {
      pattern <- rules[[group_name]]
      
      # Check if the label matches the pattern
      if (grepl(pattern, label, perl = TRUE)) { # Using perl=TRUE for potentially more complex regex
        assigned_groups[i] <- group_name
        break # Stop checking rules for this label once a match is found
      }
    }
  }
  
  return(assigned_groups)
}