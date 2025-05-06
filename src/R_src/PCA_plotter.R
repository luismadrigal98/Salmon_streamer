PCA_plotter <- function(PCA_data, x = 'PC1', y = 'PC2', color = NULL,
                        group = NULL, filename = "PCA_out.pdf")
{
  #' Plot PCA results
  #' 
  #' @param PCA_data A data frame containing PCA results.
  #' @param x The x-axis variable (default is 'PC1').
  #' @param y The y-axis variable (default is 'PC2').
  #' @param color The variable to color the points by (default is NULL).
  #' @param group The variable to group the points by (default is NULL).
  #' @param filename The name of the output file (default is "PCA_out.pdf").
  #' 
  #' @return None
  #' ___________________________________________________________________________
  
  p <- ggplot(data = PCA_data, mapping = aes(x = .data[[x]], y = .data[[y]]))
  
  # Add points with conditional aesthetics
  if (!is.null(color) && !is.null(group)) {
    p <- p + geom_point(aes(color = .data[[color]], shape = .data[[group]]), size = 3)
  } else if (!is.null(color)) { # Only color is specified
    p <- p + geom_point(aes(color = .data[[color]]), size = 3)
  } else if (!is.null(group)) { # Only group (for shape) is specified
    p <- p + geom_point(aes(shape = .data[[group]]), size = 3)
  } else { # Neither color nor group is specified
    p <- p + geom_point(size = 3)
  }
  
  p <- p + theme_bw() +
    labs(title = "PCA Plot", x = x, y = y) +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_text(color = 'black'))
  
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300,
         device = cairo_pdf, limitsize = FALSE)
}