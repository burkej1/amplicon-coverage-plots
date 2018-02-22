# Functions used in amplicon coverage plotting

create_plots <- function(df, basename, singlegene = FALSE, colourmode = "discrete_na", plot_type = "depth") {
  # Creates plots from the given dataframe using a bunch of defaults
  # Creating side annotation dataframes (plate number extraction specific to BRASTRAP)
  plate_vector <- data.frame(factor(gsub("^.+?_(.+?_[0-9]+)-.+", "\\1", colnames(df))))  # Coerced vector to df to add label
  colnames(plate_vector) <- "Plate Number"
  
  if (singlegene == TRUE) {
    # If a single gene is being plotted colour the rows by exon
    row_vector <- data.frame(factor(gsub("^.+?_(.+)_.+", "\\1", row.names(df))))
    colnames(row_vector) <- "Exon"
  } else {
    # Otherwise colour the rows by gene
    row_vector <- data.frame(factor(gsub("_.+_.+", "", row.names(df))))
    colnames(row_vector) <- "Gene"
  }
  
  # Creating hovertext
  plot_hovertext <- generate_hoverframe(df)
  
  if (plot_type == "depth" || plot_type == "all") {
    # Creating colour scale
    if (colourmode == "discrete_na") {
      # Pseudo-discrete scale with white for 0 coverage
      depth_plot_colourscale <- create_plot_colourscale(max(df), c(1, 50, 100, 200, 500, 1000), c("white", viridis_pal()(6)))
    } else if (colourmode == "discrete") {
      # Pseudo-discrete scale
      depth_plot_colourscale <- create_plot_colourscale(max(df), c(50, 100, 200, 500, 1000), viridis_pal()(6))
    } else if (colourmode == "continuous") {
      # Continuous colourscale
      depth_plot_colourscale <- create_plot_colourscale(max(df), c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), 
                                                                          viridis_pal()(11), type = "continuous")
    }
    # Creating the depth heatmap
    depth_heatmap <- heatmaply(df, 
                               main = paste(basename, "- Depth"), 
                               fontsize_col = 8,
                               fontsize_row = 8,
                               hide_colorbar = TRUE,
                               scale_fill_gradient_fun = depth_plot_colourscale,
                               col_side_colors = plate_vector,
                               row_side_colors = row_vector,
                               file = paste("coverage_", basename, ".html", sep=""),
                               margins = c(200, 150),
                               custom_hovertext = plot_hovertext)
  }
  if (plot_type == "amplicon_norm" || plot_type == "all") {
    ampnorm_plotname <- paste(basename, "_amplicon-normalised_coverage.html", sep="")
    ampnorm_matrix <- normalize_matrix(df, by = "amplicon")
    ampnorm_colourscale <- create_plot_colourscale(max(ampnorm_matrix), breaks = c(0.5, 1.0, 1.5, 2.0), colours = viridis_pal()(5))
    geneplot_ampnorm <- heatmaply(ampnorm_matrix,
                                  main = paste(basename, "- Amplicon Normalised"),
                                  hide_colorbar = TRUE,
                                  fontsize_col = 8,
                                  fontsize_row = 8,
                                  scale_fill_gradient_fun = ampnorm_colourscale,
                                  col_side_colors = plate_vector,
                                  row_side_colors = row_vector,
                                  margins = c(200, 150),
                                  file = ampnorm_plotname,
                                  custom_hovertext = plot_hovertext)
  }
  if (plot_type == "sample_norm" || plot_type == "all") {
    sampnorm_plotname <- paste(basename, "_sample-normalised_coverage.html", sep="")
    sampnorm_matrix <- normalize_matrix(df, by = "sample")
    sampnorm_colourscale <- create_plot_colourscale(max(sampnorm_matrix), breaks = c(0.5, 1.0, 1.5, 2.0), colours = viridis_pal()(5))
    geneplot_sampnorm <- heatmaply(sampnorm_matrix,
                                   main = paste(basename, "- Sample Normalised"),
                                   hide_colorbar = TRUE,
                                   fontsize_col = 8,
                                   fontsize_row = 8,
                                   scale_fill_gradient_fun = sampnorm_colourscale,
                                   col_side_colors = plate_vector,
                                   row_side_colors = row_vector,
                                   margins = c(200, 150),
                                   file = sampnorm_plotname,
                                   custom_hovertext = plot_hovertext)
  }
}

create_plot_colourscale <- function(maxvalue, breaks = c(50, 100, 200, 1000), colours = viridis_pal()(5), type = "discrete") {
  # Creates colour scales based on the given max value, breaks and colour palette.
  if (type == "discrete") {
    # Creates a pseudo discrete continuous scale from the given breaks
    # Currently 0 is assumed to be the start of the scale. Number of colours supplied should be equal to n of breaks + 1.
    breaks_formatted <- unlist(lapply(breaks, function(x) c(x, x + 0.0001)))
    colourbreaks <- rescale(c(0, breaks_formatted, maxvalue), to = c(0, 1), from = c(0, maxvalue))
    colourscale <- rep(colours, each = 2)
    colourfunction <- scale_fill_gradientn(colors = colourscale, values = colourbreaks)
  } else if (type == "continuous") {
    # Creates a continous scale capped at the final break
    colourbreaks <- rescale(c(0, breaks, maxvalue), to = c(0, 1), from = c(0, maxvalue))
    colourscale <- c(colours, tail(colours, n=1))
    colourfunction <- scale_fill_gradientn(colors = colourscale, values = colourbreaks)
  }
  return(colourfunction)
}

to_heatmap_matrix <- function(df) {
  # Function to restructure a dataframe for heatmap plotting
  df_restructure <- dcast(df, amplicon_f~sample, fill=0)
  fix_rownames <- df_restructure[,-1, drop=FALSE]  # Drop option allows it to work with single samples
  rownames(fix_rownames) <- df_restructure[,1]
  return(data.matrix(fix_rownames))
}


label_values <- function(x, value) {
  # Labels values in a dataframe (for the overlayed dataframes)
  gsub("^", paste(value, ": ", sep=""), x)
}

hoverframe_metrics <- function(dm, function_name, by = "row") {
  # Replaces all values in rows/columns with the output from a function called on those rows/columns
  FUN <- match.fun(function_name)
  if (by == "column") {
    # Iterate through columns
    for (n in 1:length(colnames(dm))) {
      dm[,n] <- FUN(dm[,n])  # Replacing each column with the function output
    }
  } else if (by == "row") {
    for (n in 1:length(row.names(dm))) {
      dm[n,] <- FUN(dm[n,])
    }
  } else {
    print("Unknown argument passed for 'by = '")
  }
  return(dm)
}

normalize_matrix <- function(dm, by = "amplicon") {
  # Normalizes the value of each cell in a matrix either by dividing by the corresponding sample average (by sample) or
  # amplicon average (by amplicon)
  if (by == "amplicon") {
    for (n in 1:length(row.names(dm))) {
      rowmean <- mean(dm[n,])
      if (rowmean == 0) {
        normalized_row <- sapply(dm[n,], function(x) x * 0)
      } else {
        normalized_row <- sapply(dm[n,], function(x, m) x / m, rowmean)
      }
      dm[n,] <- normalized_row
    }
  } else if (by == "sample") {
    for (n in 1:length(colnames(dm))) {
      colmean <- mean(dm[,n])
      if (colmean == 0) {
        normalized_col <- sapply(dm[,n], function(x) x * 0, colmean)
      } else {
      normalized_col <- sapply(dm[,n], function(x, m) x / m, colmean)
      }
      dm[,n] <- normalized_col
    }
  } else {
    print("Invalid option")
  }
  return(dm)
}

generate_hoverframe <- function(dm) {
  # Creates dataframe with extra information for the heatmap plots
  rowmeans <- sapply(data.frame(hoverframe_metrics(dm, "mean", by = "row")), label_values, "Amplicon Mean")
  rowsd <- sapply(data.frame(hoverframe_metrics(dm, "sd", by = "row")), label_values, "Amplicon SD")
  colmeans <- sapply(data.frame(hoverframe_metrics(dm, "mean", by = "column")), label_values, "Sample Mean")
  colsd <- sapply(data.frame(hoverframe_metrics(dm, "sd", by = "column")), label_values, "Sample SD")
  original_depths <- sapply(data.frame(dm), label_values, "Amplicon Depth")
  
  if (exists("new_df")) {
    rm("new_df")
  }
  for (n in 1:length(colnames(dm))) {
    newcolumn <- data.frame(paste(original_depths[, n], colmeans[, n], colsd[, n], rowmeans[, n], rowsd[, n], sep="\n"))
    if (exists("new_df")) {
      new_df <- cbind(new_df, newcolumn)
    } else {
      new_df <- data.frame(newcolumn)
    }
  }
  return(new_df)
}
