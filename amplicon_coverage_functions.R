# Functions used in amplicon coverage plotting
breaks = c(50, 100, 200, 1000)
maxvalue = 10000
create_plot_colourscale <- function(maxvalue, breaks = c(50, 100, 200, 1000), colours = viridis_pal()(5)) {
  # Creates colour scales based on the given max value, breaks and colour palette.
  # Currently 0 is assumed to be the start of the scale. Number of colours supplied should be equal to n of breaks + 1.
  breaks_formatted <- unlist(lapply(breaks, function(x) c(x, x + 0.0001)))
  colourbreaks <- rescale(c(0, breaks_formatted, maxvalue), to = c(0, 1), from = c(0, maxvalue))
  colourscale <- rep(colours, each = 2)
  colourfunction <- scale_fill_gradientn(colors = colourscale, values = colourbreaks)
  return(colourfunction)
}

to_heatmap_matrix <- function(df, replace = FALSE, highvalue = 1000, lowvalue = 0, lowfloor = 0) {
  # Function to restructure a dataframe for heatmap plotting (note: will throw an error with a single sample)
  df_restructure <- dcast(df, amplicon_f~sample, fill=0)
  fix_rownames <- df_restructure[,-1]
  rownames(fix_rownames) <- df_restructure[,1]
  if (replace == TRUE) {
    fix_rownames <- replace(fix_rownames, fix_rownames > highvalue, highvalue)
    fix_rownames <- replace(fix_rownames, fix_rownames < lowvalue, lowfloor)
  }
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
