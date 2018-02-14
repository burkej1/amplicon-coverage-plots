# Functions used in amplicon coverage plotting

# Function to create colour scales based on the max value on the given plot
create_plot_colourscale <- function(maxvalue) {
  colourbreaks <- rescale(c(0, 50, 51, 100, 101, 200, 201, 1000, 1001, maxvalue), to = c(0, 1), from = c(0, maxvalue))
  colourscale <- rep(viridis_pal()(5), each = 2)
  # colourscale <- c("black", colourscale)
  colourfunction <- scale_fill_gradientn(colors = colourscale, values = colourbreaks)
  return(colourfunction)
}

# Function to restructure a dataframe for heatmap plotting (note: will throw an error with a single sample)
to_heatmap_matrix <- function(df, replace = FALSE, highvalue = 1000, lowvalue = 0, lowfloor = 0) {
  df_restructure <- dcast(df, amplicon_f~sample, fill=0)
  fix_rownames <- df_restructure[,-1]
  rownames(fix_rownames) <- df_restructure[,1]
  if (replace == TRUE) {
    fix_rownames <- replace(fix_rownames, fix_rownames > highvalue, highvalue)
    fix_rownames <- replace(fix_rownames, fix_rownames < lowvalue, lowfloor)
  }
  return(data.matrix(fix_rownames))
}


# Labels values in a dataframe (for the overlayed dataframes)
label_values <- function(x, value) {
  gsub("^", paste(value, ": ", sep=""), x)
}

# Replaces all values in rows/columns with the output from a function called on those rows/columns
hoverframe_metrics <- function(dm, function_name, by = "row") {
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
  if (by == "amplicon") {
    for (n in 1:length(row.names(dm))) {
      rowmean <- mean(dm[n,])
      normalized_row <- sapply(dm[n,], function(x, m) x / m, rowmean)
      dm[n,] <- normalized_row
    }
  } else if (by == "sample") {
    for (n in 1:length(colnames(dm))) {
      colmean <- mean(dm[,n])
      normalized_col <- sapply(dm[,n], function(x, m) x / m, colmean)
      dm[,n] <- normalized_col
    }
  } else {
    print("Invalid option")
  }
  return(dm)
}

# dm <- plotdata_all  # Testing
# Creates dataframe with extra information for the heatmap plots
generate_hoverframe <- function(dm) {
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
