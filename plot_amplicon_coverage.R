library(ggplot2)
library(readr)
library(gplots)
library(RColorBrewer) 
library(reshape2) 
library(heatmaply)


# Function to restructure a dataframe for heatmap plotting (note: will throw an error with a single sample)
to_heatmap_matrix <- function(df, replace = FALSE) {
  df_restructure <- dcast(df, amplicon_f~sample, fill=0)
  fix_rownames <- df_restructure[,-1]
  rownames(fix_rownames) <- df_restructure[,1]
  if (replace == TRUE) {
    fix_rownames <- replace(fix_rownames, fix_rownames > 1000, 1000)
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

# Input file (should write to take a command line input)
amplicon_metrics <- read_delim("/Users/Jared/Documents/Code/amplicon-coverage-plots/run1-7_amplicon_metrics_all.tsv", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)

headers <- c("chrom", "start", "end", "amplicon_f", "amplicon_r", "reads", "bases_c", "len", "pct", "sample")
colnames(amplicon_metrics) <- headers
plotdata <- amplicon_metrics[,c("amplicon_f", "sample", "reads")]  # Subset dataframe
# plotdata <- plotdata[1:2386,]  # Subsetting for testing purposes

# Creating matrices
plotdata_all_capped <- to_heatmap_matrix(plotdata, TRUE)  # All samples and all amplicons, capped for visualisation
plotdata_all <- to_heatmap_matrix(plotdata, FALSE)  # Uncapped amplicon data, use this for overlay matrices

# Generating overlay matrix
full_hovertext <- generate_hoverframe(plotdata_all)

allplot <- heatmaply(plotdata_all_capped, 
                     fontsize_col = 8, 
                     fontsize_row = 8, 
                     file = "all_coverage.html", 
                     custom_hovertext = full_hovertext)


# Extract the gene names from the amplicons and create a heatmap for each gene
gene_names <- unique(gsub("_.+_.+", "", plotdata$amplicon_f))
for (gene in gene_names) {
  main_frame <- plotdata[grep(gene, plotdata$amplicon_f),]
  gene_subset_capped <- to_heatmap_matrix(main_frame, TRUE)
  gene_subset <- to_heatmap_matrix(main_frame, FALSE)
  hovertext <- generate_hoverframe(gene_subset)
  plotname <- paste(gene, "_coverage.html", sep="")
  geneplot <- heatmaply(gene_subset_capped, 
                        fontsize_col = 8, 
                        fontsize_row = 8, 
                        custom_hovertext = hovertext, 
                        file = plotname)
}




# # Testing
# plotdata_brca1 <- to_heatmap_matrix(plotdata[grep("BRCA1", plotdata$amplicon_f),], TRUE)  # Subsetting to amplicons in a single gene
# 
# plot <- heatmaply(plotdata_all, 
#           fontsize_col = 8, fontsize_row = 8, 
#           margins = c(180, 100), 
#           file = "BRCA1_coverage.html"  # Save as HTML
#           )




