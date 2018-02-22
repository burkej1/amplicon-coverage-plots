library(ggplot2)  # General dependency for heatmaply/shinyheatmaply
library(readr)  # For importing tsvs
library(gplots)  # Used for something
library(RColorBrewer)  # Colour palettes
library(reshape2)  # For manipulating dataframes/matrices I think
library(heatmaply)  # Interactive heatmaps
library(scales)  # Rescaling scales
library(colorspace) # For generating colourspaces for rowside annotation

source("/Users/Jared/Documents/Code/amplicon-coverage-plots/amplicon_coverage_functions.R")  # This may not work depending on the R working directory

# # Manual input file paths
# inputpath <- "/Users/Jared/Documents/Project_Results_and_Data/BSTP/run1-7_amplicon_coverage_plotting/run1-7_amplicon_metrics_clean_190218.txt"
# inputpath <- "/Users/Jared/Documents/Jobs/CONFIRM_PROSTRAP_test_plotting/amplicon_metrics_files/CNFRM-column_S3_L001.amplicon-metrics.txt"
inputpath <- "/Users/Jared/Documents/Jobs/CONFIRM_PROSTRAP_test_plotting/amplicon_metrics_files/PSTP-column_S2_L001.amplicon-metrics.txt"

# Reading input
amplicon_metrics <- read_delim(inputpath, 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)
# Adding headers
headers <- c("chrom", "start", "end", "amplicon_f", "amplicon_r", "reads", "bases_c", "len", "pct", "sample")
colnames(amplicon_metrics) <- headers
# Subset dataframe to necessary columns
plotdata <- amplicon_metrics[,c("amplicon_f", "sample", "reads")]  

# Creating main matrix
plotdata_all <- to_heatmap_matrix(plotdata)


# # Plotting single samples
# Clustering function throws a fit with a single column dataframe (understandably) expanding the df to avoid this
plotdata_all_expanded <- cbind(plotdata_all, plotdata_all)
colnames(plotdata_all_expanded)[2] <- c("X")
# Plotting (note, amplicon_norm will not work with a single sample duplicated)
create_plots(plotdata_all_expanded, "PROSTRAP", plot_type = "depth")
gene_names <- unique(gsub("_.+_.+", "", plotdata$amplicon_f))
for (gene in gene_names) {
  print(gene)
  gene_subset <- data.matrix(plotdata_all_expanded[grep(gene, row.names(plotdata_all)),])
  create_plots(gene_subset, gene, singlegene = TRUE, plot_type = "depth")
}


# # # # # # # # Gene plots # # # # # # # # 
# Extract the gene names from the amplicons and create a heatmap for each gene
gene_names <- unique(gsub("_.+_.+", "", plotdata$amplicon_f))
for (gene in gene_names) {
  print(gene)
  gene_subset <- data.matrix(plotdata_all[grep(gene, row.names(plotdata_all)),])
  create_plots(gene_subset, gene, singlegene = TRUE, plot_type = "all")
}


