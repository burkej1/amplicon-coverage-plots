library(ggplot2)
library(readr)
library(gplots)
library(RColorBrewer) 
library(reshape2) 
library(heatmaply)

# Function to restructure a dataframe for heatmap plotting
to_heatmap_matrix <- function(df) {
  df_restructure <- dcast(df, amplicon_f~sample, fill=0)
  fix_rownames <- df_restructure[,-1]
  rownames(fix_rownames) <- df_restructure[,1]
  return(data.matrix(fix_rownames))
}

# Input file (should write to take a command line input)
amplicon_metrics <- read_delim("~/Documents/Code/amplicon-coverage-plots/amplicon_metrics.tsv", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)

headers <- c("chrom", "start", "end", "amplicon_f", "amplicon_r", "reads", "bases_c", "len", "pct", "sample")
colnames(amplicon_metrics) <- headers
plotdata <- amplicon_metrics[,c("amplicon_f", "sample", "reads")]  # Subset dataframe

plotdata_f <- to_heatmap_matrix(plotdata)  # All samples and all amplicons

plotdata_brca1 <- to_heatmap_matrix(plotdata[grep("BRCA1", plotdata$amplicon_f),])  # Subsetting to amplicons in a single gene
truncated <- replace(plotdata_brca1, plotdata_brca1 > 200, 200)  # Capping read values to make poorly performing amplicons clear

heatmaply(truncated, 
          fontsize_col = 8, fontsize_row = 8, 
          margins = c(180, 100), 
          file = "BRCA1_coverage.html"  # Save as HTML
          )



