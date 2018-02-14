library(ggplot2)  # General dependency for heatmaply/shinyheatmaply
library(readr)  # For importing tsvs
library(gplots)  # Used for something
library(RColorBrewer)  # Colour palettes
library(reshape2)  # For manipulating dataframes/matrices I think
library(heatmaply)  # Interactive heatmaps
library(scales)  # Rescaling scales


# Input file (should write to take a command line input)
amplicon_metrics <- read_delim("/Users/Jared/Documents/Code/amplicon-coverage-plots/run1-7_amplicon_metrics_all.tsv", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)

headers <- c("chrom", "start", "end", "amplicon_f", "amplicon_r", "reads", "bases_c", "len", "pct", "sample")
colnames(amplicon_metrics) <- headers
plotdata <- amplicon_metrics[,c("amplicon_f", "sample", "reads")]  # Subset dataframe
# plotdata <- plotdata[1:2386,]  # Subsetting for testing purposes

# Creating matrices
# plotdata_all_capped <- to_heatmap_matrix(plotdata, TRUE, lowvalue = 100, lowfloor = 0, highvalue = 1000)  # All samples and all amplicons, capped for visualisation
plotdata_all <- to_heatmap_matrix(plotdata, FALSE)  # Uncapped amplicon data, use this for overlay matrices
full_hovertext <- generate_hoverframe(plotdata_all)  # Overlay matrix
all_plot_colourscale <- create_plot_colourscale(max(plotdata_all))

all_heatmap <- heatmaply(plotdata_all,
                         fontsize_col = 8,
                         fontsize_row = 8,
                         scale_fill_gradient_fun = all_plot_colourscale, 
                         file = "all_coverage.html",
                         margins = c(200, 150), 
                         custom_hovertext = full_hovertext)


# Extract the gene names from the amplicons and create a heatmap for each gene
gene_names <- unique(gsub("_.+_.+", "", plotdata$amplicon_f))
for (gene in gene_names) {
  gene_subset <- data.matrix(plotdata_all[grep(gene, row.names(plotdata_all)),])
  # gene_subset_capped <- data.matrix(plotdata_all_capped[grep(gene, row.names(plotdata_all)),])
  hovertext <- generate_hoverframe(gene_subset)
  plotname <- paste(gene, "_coverage.html", sep="")
  plot_colourscale <- create_plot_colourscale(max(gene_subset))
  geneplot <- heatmaply(gene_subset, 
                        fontsize_col = 8, 
                        fontsize_row = 8, 
                        scale_fill_gradient_fun = plot_colourscale,
                        margins = c(200, 150), 
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


plotdata_all_amplicon_norm <- normalize_matrix(plotdata_all, by = "sample")
plotdata_all_sample_norm <- normalize_matrix(plotdata_all, by = "amplicon")

colourbreaks <- rescale(c(0, 0.5, 1, 2, 10), to = c(0, 1))
colourscale <- brewer.pal(5, "Blues")
all_heatmap <- heatmaply(plotdata_all_amplicon_norm,
                         fontsize_col = 8,
                         fontsize_row = 8,
                         scale_fill_gradient_fun = scale_fill_gradientn(colors = c("red", "blue", "green", "orange", "purple"),
                         values = colourbreaks),
                         file = "all_coverage_amplicon_norm.html",
                         custom_hovertext = full_hovertext)

