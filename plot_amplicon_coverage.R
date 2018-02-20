library(ggplot2)  # General dependency for heatmaply/shinyheatmaply
library(readr)  # For importing tsvs
library(gplots)  # Used for something
library(RColorBrewer)  # Colour palettes
library(reshape2)  # For manipulating dataframes/matrices I think
library(heatmaply)  # Interactive heatmaps
library(scales)  # Rescaling scales
library(colorspace) # For generating colourspaces for rowside annotation
library(htmlwidgets)  # For saving non-self-contained plots


source("amplicon_coverage_functions.R")  # This may not work depending on the R working directory


# Input file (should rewrite to take a command line input)
amplicon_metrics <- read_delim("/Users/Jared/Documents/Project_Results_and_Data/BSTP/run1-7_amplicon_coverage_plotting/run1-7_amplicon_metrics_clean_190218.txt", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)

headers <- c("chrom", "start", "end", "amplicon_f", "amplicon_r", "reads", "bases_c", "len", "pct", "sample")
colnames(amplicon_metrics) <- headers
plotdata <- amplicon_metrics[,c("amplicon_f", "sample", "reads")]  # Subset dataframe

# # # # # # # # Overall plots # # # # # # # # 
# Creating matrices for overall plot
plotdata_all <- to_heatmap_matrix(plotdata, FALSE)  # Needed for other plot stuff
# full_hovertext <- generate_hoverframe(plotdata_all)  # Overlay matrix
# all_plot_colourscale <- create_plot_colourscale(max(plotdata_all), c(50, 100, 200, 500, 1000), viridis_pal()(6))


# Factored vectors of plate numbers (used for all heatmaps) and gene names (used for the overall heatmaps)
plate_vector <- data.frame(factor(gsub("^.+?_.+?_([0-9]+)-.+", "\\1", colnames(plotdata_all))))  # Coerced vector to df to add label
colnames(plate_vector) <- "Plate Number"
gene_names_vector <- data.frame(factor(gsub("_.+_.+", "", row.names(plotdata_all))))
colnames(gene_names_vector) <- "Gene"

# all_heatmap <- heatmaply(plotdata_all, 
#                          main = "All - Depth", 
#                          fontsize_col = 8,
#                          fontsize_row = 8,
#                          hide_colorbar = TRUE, 
#                          scale_fill_gradient_fun = all_plot_colourscale,
#                          col_side_colors = plate_vector, 
#                          row_side_colors = gene_names_vector, 
#                          file = "coverage_all.html",
#                          margins = c(200, 150),
#                          custom_hovertext = full_hovertext)
# 
# # Amplicon normalized all heatmap
# plotdata_all_amplicon_norm <- normalize_matrix(plotdata_all, by = "amplicon")
# plotdata_all_amplicon_norm_colourscale <- create_plot_colourscale(max(plotdata_all_amplicon_norm), 
#                                                                   breaks = c(0.5, 1.0, 1.5, 2.0), 
#                                                                   colours = viridis_pal()(5))
# all_heatmap_amplicon_norm <- heatmaply(plotdata_all_amplicon_norm, 
#                                        main = "All - Amplicon Normalised", 
#                                        fontsize_col = 8,
#                                        fontsize_row = 8,
#                                        hide_colorbar = TRUE, 
#                                        scale_fill_gradient_fun = plotdata_all_amplicon_norm_colourscale,
#                                        col_side_colors = plate_vector, 
#                                        row_side_colors = gene_names_vector, 
#                                        file = "amplicon_normalised_all.html",
#                                        margins = c(200, 150),
#                                        custom_hovertext = full_hovertext)
# 
# # Sample normalized all heatmap
# plotdata_all_sample_norm <- normalize_matrix(plotdata_all, by = "sample")
# plotdata_all_sample_norm_colourscale <- create_plot_colourscale(max(plotdata_all_sample_norm), 
#                                                                     breaks = c(0.5, 1.0, 1.5, 2.0), 
#                                                                     colours = viridis_pal()(5))
# all_heatmap_sample_norm <- heatmaply(plotdata_all_sample_norm, 
#                                      main = "All - Sample Normalised", 
#                                      fontsize_col = 8,
#                                      fontsize_row = 8,
#                                      hide_colorbar = TRUE, 
#                                      scale_fill_gradient_fun = plotdata_all_sample_norm_colourscale,
#                                      col_side_colors = plate_vector, 
#                                      row_side_colors = gene_names_vector, 
#                                      file = "sample_normalised_all.html",
#                                      margins = c(200, 150),
#                                      custom_hovertext = full_hovertext)
# 
# cat(paste("All - Sample Normalised", "0-50 | 50-100 | 100-200 | 200-500 | 500-1000 | 1000+", sep = "\n"))

# # # # # # # # Subsampled plots # # # # # # # # 
# General subsampling
plotdata_all_sampled <- plotdata_all[, sample(ncol(plotdata_all), 10)]  # Randomly choosing samples from the original dataframe
create_plots(plotdata_all_sampled, "subsample")

# Plate specific subsampling
unique_plates <- unique(plate_vector[,1])
for (plate in unique_plates) {
  plate_grep <- paste("_", plate, "-", sep="")
  plate_subset <- plotdata_all[,grep(plate_grep, colnames(plotdata_all))]
}

# # # # # # # # Gene plots # # # # # # # # 
# Extract the gene names from the amplicons and create a heatmap for each gene
gene_names <- unique(gsub("_.+_.+", "", plotdata$amplicon_f))
for (gene in gene_names) {
  print(gene)
  gene_subset <- data.matrix(plotdata_all[grep(gene, row.names(plotdata_all)),])

  # Factored vector of amplicon exons
  exon_vector <- data.frame(factor(gsub("^.+?_(.+)_.+", "\\1", row.names(gene_subset))))  # Coerced vector to df to add label
  colnames(exon_vector) <- "Exon Number"
  
  print("Depth")
  hovertext <- generate_hoverframe(gene_subset)
  plotname <- paste(gene, "_coverage.html", sep="")
  plot_colourscale <- create_plot_colourscale(max(gene_subset), c(50, 100, 200, 500, 1000), viridis_pal()(6))
  geneplot <- heatmaply(gene_subset, 
                        main = paste(gene, "- Depth"), 
                        hide_colorbar = TRUE,
                        fontsize_col = 8, 
                        fontsize_row = 8, 
                        scale_fill_gradient_fun = plot_colourscale,
                        col_side_colors = plate_vector,
                        row_side_colors = exon_vector,
                        plot_method = "ggplot", 
                        margins = c(200, 150), 
                        file = plotname,
                        custom_hovertext = hovertext)
  
  # # Normalised plots
  # # # Amplicon
  # print("Amplicon")
  # ampnorm_plotname <- paste(gene, "_amplicon-normalised_coverage.html", sep="")
  # ampnorm_matrix <- normalize_matrix(gene_subset, by = "amplicon")
  # ampnorm_colourscale <- create_plot_colourscale(max(ampnorm_matrix), breaks = c(0.5, 1.0, 1.5, 2.0), colours = viridis_pal()(5))
  # geneplot_ampnorm <- heatmaply(ampnorm_matrix, 
  #                               main = paste(gene, "- Amplicon Normalised"), 
  #                               hide_colorbar = TRUE,
  #                               fontsize_col = 8, 
  #                               fontsize_row = 8, 
  #                               scale_fill_gradient_fun = ampnorm_colourscale,
  #                               col_side_colors = plate_vector,
  #                               row_side_colors = exon_vector, 
  #                               margins = c(200, 150),
  #                               file = ampnorm_plotname,
  #                               custom_hovertext = hovertext)
  # # # Sample
  # print("Sample")
  # sampnorm_plotname <- paste(gene, "_sample-normalised_coverage.html", sep="")
  # sampnorm_matrix <- normalize_matrix(gene_subset, by = "sample")
  # sampnorm_colourscale <- create_plot_colourscale(max(sampnorm_matrix), breaks = c(0.5, 1.0, 1.5, 2.0), colours = viridis_pal()(5))
  # geneplot_sampnorm <- heatmaply(sampnorm_matrix, 
  #                                main = paste(gene, "- Sample Normalised"), 
  #                                hide_colorbar = TRUE,
  #                                fontsize_col = 8, 
  #                                fontsize_row = 8, 
  #                                scale_fill_gradient_fun = sampnorm_colourscale,
  #                                col_side_colors = plate_vector,
  #                                row_side_colors = exon_vector, 
  #                                margins = c(200, 150),  
  #                                file = sampnorm_plotname,
  #                                custom_hovertext = hovertext)
}


# # Testing
# plotdata_brca1 <- to_heatmap_matrix(plotdata[grep("BRCA1", plotdata$amplicon_f),], TRUE)  # Subsetting to amplicons in a single gene
# 
# plot <- heatmaply(plotdata_all, 
#           fontsize_col = 8, fontsize_row = 8, 
#           margins = c(180, 100), 
#           file = "BRCA1_coverage.html"  # Save as HTML
#           )
# 
# 
# plotdata_all_amplicon_norm <- normalize_matrix(plotdata_all, by = "sample")
# plotdata_all_sample_norm <- normalize_matrix(plotdata_all, by = "amplicon")
# 
# colourbreaks <- rescale(c(0, 0.5, 1, 2, 10), to = c(0, 1))
# colourscale <- brewer.pal(5, "Blues")
# all_heatmap <- heatmaply(plotdata_all_amplicon_norm,
#                          fontsize_col = 8,
#                          fontsize_row = 8,
#                          scale_fill_gradient_fun = scale_fill_gradientn(colors = c("red", "blue", "green", "orange", "purple"),
#                          values = colourbreaks),
#                          file = "all_coverage_amplicon_norm.html",
#                          custom_hovertext = full_hovertext)
# gene_subset_capped <- data.matrix(plotdata_all_capped[grep(gene, row.names(plotdata_all)),])
# plotdata_all_capped <- to_heatmap_matrix(plotdata, TRUE, lowvalue = 100, lowfloor = 0, highvalue = 1000)  # All samples and all amplicons, capped for visualisation

# library(graphics)
# colourpalette <- viridis_pal()(6)
# colourmatrix <- matrix(data = 1:6, nrow = 1, ncol = 6)
# image(z = colourmatrix, col = colourpalette)
