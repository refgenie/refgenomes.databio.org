setwd("/home/claudehu/Desktop/trials/scatlas_embedding_eval")
# equivalence to /project/shefflab/brickyard/results_analysis/scatlas/embedding_eval

library(GenomicDistributionsData)
library(GenomicDistributions)
library(LOLA)
library(ggplot2)
library(reshape2)
library(gtools)
library(pheatmap)
library(viridis)


plotSignalHeatMap = function(signal_stats_df){
  # Add rownames as a column
  signal_stats_df$cluster = rownames(signal_stats_df)
  
  # Reshape to long format
  df_long = reshape2::melt(signal_stats_df, id.vars = "cluster")
  colnames(df_long) = c("cluster", "cell", "value")
  
  df_long$cluster = factor(df_long$cluster, levels = mixedsort(unique(df_long$cluster)))
  
  ggplot(df_long, aes(x = cell, y = cluster, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(size = 8)
    )
}

plotSignalPHeatMap = function(signal_stats_df, out_path = NULL){
  # Convert to matrix
  mat <- as.matrix(signal_stats_df)
  
  # Add rownames if missing
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("row_", seq_len(nrow(mat)))
  }
  
  # Use pheatmap for hierarchical clustering and dendrograms
  pheatmap::pheatmap(
    mat,
    scale = "none",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    color = viridis::viridis(100),
    angle_col = 90,
    fontsize_row = 8,
    fontsize_col = 8,
    treeheight_row = 30,
    treeheight_col = 30,
    border_color = NA,
    filename = out_path,
    width = 10,
    height = 12,
    dpi = 300
  )
}

signalStatsList = function(bed_folder, genome_signal, key_stat="median"){
  # hg38_sig = openSignalMatrix_hg38()
  cluster_bed_paths = list.files(path=bed_folder, full.names = TRUE)
  
  median_dfs = list()
  labels = c()
    
  for (bed in cluster_bed_paths){
    # read bed cluster and get signal summary
    
    regions = readBed(bed)
    cluster_signal_summary = calcSummarySignal(regions, genome_signal)
    signal_medians = cluster_signal_summary$matrixStats[key_stat, , drop=FALSE]
    
    # extract medians and label with cluster name
    cluster_name = tools::file_path_sans_ext(basename(bed))
    rownames(signal_medians) = cluster_name
    
    median_dfs[[cluster_name]] = signal_medians
    labels = c(labels, cluster_name)
  }
  return(median_dfs)
  # final_df = do.call(rbind, median_dfs)
  # final_df = final_df[mixedsort(rownames(final_df)), ]
  # return(final_df)
}

# clusters_med_signal = signalStatsForHeatMap("leiden_cluster_beds")
hg38_sig = openSignalMatrix_hg38()
clusters_med_signals = signalStatsList("kmeans_73_clusters_beds", hg38_sig)
clusters_med_df = do.call(rbind, clusters_med_signals)
clusters_med_df = clusters_med_df[mixedsort(rownames(clusters_med_df)), ]


png("/home/claudehu/Downloads/kmeans_cell_signal_median_heatmap.png", width = 8, height = 9, units = "in", res = 300)
plotSignalHeatMap(clusters_med_df)
dev.off()


plotSignalPHeatMap(clusters_med_df, "/home/claudehu/Downloads/kmeans_cell_signal_median_pheatmap.png")

