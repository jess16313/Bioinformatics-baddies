library(org.Mm.eg.db)
library(AnnotationDbi)
library(ConsensusClusterPlus)
library(pheatmap)

expr_df <- read.delim("CGS4144_DataSet_Reconfig/data/SRP181622/SRP181622.tsv", check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_df) <- expr_df[[1]]
expr_df <- expr_df[, -1, drop = FALSE]
expr <- data.matrix(expr_df)
rownames(expr) <- gsub("\\..*", "", rownames(expr))
ens <- rownames(expr)
sym <- mapIds(org.Mm.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(expr) <- ifelse(is.na(sym) | sym == "", ens, sym)
expr <- aggregate(expr, by = list(Gene = rownames(expr)), FUN = mean)
rownames(expr) <- expr$Gene
expr$Gene <- NULL

mad_values <- apply(expr, 1, mad, na.rm = TRUE)
mad_values[is.na(mad_values)] <- 0
top5000 <- names(sort(mad_values, decreasing = TRUE))[seq_len(min(5000, nrow(expr)))]
data_5000 <- data.matrix(expr[top5000, , drop = FALSE])
write.csv(data_5000, "data_5000.csv", row.names = TRUE)

res <- ConsensusClusterPlus(d = data_5000, maxK = 4, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg = "km", distance = "euclidean", seed = 123, plot = "none", title = "consensus_k4")
cl <- res[[4]]$consensusClass
write.csv(data.frame(Sample = names(cl), Cluster = as.integer(cl)), "consensus_clusters_k4.csv", row.names = FALSE)

paths <- c(KMeans = "kmeans_clusters.csv", PAM = "pam_cluster_data.csv", Consensus = "consensus_clusters_k4.csv", Gabe = "gabe_cluster_with_samplenames.csv")
read_two <- function(path, label){df <- tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE), error = function(e) read.csv(path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)); keep <- vapply(df, function(x) any(!is.na(x) & nzchar(trimws(as.character(x)))), logical(1)); df <- df[, keep, drop = FALSE]; stopifnot(ncol(df) >= 2); out <- data.frame(Sample = trimws(as.character(df[[1]])), tmp = trimws(as.character(df[[2]])), stringsAsFactors = FALSE); names(out)[2] <- label; out[nzchar(out$Sample), , drop = FALSE]}
lst <- lapply(names(paths), function(m) read_two(paths[[m]], m)); names(lst) <- names(paths)
ann <- Reduce(function(x, y) merge(x, y, by = "Sample", all = FALSE), lst)
names(ann)[names(ann) == "Gabe"] <- "Gaussian"
write.csv(ann, "combined_cluster_annotations_FINAL.csv", row.names = FALSE)

expr <- as.matrix(read.csv("data_5000.csv", row.names = 1, check.names = FALSE))
ann <- read.csv("combined_cluster_annotations_FINAL.csv", row.names = 1, check.names = FALSE)
common <- intersect(colnames(expr), rownames(ann))
expr <- expr[, common, drop = FALSE]
ann <- ann[common, , drop = FALSE]
pheatmap(expr, scale = "row", color = colorRampPalette(c("navy","white","firebrick3"))(50), annotation_col = ann, show_rownames = FALSE, show_colnames = FALSE, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2", fontsize = 10, border_color = NA, main = "Heatmap of 5000 Most Variable Genes with Cluster Annotations")
pdf("heatmap_5000genes_clusters.pdf", width = 12, height = 10); pheatmap(expr, scale = "row", color = colorRampPalette(c("navy","white","firebrick3"))(50), annotation_col = ann, show_rownames = FALSE, show_colnames = FALSE, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2", fontsize = 10, border_color = NA, main = "Heatmap of 5000 Most Variable Genes with Cluster Annotations"); dev.off()
