load("C:/Users/Jess/Downloads/Bioinformatics Baddies/project_assignment3workspace.RData")
savehistory("analysis_history.Rhistory")
cd ..
install.packages(RandomForest)
install.packages(randomForest)
install.packages("randomForest")
library(randomForest)
data_for_rf <- as.data.frame(t(data_5000))
load("C:/Users/Jess/Downloads/Bioinformatics Baddies/project_assignment3workspace.RData")
data_for_rf <- as.data.frame(t(data_5000))
data_for_rf$sample_id <- rownames(data_for_rf)
data_for_rf <- data_for_rf %>% mutate(treatment = ifelse(grepl("Control", sample_id), "Control", "Stressed"))
library(tidymodels)
data_for_rf <- data_for_rf %>% mutate(treatment = ifelse(grepl("Control", sample_id), "Control", "Stressed"))
view(data_for_rf)
View(data_for_rf)
View(data_for_clustering)
rf_data <- data_for_rf %>%
  select(-sample_id, -refinebio_title, -refinebio_specimen_part) %>%  # Remove metadata columns
  select(treatment, everything())
colnames(data_for_rf) <- make.names(colnames(data_for_rf))
gene_names <- colnames(data_for_rf)[1:(ncol(data_for_rf)-3)]  # Exclude sample_id, treatment, tissue
colnames(data_for_rf)[1:(ncol(data_for_rf)-3)] <- paste0("gene_", 1:length(gene_names))
data_for_rf <- data_for_rf %>%
  left_join(
    metadata %>% select(refinebio_accession_code, refinebio_title, refinebio_specimen_part),
    by = c("sample_id" = "refinebio_accession_code")
  )
data_for_rf <- data_for_rf %>%
  mutate(
    treatment = factor(ifelse(grepl("Control", refinebio_title), "Control", "Stressed")),
    tissue = factor(refinebio_specimen_part)
  )
rf_data <- data_for_rf %>%
  select(-sample_id, -refinebio_title, -refinebio_specimen_part) %>%  # Remove metadata columns
  select(treatment, everything())
cat("Final RF data dimensions:", dim(rf_data), "\n")
cat("Outcome variable:\n")
print(table(rf_data$treatment))
set.seed(1223)
set.seed(123)
train_indices <- sample(1:nrow(data_for_rf), size = 0.7 * nrow(data_for_rf))
train_data <- data_for_rf[train_indices, ]
test_data <- data_for_rf[-train_indices, ]
rf_model <- randomForest(
  treatment ~ .,           # Predict treatment using all genes
  data = train_data,
  ntree = 500,            # Number of trees
  mtry = sqrt(ncol(train_data) - 1),  # Number of variables tried at each split
  importance = TRUE,       # Calculate variable importance
  proximity = TRUE         # Calculate sample proximity
)
print(rf_model)
predictions <- predict(rf_model, test_data)
confusion_matrix <- table(Predicted = predictions, Actual = test_data$treatment)
print(confusion_matrix)
data_two_rf <- readr("genes_5000/genes_5000.k=4.consensusClass.csv")
data_two_rf <- read.csv("genes_5000/genes_5000.k=4.consensusClass.csv")
data_two_rf <- as.data.frame(t(data_two_rf))
colnames(data_two_rf) <- c("sample_id", "cluster")
#oops
rf_two <- readr("genes_5000/genes_5000.k=4.consensusClass.csv", header=FALSE)
rf_two <- read.csv("genes_5000/genes_5000.k=4.consensusClass.csv", header=FALSE)
save()
save.image("project_assignment4workspace")
library(randomForest)
library(tidyverse)
load("C:/Users/Jess/Downloads/Bioinformatics Baddies/project_assignment4workspace.RData")
View(rf_model)
predictions <- predict(rf_model, test_data)
conf_matrix <- table(Predicted = predictions, Actual = test_data$treatment)
predictions <- predict(rf_model, test_data)
confusion_matrix <- table(Predicted = predictions, Actual = test_data$cluster)
std_label <- function(x) {
  x0 <- trimws(tolower(as.character(x)))
  # common mappings; add any your team used
  x0 <- gsub("[^a-z0-9]+", "_", x0)        # normalize punctuation/spaces
  d <- c(
    "control"="Control","ctrl"="Control","ctl"="Control","unstressed"="Control",
    "vehicle"="Control","sham"="Control",
    "stress"="Stress","stressed"="Stress","treated"="Stress","treatment"="Stress"
  )
  out <- ifelse(x0 %in% names(d), d[x0], x0)
  # fall back: capitalize first letter
  out <- ifelse(out %in% c("control","stress"), tools::toTitleCase(out), out)
  factor(out, levels=c("Control","Stress"))
}
preds_A <- read.csv("gabe_prediction_without_actual.csv")   # sample_id, pred_class
preds_B <- read.csv("dataForHeatmapAZAssign4.csv")
clean_id <- function(v) trimws(gsub("\\s+", "", v))
norm_id <- function(v) make.names(clean_id(v), unique=TRUE)
# preds_A$sample_id <- norm_id(preds_A$sample_id)
# preds_B$sample_id <- norm_id(preds_B$sample_id)
# preds_A$pred_class <- std_label(preds_A$pred_class)
# preds_B$pred_class <- std_label(preds_B$pred_class)
confusion_matrix$sample_id <-norm_id(confusion_matrix$sample_id)
rf_model$samle_id <- nrom_id(rf_model$sample_id)
rf_model$samle_id <- norm_id(rf_model$sample_id)
library(pheatmap)
# expr: samples x genes (you already have this)
genes_in <- intersect(colnames(expr), all_sig)
You said:
  my_sig <- unique(unlist(lapply(rf_results, `[[`, "genes")))
all_sig <- unique(my_sig)
library(pheatmap)
# expr: samples x genes (you already have this)
genes_in <- intersect(colnames(expr), all_sig)
mat <- expr[, genes_in, drop=FALSE]
# (Speed tip) If union is huge, cap to e.g. top 500 by variance:
if (ncol(mat) > 800) {
  v <- apply(mat, 2, var)
  genes_keep <- names(sort(v, decreasing=TRUE))[1:800]
  mat <- mat[, genes_keep, drop=FALSE]
}
# z-score genes for visualization
mat_scaled <- scale(mat)
# Sidebar annotation (pick ONE: A1 groups OR clusters)
# Using clusters you already built:
ann <- data.frame(Cluster = labels[rownames(mat_scaled)])
rownames(ann) <- rownames(mat_scaled)
# If using Assignment 1 Control/Stress instead:
# ann <- merge(data.frame(sample_id=rownames(mat_scaled)),
#              labels_A1[,c("sample_id","group")], by="sample_id", all.x=TRUE)
# rownames(ann) <- ann$sample_id; ann$sample_id <- NULL
pheatmap(
  t(mat_scaled),
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = TRUE,    # row dendrogram (genes)
  cluster_cols = TRUE,    # column dendrogram (samples)
  clustering_method = "complete",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  annotation_col = ann,   # sidebar on columns (samples)
  main = "Predictive Gene Signatures Heatmap",
  legend = TRUE,
  border_color = NA
)
save.image(file="proj_assignment4workspace.RData")
