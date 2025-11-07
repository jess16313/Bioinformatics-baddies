# loading and set up and predicting groups from assignment 1

install.packages("readr")
library(readr)

data <- read_tsv("results/SRP181622_Symbols.tsv")
metadata <- read_tsv("data/SRP181622/metadata_SRP181622.tsv")

install.packages("cluster")
install.packages("factoextra")
library(cluster)
library

install.packages("stringr")
library(stringr)

row_variances <- apply(data[, sapply(data, is.numeric)], 1, var)
sorted_indices <- order(row_variances, decreasing = TRUE)

most_variable_rows <- data[sorted_indices[1:5000], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- t(most_variable_rows)

install.packages('e1071') 
install.packages('caTools')
install.packages('ggplot2')
install.packages('caret')

library(caret)
library(e1071) 
library(caTools)
library(ggplot2)

install.packages("pROC")
library(pROC)

titles <- metadata$refinebio_title
samples <- metadata$refinebio_accession_code

for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (det_res == TRUE) {
    titles[i] <- "Control"
  } else {
    titles[i] <- "Stressed"
  }
}

set.seed(1234)

numTitles <- as.numeric(factor(titles, levels = c("Control", "Stressed"), labels = c(0, 1)))

scaleData <- scale(most_variable_rows)

for (i in seq_along(numTitles)) {
  if (numTitles[i] == 1) {
    numTitles[i] <- 0
  } else {
    numTitles[i] <- 1
  }
}

split <- sample.split(numTitles, SplitRatio = 0.8)
training_set <- subset(scaleData, split == TRUE)
test_set <- subset(scaleData, split == FALSE)

classifier <- svm(numTitles ~ ., 
                  data = scaleData, 
                  type = 'C-classification', 
                  kernel = 'radial', 
                  gamma = 0.1,
                  probability = TRUE)

pred <- predict(classifier, newdata = test_set)
pred_probs<- predict(classifier, newdata = test_set, probability = TRUE)

write.csv(pred, "dataForHeatmapAZAssign4.csv")

probabilities <- attr(pred_probs, "probabilities")

test_samples <- rownames(test_set)

actual_classes <- vector()

for (i in seq_along(test_samples)) {
  for (k in seq_along(samples)) {
    if (test_samples[i] == samples[k]) {
    actual_classes[i] <- numTitles[k]
    }
  }
}

pred
actual_classes
table(actual_classes, pred)

write.csv(actual_classes, "actualAZ.csv")

predictions <- as.numeric(pred)
for (i in seq_along(predictions)) {
  if (predictions[i] == 1) {
    predictions[i] <- 0
  } else {
    predictions[i] <- 1
  }
}

probs <- vector()
for (i in seq_along(test_samples)) {
  probs[i] <- probabilities[i, predictions[i]+1]
}

roc_obj <- roc(response = actual_classes, predictor = probs)
auc_value <- auc(roc_obj)
print(auc_value)

sample_model_matrix <- cbind(test_samples, predictions)

# predicting clusters from assignment 3

pam_res <- pam(most_variable_rows, k = 2)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

for (i in seq_along(cluster_assignments)) {
  if (cluster_assignments[i] == 1) {
    cluster_assignments[i] <- 0
  } else {
    cluster_assignments[i] <- 1
  }
}

clusters = as.factor(cluster_assignments)

split <- sample.split(clusters, SplitRatio = 0.8)
training_set <- subset(scaleData, split == TRUE)
test_set <- subset(scaleData, split == FALSE)

classifier <- svm(numTitles ~ ., 
                  data = scaleData, 
                  type = 'C-classification', 
                  kernel = 'radial', 
                  gamma = 0.1, 
                  probability = TRUE)

pred <- predict(classifier, newdata = test_set)
pred_probs<- predict(classifier, newdata = test_set, probability = TRUE)

probabilities <- attr(pred_probs, "probabilities")

test_samples <- rownames(test_set)

actual_clusters <- vector()

for (i in seq_along(test_samples)) {
  for (k in seq_along(samples)) {
    if (test_samples[i] == samples[k]) {
      actual_clusters[i] <- cluster_assignments[k]
    }
  }
}

pred
actual_clusters
table(actual_clusters, pred)

predictions <- as.numeric(pred)
for (i in seq_along(predictions)) {
  if (predictions[i] == 1) {
    predictions[i] <- 0
  } else {
    predictions[i] <- 1
  }
}

probs <- vector()
for (i in seq_along(test_samples)) {
  probs[i] <- probabilities[i, predictions[i]+1]
}

roc_obj <- roc(response = actual_clusters, predictor = probs)
auc_value <- auc(roc_obj)
print(auc_value)

sample_model_matrix <- cbind(sample_model_matrix, test_samples, predictions)

# retraining with different numbers of genes

# groups from assignment 1

# 10 genes

most_variable_rows <- data[sorted_indices[1:10], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- t(most_variable_rows)
scaleData <- scale(most_variable_rows)

split <- sample.split(numTitles, SplitRatio = 0.8)
training_set <- subset(scaleData, split == TRUE)
test_set <- subset(scaleData, split == FALSE)

classifier <- svm(numTitles ~ ., 
                  data = scaleData, 
                  type = 'C-classification', 
                  kernel = 'radial', 
                  gamma = 0.1,
                  probability = TRUE)

pred <- predict(classifier, newdata = test_set)
pred_probs<- predict(classifier, newdata = test_set, probability = TRUE)

probabilities <- attr(pred_probs, "probabilities")

test_samples <- rownames(test_set)

actual_classes <- vector()

for (i in seq_along(test_samples)) {
  for (k in seq_along(samples)) {
    if (test_samples[i] == samples[k]) {
      actual_classes[i] <- numTitles[k]
    }
  }
}

pred
actual_classes
table(actual_classes, pred)

predictions <- as.numeric(pred)
for (i in seq_along(predictions)) {
  if (predictions[i] == 1) {
    predictions[i] <- 0
  } else {
    predictions[i] <- 1
  }
}

probs <- vector()
for (i in seq_along(test_samples)) {
  probs[i] <- probabilities[i, predictions[i]+1]
}

roc_obj <- roc(response = actual_classes, predictor = probs)
auc_value <- auc(roc_obj)
print(auc_value)

sample_model_matrix <- cbind(sample_model_matrix, test_samples, predictions)

# 100 genes 

most_variable_rows <- data[sorted_indices[1:100], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- t(most_variable_rows)
scaleData <- scale(most_variable_rows)

training_set <- subset(scaleData, split == TRUE)
test_set <- subset(scaleData, split == FALSE)

classifier <- svm(numTitles ~ ., 
                  data = scaleData, 
                  type = 'C-classification', 
                  kernel = 'radial', 
                  gamma = 0.1,
                  probability = TRUE)

pred <- predict(classifier, newdata = test_set)
pred_probs<- predict(classifier, newdata = test_set, probability = TRUE)

probabilities <- attr(pred_probs, "probabilities")

test_samples <- rownames(test_set)

actual_classes <- vector()

for (i in seq_along(test_samples)) {
  for (k in seq_along(samples)) {
    if (test_samples[i] == samples[k]) {
      actual_classes[i] <- numTitles[k]
    }
  }
}

pred
actual_classes
table(actual_classes, pred)

predictions <- as.numeric(pred)
for (i in seq_along(predictions)) {
  if (predictions[i] == 1) {
    predictions[i] <- 0
  } else {
    predictions[i] <- 1
  }
}

probs <- vector()
for (i in seq_along(test_samples)) {
  probs[i] <- probabilities[i, predictions[i]+1]
}

roc_obj <- roc(response = actual_classes, predictor = probs)
auc_value <- auc(roc_obj)
print(auc_value)

sample_model_matrix <- cbind(sample_model_matrix, test_samples, predictions)

# 1000 genes

most_variable_rows <- data[sorted_indices[1:1000], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- t(most_variable_rows)
scaleData <- scale(most_variable_rows)

training_set <- subset(scaleData, split == TRUE)
test_set <- subset(scaleData, split == FALSE)

classifier <- svm(numTitles ~ ., 
                  data = scaleData, 
                  type = 'C-classification', 
                  kernel = 'radial', 
                  gamma = 0.1,
                  probability = TRUE)

pred <- predict(classifier, newdata = test_set)
pred_probs<- predict(classifier, newdata = test_set, probability = TRUE)

probabilities <- attr(pred_probs, "probabilities")

test_samples <- rownames(test_set)

actual_classes <- vector()

for (i in seq_along(test_samples)) {
  for (k in seq_along(samples)) {
    if (test_samples[i] == samples[k]) {
      actual_classes[i] <- numTitles[k]
    }
  }
}

pred
actual_classes
table(actual_classes, pred)

predictions <- as.numeric(pred)
for (i in seq_along(predictions)) {
  if (predictions[i] == 1) {
    predictions[i] <- 0
  } else {
    predictions[i] <- 1
  }
}

probs <- vector()
for (i in seq_along(test_samples)) {
  probs[i] <- probabilities[i, predictions[i]+1]
}

roc_obj <- roc(response = actual_classes, predictor = probs)
auc_value <- auc(roc_obj)
print(auc_value)

sample_model_matrix <- cbind(sample_model_matrix, test_samples, predictions)

# 10000 genes

most_variable_rows <- data[sorted_indices[1:10000], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- t(most_variable_rows)
scaleData <- scale(most_variable_rows)

training_set <- subset(scaleData, split == TRUE)
test_set <- subset(scaleData, split == FALSE)

classifier <- svm(numTitles ~ ., 
                  data = scaleData, 
                  type = 'C-classification', 
                  kernel = 'radial', 
                  gamma = 0.1,
                  probability = TRUE)

pred <- predict(classifier, newdata = test_set)
pred_probs<- predict(classifier, newdata = test_set, probability = TRUE)

probabilities <- attr(pred_probs, "probabilities")

test_samples <- rownames(test_set)

actual_classes <- vector()

for (i in seq_along(test_samples)) {
  for (k in seq_along(samples)) {
    if (test_samples[i] == samples[k]) {
      actual_classes[i] <- numTitles[k]
    }
  }
}

pred
actual_classes
table(actual_classes, pred)

predictions <- as.numeric(pred)
for (i in seq_along(predictions)) {
  if (predictions[i] == 1) {
    predictions[i] <- 0
  } else {
    predictions[i] <- 1
  }
}

probs <- vector()
for (i in seq_along(test_samples)) {
  probs[i] <- probabilities[i, predictions[i]+1]
}

roc_obj <- roc(response = actual_classes, predictor = probs)
auc_value <- auc(roc_obj)
print(auc_value)

sample_model_matrix <- cbind(sample_model_matrix, test_samples, predictions)

# clusters from assignment 3

# 10 genes

most_variable_rows <- data[sorted_indices[1:10], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- t(most_variable_rows)
scaleData <- scale(most_variable_rows)

pam_res <- pam(most_variable_rows, k = 2)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

for (i in seq_along(cluster_assignments)) {
  if (cluster_assignments[i] == 1) {
    cluster_assignments[i] <- 0
  } else {
    cluster_assignments[i] <- 1
  }
}

clusters = as.factor(cluster_assignments)

split <- sample.split(clusters, SplitRatio = 0.8)
training_set <- subset(scaleData, split == TRUE)
test_set <- subset(scaleData, split == FALSE)

classifier <- svm(numTitles ~ ., 
                  data = scaleData, 
                  type = 'C-classification', 
                  kernel = 'radial', 
                  gamma = 0.1, 
                  probability = TRUE)

pred <- predict(classifier, newdata = test_set)
pred_probs<- predict(classifier, newdata = test_set, probability = TRUE)

probabilities <- attr(pred_probs, "probabilities")

test_samples <- rownames(test_set)

actual_clusters <- vector()

for (i in seq_along(test_samples)) {
  for (k in seq_along(samples)) {
    if (test_samples[i] == samples[k]) {
      actual_clusters[i] <- cluster_assignments[k]
    }
  }
}

pred
actual_clusters
table(actual_clusters, pred)

predictions <- as.numeric(pred)
for (i in seq_along(predictions)) {
  if (predictions[i] == 1) {
    predictions[i] <- 0
  } else {
    predictions[i] <- 1
  }
}

probs <- vector()
for (i in seq_along(test_samples)) {
  probs[i] <- probabilities[i, predictions[i]+1]
}

roc_obj <- roc(response = actual_clusters, predictor = probs)
auc_value <- auc(roc_obj)
print(auc_value)

sample_model_matrix <- cbind(sample_model_matrix, test_samples, predictions)

# 100 genes

most_variable_rows <- data[sorted_indices[1:100], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- t(most_variable_rows)
scaleData <- scale(most_variable_rows)

pam_res <- pam(most_variable_rows, k = 2)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

for (i in seq_along(cluster_assignments)) {
  if (cluster_assignments[i] == 1) {
    cluster_assignments[i] <- 0
  } else {
    cluster_assignments[i] <- 1
  }
}

clusters = as.factor(cluster_assignments)

split <- sample.split(clusters, SplitRatio = 0.8)
training_set <- subset(scaleData, split == TRUE)
test_set <- subset(scaleData, split == FALSE)

classifier <- svm(numTitles ~ ., 
                  data = scaleData, 
                  type = 'C-classification', 
                  kernel = 'radial', 
                  gamma = 0.1, 
                  probability = TRUE)

pred <- predict(classifier, newdata = test_set)
pred_probs<- predict(classifier, newdata = test_set, probability = TRUE)

probabilities <- attr(pred_probs, "probabilities")

test_samples <- rownames(test_set)

actual_clusters <- vector()

for (i in seq_along(test_samples)) {
  for (k in seq_along(samples)) {
    if (test_samples[i] == samples[k]) {
      actual_clusters[i] <- cluster_assignments[k]
    }
  }
}

pred
actual_clusters
table(actual_clusters, pred)

predictions <- as.numeric(pred)
for (i in seq_along(predictions)) {
  if (predictions[i] == 1) {
    predictions[i] <- 0
  } else {
    predictions[i] <- 1
  }
}

probs <- vector()
for (i in seq_along(test_samples)) {
  probs[i] <- probabilities[i, predictions[i]+1]
}

roc_obj <- roc(response = actual_clusters, predictor = probs)
auc_value <- auc(roc_obj)
print(auc_value)

sample_model_matrix <- cbind(sample_model_matrix, test_samples, predictions)

# 1000 genes

most_variable_rows <- data[sorted_indices[1:1000], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- t(most_variable_rows)
scaleData <- scale(most_variable_rows)

pam_res <- pam(most_variable_rows, k = 2)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

for (i in seq_along(cluster_assignments)) {
  if (cluster_assignments[i] == 1) {
    cluster_assignments[i] <- 0
  } else {
    cluster_assignments[i] <- 1
  }
}

clusters = as.factor(cluster_assignments)

split <- sample.split(clusters, SplitRatio = 0.8)
training_set <- subset(scaleData, split == TRUE)
test_set <- subset(scaleData, split == FALSE)

classifier <- svm(numTitles ~ ., 
                  data = scaleData, 
                  type = 'C-classification', 
                  kernel = 'radial', 
                  gamma = 0.1, 
                  probability = TRUE)

pred <- predict(classifier, newdata = test_set)
pred_probs<- predict(classifier, newdata = test_set, probability = TRUE)

probabilities <- attr(pred_probs, "probabilities")

test_samples <- rownames(test_set)

actual_clusters <- vector()

for (i in seq_along(test_samples)) {
  for (k in seq_along(samples)) {
    if (test_samples[i] == samples[k]) {
      actual_clusters[i] <- cluster_assignments[k]
    }
  }
}

pred
actual_clusters
table(actual_clusters, pred)

predictions <- as.numeric(pred)
for (i in seq_along(predictions)) {
  if (predictions[i] == 1) {
    predictions[i] <- 0
  } else {
    predictions[i] <- 1
  }
}

probs <- vector()
for (i in seq_along(test_samples)) {
  probs[i] <- probabilities[i, predictions[i]+1]
}

roc_obj <- roc(response = actual_clusters, predictor = probs)
auc_value <- auc(roc_obj)
print(auc_value)

sample_model_matrix <- cbind(sample_model_matrix, test_samples, predictions)

# 10000 genes

most_variable_rows <- data[sorted_indices[1:10000], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- t(most_variable_rows)
scaleData <- scale(most_variable_rows)

pam_res <- pam(most_variable_rows, k = 2)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

for (i in seq_along(cluster_assignments)) {
  if (cluster_assignments[i] == 1) {
    cluster_assignments[i] <- 0
  } else {
    cluster_assignments[i] <- 1
  }
}

clusters = as.factor(cluster_assignments)

split <- sample.split(clusters, SplitRatio = 0.8)
training_set <- subset(scaleData, split == TRUE)
test_set <- subset(scaleData, split == FALSE)

classifier <- svm(numTitles ~ ., 
                  data = scaleData, 
                  type = 'C-classification', 
                  kernel = 'radial', 
                  gamma = 0.1, 
                  probability = TRUE)

pred <- predict(classifier, newdata = test_set)
pred_probs<- predict(classifier, newdata = test_set, probability = TRUE)

probabilities <- attr(pred_probs, "probabilities")

test_samples <- rownames(test_set)

actual_clusters <- vector()

for (i in seq_along(test_samples)) {
  for (k in seq_along(samples)) {
    if (test_samples[i] == samples[k]) {
      actual_clusters[i] <- cluster_assignments[k]
    }
  }
}

pred
actual_clusters
table(actual_clusters, pred)

predictions <- as.numeric(pred)
for (i in seq_along(predictions)) {
  if (predictions[i] == 1) {
    predictions[i] <- 0
  } else {
    predictions[i] <- 1
  }
}

probs <- vector()
for (i in seq_along(test_samples)) {
  probs[i] <- probabilities[i, predictions[i]+1]
}

roc_obj <- roc(response = actual_clusters, predictor = probs)
auc_value <- auc(roc_obj)
print(auc_value)

sample_model_matrix <- cbind(sample_model_matrix, test_samples, predictions)

write.csv(sample_model_matrix, "AZSampleModelMatrix.csv")

