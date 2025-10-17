install.packages("readr")
library(readr)

data <- read_tsv("results/SRP181622_Symbols.tsv")
metadata <- read_tsv("data/SRP181622/metadata_SRP181622.tsv")

titles <- metadata$refinebio_title

install.packages("cluster")
install.packages("factoextra")
library(cluster)
library

install.packages("stringr")
library(stringr)

row_variances <- apply(data[, sapply(data, is.numerci)], 1, var)
sorted_indices <- order(row_variances, decreasing = TRUE)

most_variable_rows <- data[sorted_indices[1:5000], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]

pam_res <- pam(t(most_variable_rows), k = 2)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 1) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 2) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

samples <- metadata$refinebio_accession_code
cluster_data <- data.frame(Sample = samples, Cluster = cluster_assignments)
write.csv(cluster_data, "results/pam_cluster_data.csv", row.names = FALSE)

pam_res <- pam(t(most_variable_rows), k = 4)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 1) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 2) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 3) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 4) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

pam_res <- pam(t(most_variable_rows), k = 6)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 1) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 2) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 3) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 4) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 5) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 6) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

pam_res <- pam(t(most_variable_rows), k = 8)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 1) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 2) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 3) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 4) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 5) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 6) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 7) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 8) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

pam_res <- pam(t(most_variable_rows), k = 10)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 1) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 2) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 3) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 4) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 5) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 6) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 7) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 8) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 3) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 4) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 5) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 6) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 7) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 9) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 3) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 4) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 5) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 6) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 7) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 10) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

most_variable_rows <- data[sorted_indices[1:10000], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]

pam_res <- pam(t(most_variable_rows), k = 2)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 1) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 2) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

most_variable_rows <- data[sorted_indices[1:1000], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]

pam_res <- pam(t(most_variable_rows), k = 2)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 1) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 2) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

most_variable_rows <- data[sorted_indices[1:100], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]

pam_res <- pam(t(most_variable_rows), k = 2)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 1) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 2) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

most_variable_rows <- data[sorted_indices[1:10], ]
most_variable_rows <- most_variable_rows[, -1]
most_variable_rows <- most_variable_rows[, -1]

pam_res <- pam(t(most_variable_rows), k = 2)
cluster_assignments <- pam_res$clustering
print(cluster_assignments)
cluster_sizes <- table(pam_res$clustering)
print(cluster_sizes)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 1) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)

count <- 0
for (i in seq_along(titles)) {
  det_res <- str_detect(titles[i], control)
  if (cluster_assignments[i] == 2) {
    if (det_res == TRUE) {
      count <- count + 1
    }
  }
}
print(count)