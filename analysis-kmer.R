library(dplyr)
library(ranger)
library(hmeasure)
library(kernlab)
library(biogram)

set.seed(2017)

holdout_perc <- 0.2

seqs <- c(lapply(1L:50, function(ith_seq) 
  sample(c("A", "B", "C", "D"), size = sample(21L:30, 1), replace = TRUE)),
  lapply(1L:50, function(ith_seq) 
    sample(c("A", "B", "C", "D"), size = sample(21L:30, 1), replace = TRUE, 
           prob = c(0.35, 0.25, 0.25, 0.15)))) %>% 
  list2matrix()

kmers <- count_ngrams(seqs, 3, c("A", "B", "C", "D"), pos = FALSE) %>% 
  binarize()

# Warning message:
#   In count_ngrams(splitted_sequences, 3, c("A", "R", "N", "D", "C",  :
#                                              'seq' contains following unigrams not present in 'u' parameter:
#                                              NA

y <- factor(c(rep("banana", 50), rep("piolin", 50)))

test_ids <- unlist(lapply(levels(y), function(ith_group) {
  belonging_to_group <- y == ith_group
  sample(which(belonging_to_group), size = round(holdout_perc * sum(belonging_to_group), 0))
}))

train_kmer <- kmers[-test_ids, ]

train_label_ids <- y[-test_ids]


all_importances <- test_features(train_label_ids == "banana", train_kmer, criterion = "ig", adjust = NULL)

imp_features <- cut(all_importances, breaks = c(0, 0.05, 1))[[1]]

full_dat_small <- data.frame(as.matrix(kmers[, imp_features]), y = y)
full_dat_large <- data.frame(as.matrix(kmers), y = y)

model_rf_small <- ranger(y ~ ., data = full_dat_small[-test_ids, ], importance = "impurity", probability = TRUE)

model_rf_large <- ranger(y ~ ., data = full_dat_large[-test_ids, ], importance = "impurity", probability = TRUE)

predictions_df <- data.frame(small_model = predict(model_rf_small, full_dat_small[test_ids, ])[["predictions"]][, 2],
                             large_model = predict(model_rf_large, full_dat_large[test_ids, ])[["predictions"]][, 2])

res <- HMeasure(full_dat_small[test_ids, "y"], predictions_df, threshold = 0.5)[["metrics"]]
