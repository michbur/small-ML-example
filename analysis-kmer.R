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

full_dat <- data.frame(y = factor(c(rep("banana", 50), rep("piolin", 50))),
                       as.matrix(kmers))

test_ids <- unlist(lapply(levels(df[["y"]]), function(ith_group) {
  belonging_to_group <- full_dat[["y"]] == ith_group
  sample(which(belonging_to_group), size = round(holdout_perc * sum(belonging_to_group), 0))
}))

train_dat <- full_dat[-test_ids, ]
test_dat <- full_dat[test_ids, ]

model_rf <- ranger(y ~ ., data = train_dat, importance = "impurity")

res <- HMeasure(test_dat[["y"]], data.frame(ranger = predict(model_rf, test_dat)[["predictions"]]))[["metrics"]]

sort(importance(model_rf), decreasing = TRUE)
