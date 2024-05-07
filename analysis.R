library(dplyr)
library(ranger)
library(hmeasure)
library(kernlab)

set.seed(151390)

n_randoms <- 20
holdout_perc <- 0.2

df <- data.frame(y = factor(c(rep("banana", 150), rep("piolin", 150))))

full_dat <- runif(300 * n_randoms) %>% 
  matrix(ncol = n_randoms) %>% 
  data.frame() %>% 
  setNames(paste0("randomFeature", 1L:n_randoms)) %>% 
  cbind(df, .) %>% 
  data.frame(., perfectVariable = c(runif(150, 0, 0.6), runif(150, 0.4, 1)))

test_ids <- unlist(lapply(levels(df[["y"]]), function(ith_group) {
  belonging_to_group <- full_dat[["y"]] == ith_group
  sample(which(belonging_to_group), size = round(holdout_perc * sum(belonging_to_group), 0))
}))

train_dat <- full_dat[-test_ids, ]
test_dat <- full_dat[test_ids, ]

model_rf <- ranger(y ~ ., data = train_dat)
model_svm <- ksvm(y ~ ., data = train_dat)

res <- HMeasure(test_dat[["y"]], 
                data.frame(ranger = predict(model_rf, test_dat)[["predictions"]],
                           svm = predict(model_svm, test_dat)))[["metrics"]]

