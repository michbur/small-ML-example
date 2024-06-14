library(dplyr)
library(ranger)
library(hmeasure)
library(kernlab)
library(biogram)

test_file_address <- "test_dataset_FullSeq.csv"
test_dataset <- read.csv(test_file_address, header = TRUE, sep = ",")

test_model <- test_dataset %>% 
  as.data.frame() %>% 
  select(-(Fold)) %>% 
  rename_with(~c("UniProt.Acc", "Full.seq", "PS")) %>% 
  mutate(actual_PS = ifelse(`PS` == "True", TRUE, FALSE)) %>% 
  select(-(PS))

train_file_address <- "train_dataset.csv"
train_dataset <- read.csv(train_file_address, header = TRUE, sep = ",")

#1. Removal of atypical AA with command Linux

#2. data preparation: Unigrams of the sequence in matrix with labeled (piolin/bananas as true/false)
#2.1 Split the sequences into character unigrams
training_seqs <- train_dataset %>% 
  select(Full.seq:actual_PS)

split_sequence <- function(seq) {
  unlist(strsplit(seq, split = ""))
}

# 2.1 List; its elements are char vectors (as in analysis-kmer)
splitted_sequences <- lapply(training_seqs$Full.seq, split_sequence) %>% 
  list2matrix()
# 
# # Create a new dataframe to store the final result
# result <- data.frame()
# 
# # Iterate over the splitted_sequences list and bind rows to the result dataframe
# for (i in 1:length(splitted_sequences)) {
#   # Create a temporary dataframe with the label and the splitted sequence
#   temp_df <- data.frame(label = training_seqs$actual_PS[i], t(as.data.frame(splitted_sequences[[i]], stringsAsFactors = FALSE)))
#   # Bind the temporary dataframe to the result
#   result <- bind_rows(result, temp_df)
# }
# 
# # Reset row names for a clean output
# rownames(result) <- NULL
# 
# 
# training_mtx <- result %>%
#   list2matrix()
# 
# transposed_training_seqs <- t(training_mtx  )
# 
# u_seq <- unique(as.vector(transposed_training_seqs ))

kmers_seq <- count_ngrams(splitted_sequences , 3, c("A", "R", "N", "D", "C", "Q", "E", "G", "H",
                                 "I", "L", "K", "M", "F", "P", "S", "T", "W", 
                                 "Y", "V"), pos = FALSE) %>% 
  binarize()



model_rf <- ranger(actual_PS ~ ., data = train_dataset[["Full.seq"]], importance = "impurity")

res <- HMeasure(test_model[["actual_PS"]], data.frame(ranger = predict(model_rf, train_dataset)[["predictions"]]))[["metrics"]]

sort(importance(model_rf), decreasing = TRUE)


#try grouping FRY / GSQN / others => sticker-spacer
detect <- train_dataset %>% 
  str_detect(Full.seq, 'U')
