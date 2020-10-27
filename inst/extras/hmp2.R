library(microbiome)
library(tidyr)
set.seed(4235421)
#BiocManager::install("HMP2Data")
library(HMP2Data)
ps <- momspi16S()
df <- meta(ps) %>%  
  filter(sample_body_site=="rectum") %>% 
  group_by(subject_id) %>% tally(n()) %>% 
  arrange(n) %>% 
  filter(n>=10)

ps.sub <- subset_samples(ps, subject_id %in% df$subject_id)
ps.sub <- subset_samples(ps.sub, sample_body_site =="rectum")
hmp2 <- prune_taxa(taxa_sums(ps.sub) > 0, ps.sub)
sample_data(hmp2) <- sample_data(hmp2)[,c("size","sample_id","subject_id", 
                                        "visit_number", "sample_body_site")]

labels <- unique(meta(hmp2)$subject_id)

sample_data(hmp2)$subject_id <- gsub(labels[1], "Participant_1",sample_data(hmp2)$subject_id)
sample_data(hmp2)$subject_id <- gsub(labels[2], "Participant_2",sample_data(hmp2)$subject_id)
sample_data(hmp2)$subject_id <- gsub(labels[3], "Participant_3",sample_data(hmp2)$subject_id)
sample_data(hmp2)$subject_id <- gsub(labels[4], "Participant_4",sample_data(hmp2)$subject_id)
sample_data(hmp2)$subject_id <- gsub(labels[5], "Participant_5",sample_data(hmp2)$subject_id)
sample_data(hmp2)$subject_id <- gsub(labels[6], "Participant_6",sample_data(hmp2)$subject_id)
sample_data(hmp2)$subject_id <- gsub(labels[7], "Participant_7",sample_data(hmp2)$subject_id)
sample_data(hmp2)$subject_id <- gsub(labels[8], "Participant_8",sample_data(hmp2)$subject_id)
sample_data(hmp2)$subject_id <- gsub(labels[9], "Participant_9",sample_data(hmp2)$subject_id)
sample_data(hmp2)$subject_id <- gsub(labels[10], "Participant_10",sample_data(hmp2)$subject_id)
sample_data(hmp2)$subject_id <- gsub(labels[11], "Participant_11",sample_data(hmp2)$subject_id)
sample_data(hmp2)$subject_id <- gsub(labels[12], "Participant_12", sample_data(hmp2)$subject_id)
sample_data(hmp2)$subject_id <- gsub(labels[13], "Participant_13", sample_data(hmp2)$subject_id)
head(meta(hmp2))
save(hmp2, file = "hmp2.rda")
