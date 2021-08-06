#!/usr/bin/env R

# logging
sink(file=file(snakemake@log[[1]], open="wt"), type="message")

library(tidyverse)
# RGI
RGI <- read.delim(file = snakemake@input[["RGI"]])
RGI <- RGI %>% filter(ORF_ID != "ORF_ID")
RGI <- RGI %>% select(1,9,13)
colnames(RGI) <- c("ORF_ID","Best_Hit_ARO","ARG_SNPs")

# DeepARG
DeepARG <- read.delim(file = snakemake@input[["DeepARG"]])
DeepARG <- DeepARG %>% filter(X.ARG != "#ARG")
DeepARG <- DeepARG %>%select(1,4)
colnames(DeepARG)<- c("DeepARG_ARG","ORF_ID")

# AMR index
PathoFact_AMR_index <- read.delim(file = snakemake@input[["AMR_index"]])
PathoFact_AMR_index$ARG <- toupper(PathoFact_AMR_index$ARG)
PathoFact_AMR_index <- unique(PathoFact_AMR_index)

# Combine
Combined <- full_join(DeepARG, RGI, by="ORF_ID")
Combined$DeepARG_ARG <- fct_explicit_na(as.character(Combined$DeepARG_ARG), na_level = "-")
Combined$Best_Hit_ARO <- fct_explicit_na(as.character(Combined$Best_Hit_ARO), na_level = "-")
Combined$ARG_SNPs <- fct_explicit_na(as.character(Combined$ARG_SNPs), na_level = "n/a")
Combined$Best_Hit_ARO <- toupper(Combined$Best_Hit_ARO)

Unique <- Combined %>% filter(Combined$DeepARG_ARG == "-" | Combined$Best_Hit_ARO == "-")
Unique$compare <- factor(rep("-", nrow(Unique)))
Combined <- Combined %>% filter(Combined$DeepARG_ARG != "-" & Combined$Best_Hit_ARO != "-")
Combined$compare <- as.character(Combined$DeepARG_ARG) == as.character(Combined$Best_Hit_ARO)
Combined <- rbind(Combined, Unique)
Combined$Prediction <- ifelse(Combined$compare == "TRUE", as.character(Combined$DeepARG_ARG),
                              ifelse(Combined$compare == "FALSE", as.character(Combined$Best_Hit_ARO), 
                                     ifelse(Combined$DeepARG_ARG == "-", as.character(Combined$Best_Hit_ARO),
                                            ifelse(Combined$Best_Hit_ARO == "-", as.character(Combined$DeepARG_ARG), "-"))))
Combined$Database <- ifelse(Combined$DeepARG_ARG == "-", "RGI",
                            ifelse(Combined$Best_Hit_ARO == "-", "DeepARG","DeepARG/RGI"))

Combined <- merge(Combined, PathoFact_AMR_index, by.x="Prediction", by.y = "ARG_Name", all.x = TRUE)

Combined <- Combined %>% select(3,8,5,9,10,11,7)

write.table(Combined, file = snakemake@output[["AMR_combined"]], sep="\t", row.names=FALSE, quote=FALSE)
