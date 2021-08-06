#!/usr/bin/env R

# logging
sink(file=file(snakemake@log[[1]], open="wt"), type="message")

library(tidyverse)

# load virulence prediction
Virulence_factors <- read.delim(file=snakemake@input[["Virulence_factor"]])
Virulence_factors <- Virulence_factors %>% select(2,5,7)
Virulence_factors$ORF_ID <- sprintf("%010d", Virulence_factors$ORF_ID)

# load toxin prediction
Toxins <- read.delim(file=snakemake@input[["Toxins"]])
Toxins <- Toxins %>% select(2,4,5,6)
Toxins$ORF_ID <- sprintf("%010d", Toxins$ORF_ID)

# load AMR prediction
AMR_MGE <- read.delim(file=snakemake@input[["AMR_MGE"]])
AMR_MGE$ORF_ID <- sprintf("%010d", AMR_MGE$ORF_ID)
AMR_MGE$Contig_ID <- sprintf("%010d", AMR_MGE$Contig_ID)
AMR_MGE <- AMR_MGE %>% select(1:9,11)

# Combine files
Predictions_dfs<-list(Toxins, Virulence_factors, AMR_MGE)
PathoFact_predictions<-Reduce(full_join, Predictions_dfs)
PathoFact_predictions <- PathoFact_predictions[,c(1,9,8,7,2,4,3,5,6,10:15)]

write.table(PathoFact_predictions, file = snakemake@output[["PathoFact_report"]], sep="\t", row.names=FALSE, quote=FALSE)

