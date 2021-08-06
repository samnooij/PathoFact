#!/usr/bin/env R

# logging
sink(file=file(snakemake@log[[1]], open="wt"), type="message")

library(tidyverse)

# HMM
HMM <- read.delim(file= snakemake@input[["hmm"]], header=FALSE)
colnames(HMM) <- c("ORF_ID","Virulence_HMM_prediction")

# classifier
classifier <- read.delim(file = snakemake@input[["classifier"]], header=FALSE)
colnames(classifier) <- c("ORF_ID","Virulence_classifier_prediction")

# ORF translation
ORF_translation<-read.delim(file = snakemake@input[["translation"]], header=FALSE)
colnames(ORF_translation)<-c("ORF_ID","ORF")
ORF_translation$ORF<- sub('>', '', ORF_translation$ORF)
ORF_translation<-ORF_translation[,c(2,1)]

# Complete Virulence prediction
Virulence <- full_join(HMM, classifier, by="ORF_ID")
Virulence <- full_join(ORF_translation, Virulence, by="ORF_ID")
Virulence$Virulence_prediction <- ifelse(Virulence$Virulence_HMM_prediction == "pathogenic" & Virulence$Virulence_classifier_prediction == "pathogenic", "pathogenic",
                                         ifelse(Virulence$Virulence_HMM_prediction == "unclassified" & Virulence$Virulence_classifier_prediction == "pathogenic", "pathogenic",
                                                ifelse(Virulence$Virulence_HMM_prediction == "negative" & Virulence$Virulence_classifier_prediction == "pathogenic", "unclassified",
                                                       ifelse(Virulence$Virulence_HMM_prediction == "pathogenic" & Virulence$Virulence_classifier_prediction == "negative", "unclassified","non_pathogenic"))))

# Signal P
SignalP <- read.delim(file = snakemake@input[["SignalP"]])
SignalP <- SignalP %>% select(1,2)
colnames(SignalP)<- c("ORF_ID","Signal_peptide")

# Combine Virulence prediction with SignalP prediction
Virulence <- full_join(Virulence, SignalP, by="ORF_ID")
Virulence$Virulence_confidence_level <- ifelse(Virulence$Virulence_prediction == "pathogenic" & Virulence$Signal_peptide == "Y", "1: Secreted Virulence factor",
                                               ifelse(Virulence$Virulence_prediction == "pathogenic"&  Virulence$Signal_peptide == "N", "2: Non-secreted Virulence factor",
                                                      ifelse(Virulence$Virulence_prediction == "unclassified" & Virulence$Signal_peptide == "Y", "3: Potential Secreted Virulence factor",
                                                             ifelse(Virulence$Virulence_prediction == "unclassified" & Virulence$Signal_peptide == "N", "4: Potential Non-secreted Virulence factor","-"))))
Virulence$ORF_ID <- sprintf("%010d",as.numeric(Virulence$ORF_ID))

write.table(Virulence, file=snakemake@output[["Virulence_report"]], sep="\t", row.names=FALSE, quote=FALSE)
