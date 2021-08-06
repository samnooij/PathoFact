#!/usr/bin/env R

# logging
sink(file=file(snakemake@log[[1]], open="wt"), type="message")

library(tidyverse)

# Input
## Toxin prediction
gene_table_HMM<-read.delim(file = snakemake@input[["input_HMM"]], comment.char="#")
gene_table_HMM <- gene_table_HMM[,c(1,2,4,3)]
colnames(gene_table_HMM) <- c("ORF_ID","HMM_Name","Score","Significance_evalue")
gene_table_HMM$ORF_ID<-sprintf("%010d",as.numeric(gene_table_HMM$ORF_ID))

## ORF translation
ORF_translation<-read.delim(file = snakemake@input[["translation"]], header=FALSE)
colnames(ORF_translation)<-c("ORF_ID","ORF")
ORF_translation$ORF<- sub('>', '', ORF_translation$ORF)
ORF_translation<-ORF_translation[,c(2,1)]
ORF_translation$ORF_ID<-sprintf("%010d",as.numeric(ORF_translation$ORF_ID))

## SignalP prediction
signalP<-read.table(file = snakemake@input[["signalP"]], sep= "\t", header=TRUE)
colnames(signalP)<- c("ORF_ID","Signal_peptide","Signal_prediction","Organism")
signalP$ORF_ID<-sprintf("%010d",as.numeric(signalP$ORF_ID))

## Gene library
library_HMM<-read.table(file = snakemake@input[["library"]], sep = ";", header = T, quote = '#')
colnames(library_HMM)<-c("HMM_ID","NAME","Alternative_name", "Database", "Description")

#Description file of the HMM
gene_table_library<-merge.data.frame(gene_table_HMM, library_HMM, by.x= "HMM_Name", by.y= "HMM_ID", all.x=TRUE)
gene_table_library<-gene_table_library[,c(2,1,3,4,5,6,7,8)]
gene_table_library<-merge.data.frame(ORF_translation, gene_table_library, by = "ORF_ID")

write.table(gene_table_library,file = snakemake@output[["gene_library"]], sep="\t", row.names=FALSE, quote=FALSE)

# Toxin prediction
gene_table_Toxic <- data.frame(table(gene_table_HMM[1]))
colnames(gene_table_Toxic) <- c("ORF_ID","Number_of_hits")
gene_table_Toxic$Toxin_prediction <- "pathogenic"
gene_table_Toxic <- left_join(ORF_translation, gene_table_Toxic, by="ORF_ID")
gene_table_Toxic$Toxin_prediction <- fct_explicit_na(gene_table_Toxic$Toxin_prediction, na_level = "non_pathogenic")
gene_table_Toxic[is.na(gene_table_Toxic)] <-0
gene_table_Toxic <- full_join(gene_table_Toxic, signalP, by="ORF_ID")

gene_table_Toxic$Toxin_confidence_level <- ifelse(gene_table_Toxic$Toxin_prediction == "pathogenic" & gene_table_Toxic$Signal_peptide == "Y", "1: Secreted Toxin",
                                                  ifelse(gene_table_Toxic$Toxin_prediction == "pathogenic" & gene_table_Toxic$Signal_peptide == "N", "2: Non-secreted Toxin", "-"))
gene_table_Toxic <- gene_table_Toxic %>% select(1,2,3,4,5,8)
write.table(gene_table_Toxic, file=snakemake@output[["gene_toxic"]], sep="\t", row.names=FALSE, quote=FALSE)
