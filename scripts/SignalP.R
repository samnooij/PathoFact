#!/usr/bin/env R

# logging
sink(file=file(snakemake@log[[1]], open="wt"), type="message")

library(tidyverse)

V5_p <- read.delim(file = snakemake@input[["SignalP_gramP"]], header=FALSE, comment.char="#")
V5_p$V5_p_SignalP <- ifelse(V5_p$V2 == "OTHER", "N","Y")
V5_p <- V5_p %>% select(1,2,8)
colnames(V5_p)<- c("ID","V5_p_Prediction","V5_p_SignalP")

V5_n <- read.delim(file = snakemake@input[["SignalP_gramN"]], header=FALSE, comment.char="#")
V5_n$V5_n_SignalP <- ifelse(V5_n$V2 == "OTHER", "N","Y")
V5_n <- V5_n %>% select(1,2,8)
colnames(V5_n)<- c("ID","V5_n_Prediction","V5_n_SignalP")

SignalP <- full_join(V5_p, V5_n)
SignalP$SignalPeptide <- ifelse(SignalP$V5_p_SignalP == "Y" | SignalP$V5_n_SignalP =="Y", "Y","N")
SignalP$Organism <- ifelse(SignalP$V5_p_SignalP == "Y" & SignalP$V5_n_SignalP == "Y", "gram+ or gram-",
                           ifelse(SignalP$V5_p_SignalP == "Y" & SignalP$V5_n_SignalP == "N", "gram+",
                                  ifelse(SignalP$V5_p_SignalP == "N" & SignalP$V5_n_SignalP == "Y", "gram-","-")))

SignalP$compare <- as.character(SignalP$V5_p_Prediction) == as.character(SignalP$V5_n_Prediction)
SignalP$SignalP_prediction <- ifelse(SignalP$Organism == "gram+ or gram-" & SignalP$compare == "FALSE",paste(SignalP$V5_p_Prediction, SignalP$V5_n_Prediction, sep = " / "),
                                     ifelse(SignalP$SignalPeptide == "Y" & SignalP$compare == "TRUE", as.character(SignalP$V5_p_Prediction), "-"))
SignalP <- SignalP %>% select(1,6,9,7)

write.table(SignalP, file = snakemake@output[["SignalP_report"]], sep="\t", row.names=FALSE, quote=FALSE)
