#!/usr/bin/env R

# logging
sink(file=file(snakemake@log[[1]], open="wt"), type="message")

library(tidyverse)
library(reshape2)

HMM_results<-read.delim(file = snakemake@input[["hmmresults"]], header = FALSE)
ID <- read.delim(file = snakemake@input[["ID"]], header = FALSE)
HMM_results<- merge(HMM_results, ID, by.x = "V1", by.y = "V1", all.y = TRUE)

positive <- read.table(file = snakemake@input[["positive"]], quote="\"", comment.char="")
negative <- read.table(file = snakemake@input[["negative"]], quote="\"", comment.char="")
shared <- read.table(file = snakemake@input[["shared"]], quote="\"", comment.char="")
positive$class <- "positive"
negative$class <- "negative"
shared$class <- "shared"

class_list<- list(positive, negative, shared)
Class <- Reduce(rbind, class_list)

HMM_results$V2.x <-gsub("_.*","",HMM_results$V2.x)
Virulence_hmm<-merge(HMM_results, Class, by.x = "V2.x", by.y = "V1", all.x = TRUE)
Virulence_hmm <- Virulence_hmm %>% select(2,5,6)
Virulence_hmm[is.na(Virulence_hmm)]<-'negative'
Virulence_hmm<-unique(Virulence_hmm)
Virulence_hmm$Prediction <- "True"
colnames(Virulence_hmm)<-c("ID","Sequence","class","Prediction")
Virulence_hmm <- dcast(Virulence_hmm, formula = ...~ class, value.var = "Prediction")

fncols <- function(Virulence_hmm, cname) {
  add <- cname[!cname%in%names(Virulence_hmm)]
  if(length(add)!=0) Virulence_hmm[add] <- NA
  Virulence_hmm
}

Virulence_hmm<-fncols(Virulence_hmm, c("positive", "negative", "shared"))
Virulence_hmm <- Virulence_hmm[, c("ID", "negative", "positive", "shared")]
Virulence_hmm[is.na(Virulence_hmm)]<-"False"
Virulence_hmm$ID<-sprintf("%010d",as.numeric(Virulence_hmm$"ID"))

write.table(Virulence_hmm, file = snakemake@output[[1]], sep="\t", row.names=F, quote=F)

