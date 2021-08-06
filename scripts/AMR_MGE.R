#!/usr/bin/env R

# logging
sink(file=file(snakemake@log[[1]], open="wt"), type="message")

library(tidyverse)

# ORF Translation
ORF_translation<-read.delim(file = snakemake@input[["ORF_translation"]], header=FALSE)
colnames(ORF_translation) <- c("ORF_ID", "ORF")
ORF_translation$ORF<- sub('>', '', ORF_translation$ORF)
ORF_translation$ORF_ID <- sprintf("%010d", ORF_translation$ORF_ID)

# Mapping File
Contig_ORFs <- Contig_ORFs <- read.delim(file = snakemake@input[["Contig_ORFs"]], header=FALSE)
colnames(Contig_ORFs) <- c("Contig","ORF")

# AMR prediction
AMR <- read.delim(file = snakemake@input[["AMR"]])
AMR$ORF_ID <- sprintf("%010d", AMR$ORF_ID)
AMR <- merge(AMR, ORF_translation, by="ORF_ID", all=TRUE)
AMR <- AMR %>% mutate_all(list(~replace(as.character(.), is.na(.), "-")))
AMR <- merge(AMR, Contig_ORFs, by="ORF", all=TRUE)
AMR <- AMR %>% select(9,1,2,3,4,5,6,7,8)

# Contig Translation
Contig_translation <- read.delim(file = snakemake@input[["Contig_translation"]], header=FALSE)
colnames(Contig_translation) <- c("Contig_ID","Contig")
Contig_translation$Contig<- sub('>', '', Contig_translation$Contig)
Contig_translation$Contig_ID <- sprintf("%010d", Contig_translation$Contig_ID)

# Plasmid Prediction
## PlasFlow
PlasFlow <- read.delim(file = snakemake@input[["PlasFlow"]], header=TRUE)
PlasFlow <- separate(data= PlasFlow, col=label, into=c("left", "right"), sep="\\.")
PlasFlow <- PlasFlow %>% select(1,2)
colnames(PlasFlow) <-c("Contig_ID","PlasFlow_prediction")
PlasFlow <- PlasFlow %>% filter(Contig_ID != "contig_name")
PlasFlow$Contig_ID <- sprintf("%010d", PlasFlow$Contig_ID)

## MOBsuite
#MOB_suite <- read.delim(file = snakemake@input[["MOB_suite"]], header=TRUE)
#MOB_suite <- MOB_suite %>% filter(sample_id != "sample_id")
#MOB_suite$mash_neighbor_distance <- as.numeric(as.character(MOB_suite$mash_neighbor_distance))
#MOB_suite$MOB_suite_prediction <- if_else(MOB_suite$mash_neighbor_distance <= 0.06 , "plasmid","-")
#MOB_suite <- MOB_suite %>% select(1,27) %>% filter(MOB_suite_prediction == "plasmid")
#colnames(MOB_suite)<- c("Contig_ID","MOB_suite_prediction")
#MOB_suite$Contig_ID <- sprintf("%010d", MOB_suite$Contig_ID)

## Combine PlasFlow and MOB_suite
#Plasmid <- merge(PlasFlow, MOB_suite, by="Contig_ID", all = TRUE)
Plasmid <- merge(PlasFlow, Contig_translation, by= "Contig_ID", all= TRUE)
Plasmid$PlasFlow_prediction <- fct_explicit_na(Plasmid$PlasFlow_prediction, na_level = "unclassified")
#Plasmid$MOB_suite_prediction <- fct_explicit_na(Plasmid$MOB_suite_prediction, na_level = "unclassified")
#Plasmid$Plasmid_prediction <- ifelse(Plasmid$PlasFlow_prediction == "plasmid" | Plasmid$MOB_suite_prediction == "plasmid", "plasmid", "-")
#Plasmid$Plasmid_database <- ifelse(Plasmid$Plasmid_prediction == "plasmid" & Plasmid$PlasFlow_prediction == "plasmid" & Plasmid$MOB_suite_prediction == "plasmid", "PlasFlow/MOB_suite",
#                                   ifelse(Plasmid$Plasmid_prediction == "plasmid" & Plasmid$PlasFlow_prediction == "plasmid" & Plasmid$MOB_suite_prediction != "plasmid", "PlasFlow",
#                                          ifelse(Plasmid$Plasmid_prediction == "plasmid" & Plasmid$PlasFlow_prediction != "plasmid" & Plasmid$MOB_suite_prediction == "plasmid", "MOB_suite", "-")))
#Plasmid$Chromosome_prediction <- ifelse(Plasmid$PlasFlow_prediction == "chromosome" & Plasmid$MOB_suite_prediction != "plasmid", "chromosome","-")
#Plasmid <- Plasmid %>% select(4,1,7,5,6)  

Plasmid$Plasmid_prediction <- ifelse(Plasmid$PlasFlow_prediction == "plasmid", "plasmid", "-")
Plasmid$Chromosome_prediction <- ifelse(Plasmid$PlasFlow_prediction == "chromosome", "chromosome","-")
Plasmid <- Plasmid %>% select(1,3,5,4)

# Phage
## VirFinder
VirFinder <- read.delim(file = snakemake@input[["DeepVirFinder"]])
VirFinder <- VirFinder %>% select(1,3,4)
colnames(VirFinder)<- c("Contig_ID","VirFinder_score","VirFinder_pvalue" )
VirFinder <- VirFinder %>% filter(Contig_ID != "name")
VirFinder$VirFinder_pvalue <- as.numeric(as.character(VirFinder$VirFinder_pvalue))
VirFinder$VirFinder_score <- as.numeric(as.character(VirFinder$VirFinder_score))
VirFinder$VirFinder_prediction <- ifelse(VirFinder$VirFinder_pvalue <=0.05 & VirFinder$VirFinder_score >=0.70, "phage", "-")
#VirFinder$Contig_ID <- sprintf("%010d", VirFinder$Contig_ID)

## VirSorter
vs.pred <- read.csv(file = snakemake@input[["VirSorter"]], quote="", header=FALSE)
vs.head <- setNames(data.frame(matrix(ncol=12, nrow=0)), c("Contig_id","Nb.genes.contigs","Fragment","Nb.genes","Category","Nb.phage.hallmark.genes","Phage.gene.enrichment.sig","Non-Caudovirales.phage.gene.enrichment.sig","Pfam.depletion.sig","Uncharacterized.enrichment.sig","Strand.switch.depletion.sig","Short.genes.enrichment.sig"))
colnames(vs.pred)<-colnames(vs.head)
colnames(vs.pred)[1]<- "vs.id"
vs.cats <- do.call(rbind,strsplit(x=as.character(vs.pred$vs.id[grep("category",vs.pred$vs.id)]),split=" - ",fixed=T))[,2]
vs.num <- grep("category",vs.pred$vs.id)
vs.pred$Category <- paste(c("",rep.int(vs.cats, c(vs.num[-1],nrow(vs.pred)) - vs.num)), vs.pred$Category)
vs.pred <- vs.pred[-grep("#",vs.pred$vs.id),]
vs.pred$node <- gsub(pattern="VIRSorter_",replacement="",x=vs.pred$vs.id)
vs.pred$node <- gsub(pattern="-circular",replacement="",x=vs.pred$node)
vs.pred$node <- gsub(pattern="cov_(\\d+)_",replacement="cov_\\1.",x=vs.pred$node,perl=F)

VirSorter <- vs.pred %>% select(5,13)
VirSorter$node <- gsub("(^|[^0-9])0+", "\\1", VirSorter$node, perl = TRUE)
VirSorter$Virsorter_prediction <- ifelse(VirSorter$Category == "Prophages 3", "-",
                                         ifelse(VirSorter$Category == "Complete phage contigs 3", "-",
                                                ifelse(VirSorter$Category == "-", "-", "phage")))
VirSorter <- VirSorter %>% select(2,3)
colnames(VirSorter)<- c("Contig_ID","VirSorter_prediction")
VirSorter$Contig_ID <- as.numeric(as.character(VirSorter$Contig_ID))
#VirSorter$Contig_ID <- sprintf("%010d", VirSorter$Contig_ID)

## Combine Phage prediction
Phage <- merge(VirFinder, VirSorter, by="Contig_ID", all = TRUE)
if(all(is.na(Phage$VirSorter_prediction))){
  Phage$VirSorter_prediction <- factor(rep("-", nrow(Phage)))
} else {
  Phage$VirSorter_prediction <- fct_explicit_na(Phage$VirSorter_prediction, na_level = "-")
}
Phage$Phage_prediction <- ifelse(Phage$VirFinder_prediction == "phage" | Phage$VirSorter_prediction == "phage", "phage", "-")
Phage$Phage_database <- ifelse(Phage$Phage_prediction == "phage" & Phage$VirFinder_prediction == "phage" & Phage$VirSorter_prediction == "phage", "VIRSorter/DeepVirFinder",
                               ifelse(Phage$Phage_prediction == "phage" & Phage$VirFinder_prediction == "phage" & Phage$VirSorter_prediction != "phage", "DeepVirFinder",
                                      ifelse(Phage$Phage_prediction == "phage" & Phage$VirFinder_prediction != "phage" & Phage$VirSorter_prediction == "phage", "VIRSorter","-")))
Phage <- Phage %>% select(1,6,7)
Phage$Contig_ID <- sprintf("%010d", Phage$Contig_ID)

# Combine MGEs
MGEs <- merge(Plasmid, Phage, by = "Contig_ID", all=TRUE)
MGEs$MGE_prediction <- ifelse(MGEs$Plasmid_prediction == "plasmid" & MGEs$Phage_prediction == "-" & MGEs$Chromosome_prediction == "-", "plasmid",
                              ifelse(MGEs$Plasmid_prediction == "-" & MGEs$Phage_prediction == "phage" & MGEs$Chromosome_prediction == "-", "phage",
                                     ifelse(MGEs$Plasmid_prediction == "-" & MGEs$Phage_prediction == "-" & MGEs$Chromosome_prediction == "chromosome", "chromosome",
                                            ifelse(MGEs$Plasmid_prediction == "plasmid" & MGEs$Phage_prediction == "-" & MGEs$Chromosome_prediction == "chromosome", "ambiguous (plasmid/chromosome)",
                                                  ifelse(MGEs$Plasmid_prediction == "-" & MGEs$Phage_prediction == "phage" & MGEs$Chromosome_prediction == "chromosome", "ambiguous (phage/chromosome)",
                                                         ifelse(MGEs$Plasmid_prediction == "plasmid" & MGEs$Phage_prediction == "phage" & MGEs$Chromosome_prediction == "-", "ambiguous (plasmid/phage)","unclassified"))))))


# Combine AMR MGEs
AMR_MGE <- merge(AMR, MGEs, by="Contig", all.x = TRUE)
AMR_MGE <- AMR_MGE %>% select(1,10,2,3,4,5,6,7,8,9,11,12,13,14,15)                                         
write.table(AMR_MGE, file = snakemake@output[["Report_1"]], sep="\t", row.names=FALSE, quote=FALSE)

AMR_MGE_final <- AMR_MGE %>% select(1:10,15)
write.table(AMR_MGE_final, file = snakemake@output[["Report_2"]], sep="\t", row.names=FALSE, quote=FALSE)
