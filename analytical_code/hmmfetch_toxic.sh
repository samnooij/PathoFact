#! /bin/bash -l

hmmfetch="/work/projects/ecosystem_biology/local_tools/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmfetch"

#KEGG
db=KEGG
hmm_file=/work/projects/ecosystem_biology/data/hmm/KEGG/KO.hmm
IDs=ID_mini_KEGG_newest_2.txt
keggout_file=$db.hmm
$hmmfetch -f $hmm_file $IDs >>$keggout_file

#META-cyc
#db=metacyc
#hmm_file=/work/projects/ecosystem_biology/data/hmm/$db/$db.hmm
#metacycout_file=$db.hmm
#IDs=ID_mini_MetaCyc.txt
#$hmmfetch -f $hmm_file $IDs >>$metacycout_file

#Pfam
db=Pfam-A
hmm_file=/work/projects/ecosystem_biology/data/hmm/$db/$db.hmm
pfamout_file=$db.hmm
IDs=ID_mini_Pfam_newest_2.txt
$hmmfetch -f $hmm_file $IDs >>$pfamout_file

#TIGR
db=tigrpfam
hmm_file=/work/projects/ecosystem_biology/data/hmm/TIGRPFAM/$db.hmm
tigrpfamout_file=$db.hmm
IDs=ID_mini_TIGR_newest_2.txt
$hmmfetch -f $hmm_file $IDs >>$tigrpfamout_file

#swissprot
db=swissprot
hmm_file=/work/projects/ecosystem_biology/data/hmm/SwissProt/$db.hmm
swissprotout_file=$db.hmm
IDs=ID_mini_Swiss_newest_2.txt
$hmmfetch -f $hmm_file $IDs >>$swissprotout_file
