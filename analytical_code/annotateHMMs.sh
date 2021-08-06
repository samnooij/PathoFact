#! /bin/bash -l

cpus=12 # run on with -l nodes=1
hmmscan_tool="/work/projects/ecosystem_biology/local_tools/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmsearch --cpu $cpus --noali --notextw "
fasta=swissprot.faa # change this to your filename

#KEGG
db=kegg
hmm_file=/work/projects/ecosystem_biology/data/hmm/KEGG/KO.hmm
keggout_file=$fasta.$db.hmmscan
time $hmmscan_tool --tblout $keggout_file $hmm_file $fasta >/dev/null
#out_file=$fasta
#perl /home/users/aheintzbuschart/myScripts/consolidate_hmmscan_results_justKEGG.pl $out_file $keggout_file

#META-cyc
db=metacyc
hmm_file=/work/projects/ecosystem_biology/data/hmm/$db/$db.hmm
metacycout_file=$fasta.$db.hmmscan
time $hmmscan_tool --tblout $metacycout_file $hmm_file $fasta >/dev/null

#Pfam
db=Pfam-A
hmm_file=/work/projects/ecosystem_biology/data/hmm/$db/$db.hmm
pfamout_file=$fasta.$db.hmmscan
time $hmmscan_tool --tblout $pfamout_file $hmm_file $fasta >/dev/null

#TIGR
db=tigrpfam
hmm_file=/work/projects/ecosystem_biology/data/hmm/TIGRPFAM/$db.hmm
tigrpfamout_file=$fasta.$db.hmmscan
time $hmmscan_tool --tblout $tigrpfamout_file $hmm_file $fasta >/dev/null

#swissprot
db=swissprot
hmm_file=/work/projects/ecosystem_biology/data/hmm/SwissProt/$db.hmm
swissprotout_file=$fasta.$db.hmmscan
time $hmmscan_tool --tblout $swissprotout_file $hmm_file $fasta >/dev/null

date

#consolidate
perl /work/projects/ecosystem_biology/local_tools/perlscripts/hmms/consolidate_hmmscan_results.pl $fasta $keggout_file $metacycout_file $pfamout_file $swissprotout_file $tigrpfamout_file

date

input=${fasta}pfam.tsv
idn=pfamID
python /home/users/aheintzbuschart/myScripts/150705_MUST_hmmParsePfam.py $input $idn -g $(grep ">" $fasta | wc -l)

input=${fasta}kegg.tsv
idn=KO
python /home/users/aheintzbuschart/myScripts/150705_MUST_hmmParse.py $input $idn -g $(grep ">" $fasta | wc -l)

input=${fasta}metacyc.tsv
idn=metaCycID
python /home/users/aheintzbuschart/myScripts/150705_MUST_hmmParse.py $input $idn -g $(grep ">" $fasta | wc -l)

input=${fasta}swissprot.tsv
idn=swissprotEC
python /home/users/aheintzbuschart/myScripts/150705_MUST_hmmParse.py $input $idn -g $(grep ">" $fasta | wc -l)

input=${fasta}tigrpfam.tsv
idn=tigrID
python /home/users/aheintzbuschart/myScripts/150705_MUST_hmmParse.py $input $idn -g $(grep ">" $fasta | wc -l)

python /home/users/aheintzbuschart/myScripts/150310_MUST_hmmBestAll.py ${fasta}kegg.tsv ${fasta}metacyc.tsv ${fasta}swissprot.tsv ${fasta}pfam.tsv ${fasta}tigrpfam.tsv -g $(grep ">" $fasta | wc -l)

date
