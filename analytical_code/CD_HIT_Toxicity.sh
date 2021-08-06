#! /bin/bash -l



cdhit="/work/projects/ecosystem_biology/local_tools/cd-hit-v4.6.1-2012-08-27_OpenMP/cd-hit-2d"
threads="1"
dbA=/work/users/ldenies/Toxin_HMM_database/ final_T3DB.faa
dbB=/work/users/ldenies/Toxin_HMM_database/T3DB_selected.faa

$cdhit -i $dbA -i2 $dbB -o clustered.DBs -c 0.9 -n 5 -G 0 -aS 0.70 -g 1 -s2 0.1 -S2 50000 -T $threads -M 10000

#output clustered.DBs contains all sequences in dbB that are not in dbA -> dbA keeps all sequences and the new file is non-overlapping

grep ">" clustered.DBs | wc -l
grep ">" $dbB | wc -l


