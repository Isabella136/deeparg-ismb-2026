#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate cd-hit

db1=../../data/database/v1_features.fasta
db2=../../data/database/v2_features.fasta

cd-hit -i $db1 -o db_v1_40 -T 0 -M 0 -d 0 -c 0.4 -n 2
cd-hit -i $db2 -o db_v2_40 -T 0 -M 0 -d 0 -c 0.4 -n 2
cd-hit -i db_rep_40.fasta -o db_rep_40 -T 0 -M 0 -d 0 -c 0.4 -n 2
clstr_renumber.pl db_all_40.clstr > db_all_renumbered_40.clstr
clstr_rev.pl db_all_renumbered_40.clstr db_rep_40.clstr > db_all_final_40.clstr