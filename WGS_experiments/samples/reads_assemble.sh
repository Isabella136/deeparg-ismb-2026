#!/bin/bash

#SBATCH --job-name=spades
#SBATCH --output=spades_log.out
#SBATCH --error=spades_log.err
#SBATCH --time=4:00:00
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --qos=highmem
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128g

module load spades
mapfile -t biosamples < ../real_samples.txt

for biosample in "${biosamples[@]}"; do
    cd $biosample
    fasta1=$(ls ./*_1.quality_trimmed.fasta)
    fasta2=$(ls ./*_2.quality_trimmed.fasta)
    rm -r spades
    mkdir spades
    srun --exclusive --nodes=1 --mem=128g --output=./spades_log.out --error=./spades_log.err bash -c \
        "spades.py -1 $fasta1 -2 $fasta2 -o spades --isolate -t 8 -m 128"&
    cd ..
done
echo "done with loop"
wait