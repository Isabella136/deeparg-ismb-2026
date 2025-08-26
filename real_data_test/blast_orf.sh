#!/bin/bash
#SBATCH --job-name=orf_blast
#SBATCH --output=orf_blast_log.out
#SBATCH --error=orf_blast_log.err
#SBATCH --time=12:00:00
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --qos=highmem
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128g

eval "$(conda shell.bash hook)"
conda activate blast
mapfile -t biosamples < real_samples.txt

for biosample in "${biosamples[@]}"; do
    echo $biosample
    cd $biosample/spades
    srun --exclusive --nodes=1 --mem=128g --output=blast_log.out --error=blast_log.err bash -c \
        """ blastx -task blastx-fast -query orf.fasta -db ../../bacteria_refseq_protein_dedup -out blast_results.txt \
         -qcov_hsp_perc 75 -max_target_seqs 50 -outfmt 6 -num_threads 8 -evalue 0.01"""&
    cd ../..
done
wait