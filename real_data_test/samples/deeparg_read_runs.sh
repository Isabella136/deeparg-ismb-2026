#!/bin/bash

#SBATCH --job-name=deeparg_read
#SBATCH --output=deeparg_read_log.out
#SBATCH --error=deeparg_read_log.err
#SBATCH --time=12:00:00
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --qos=highmem
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128g

MAX_JOBS=10
running_jobs=0

eval "$(conda shell.bash hook)"
conda activate ../../deeparg/deeparg_env/

cd ../../
mapfile -t biosamples < real_data_test/real_samples.txt

for biosample in "${biosamples[@]}"; do
    fasta=$(ls ./real_data_test/sample/$biosample/*_interleaved.quality_trimmed.fasta)
    echo $fasta
    mkdir -p ./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_30
    mkdir -p ./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_50
    mkdir -p ./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_80
    if ! test -f ./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_30/X.mapping.ARG; then
        ((running_jobs++))
        echo "running_jobs = $running_jobs"
        srun --exclusive --nodes=1 --mem=128g --output=./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_30/log.out \
            --error=./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_30/log.err bash -c \
            "deeparg predict \
                --model SS \
                -i $fasta \
                -o ./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_30/X \
                -d data/ \
                --type nucl \
                --min-prob 0.8 \
                --arg-alignment-identity 30 \
                --arg-alignment-evalue 1e-10 \
                --arg-num-alignments-per-entry 1000"&
        if (( running_jobs >= MAX_JOBS )); then
            echo "waiting"
            wait -n
            ((running_jobs--))
        fi
    fi
    
    if ! test -f ./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_50/X.mapping.ARG; then
        ((running_jobs++))
        echo "running_jobs = $running_jobs"
        srun --exclusive --nodes=1 --mem=128g --output=./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_50/log.out \
            --error=./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_50/log.err bash -c \
            "deeparg predict \
                --model SS \
                -i $fasta \
                -o ./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_50/X \
                -d data/ \
                --type nucl \
                --min-prob 0.8 \
                --arg-alignment-identity 50 \
                --arg-alignment-evalue 1e-10 \
                --arg-num-alignments-per-entry 1000"&
        if (( running_jobs >= MAX_JOBS )); then
            echo "waiting"
            wait -n
            ((running_jobs--))
        fi
    fi
    if ! test -f ./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_80/X.mapping.ARG; then
        ((running_jobs++))
        echo "running_jobs = $running_jobs"
        srun --exclusive --nodes=1 --mem=128g --output=./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_80/log.out \
            --error=./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_80/log.err bash -c \
            "deeparg predict \
                --model SS \
                -i $fasta \
                -o ./real_data_test/sample/$biosample/deeparg_results/arg_alignment_identity_80/X \
                -d data/ \
                --type nucl \
                --min-prob 0.8 \
                --arg-alignment-identity 80 \
                --arg-alignment-evalue 1e-10 \
                --arg-num-alignments-per-entry 1000"&
        if (( running_jobs >= MAX_JOBS )); then
            echo "waiting"
            wait -n
            ((running_jobs--))
        fi
    fi
done
echo "done with loop"
wait