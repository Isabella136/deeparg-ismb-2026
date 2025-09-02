#!/bin/bash

#SBATCH --job-name=preprocess
#SBATCH --output=pre_log.out
#SBATCH --error=pre_log.err
#SBATCH --time=1:00:00
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --qos=highmem
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128g

module load bbmap
module load fastqc
mapfile -t biosamples < ../real_samples.txt

for biosample in "${biosamples[@]}"; do
    cd $biosample
    fastq1=$(ls ./*_1.fastq)
    fastq2=$(ls ./*_2.fastq)
    echo $fastq1
    srr=$(echo $fastq1 | awk -F "_" '{print $1}')
    echo $srr
    srun --exclusive --nodes=1 --mem=128g --output=./pre_log.out --error=./pre_log.err bash -c \
        """bbduk.sh in1=$fastq1 in2=$fastq2 out1=${srr}_1.quality_filtered.fastq out2=${srr}_2.quality_filtered.fastq 
            rm $fastq1 & rm $fastq2
            bbduk.sh in1=${srr}_1.quality_filtered.fastq in2=${srr}_2.quality_filtered.fastq \
                out1=${srr}_1.quality_trimmed.fastq out2=${srr}_2.quality_trimmed.fastq \
                qtrim=rl trimq=10 minlen=30 
            rm *.quality_filtered.fastq
            fastqc ${srr}_1.quality_trimmed.fastq ${srr}_2.quality_trimmed.fastq
            reformat.sh in1=${srr}_1.quality_trimmed.fastq in2=${srr}_2.quality_trimmed.fastq \
                out1=${srr}_1.quality_trimmed.fasta out2=${srr}_2.quality_trimmed.fasta
            rm *quality_trimmed.fastq
            reformat.sh in1=${srr}_1.quality_trimmed.fasta in2=${srr}_2.quality_trimmed.fasta \
                out=${srr}_interleaved.quality_trimmed.fasta
        """&
    cd ..
done
echo "done with loop"
wait