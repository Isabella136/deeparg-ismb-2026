# Description
This repository includes all the figures and experiments for "Database Description" and "Examining the Interaction between Alignment and Neural Networks". All figures require Python 3 to be installed

# Getting the Database
DeepARG's database can be downloaded by running

    deeparg download_data -o ./data

# Recreating Figures 1 and 2
Install Biopython (>= v1.85), pandas (>= v2.3.2), seaborn (v0.13.2) and networkx (v3.6.1)

    pip install biopython
    pip install pandas
    pip install seaborn
    pip install networkx

Run the following:

    cd database/database clustering
    ./cluster.sh
    cd ../
    python concat_ver_to_seq.py
    python feature_data_extraction.py
    python database_distribution_visualization.py 2

Alternatively, you can create figures for DeepARG version 1 by also running

    python database_distribution_visualization.py 1

# Running the WGS experiments (Slurm is recommented)
Install NCBI's Entrez Direct (installation link [here](https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Getting_Started)).  
Download and install [bbtools](https://github.com/bbushnell/BBTools?tab=readme-ov-file#-installation), [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc), [spades](https://github.com/ablab/spades/releases/tag/v4.2.0), and [prodigal](https://github.com/hyattpd/Prodigal/wiki/installation).  
Create a virtual conda environment to run DeepARG:
- Create a virtual environment with conda:

    ```conda create -n deeparg_env python=2.7.18
    source activate deeparg_env
- Install diamond with conda (inside virtual environment): 

    ```conda install -c bioconda diamond==0.9.24
- Optional (used for short reads pipeline): 

    ```conda install -c bioconda trimmomatic
    conda install -c bioconda vsearch
    conda install -c bioconda bedtools==2.29.2
    conda install -c bioconda bowtie2==2.3.5.1
    conda install -c bioconda samtools
- Activate virtual environment

    ```conda activate deeparg_env
Fetch the raw reads by running the `get_reads.sh` script:

    cd WGS_experiments/samples
    ./get_reads.sh

Once done, start a preprocessing job with sbatch by running

    sbatch ./preprocess_reads.sh

After the job is completed, you can now run DeepARG-SS:

    sbatch ./deeparg_read_runs.sh

Meanwhile, you can also assembled the reads into contigs with

    sbatch ./reads_assemble.sh

Once the assembly is done, you can now find ORFs and run DeepARG-LS:

    sbatch ./deeparg_contig_runs.sh




