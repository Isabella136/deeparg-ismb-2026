#!/bin/bash
#SBATCH --job-name=fetch_info
#SBATCH --output=fetch_info.out
#SBATCH --error=fetch_info.err
#SBATCH --time=12:00:00
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --qos=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g

grep "^N" ../samples/SAMN*/spades/blast_results.txt| cut -f2  | sort | uniq -d > protein_ids.txt
mapfile -t protein_ids < protein_ids.txt

echo -e "Protein ID\tProtein Length\tGene Name|Protein Name\tCDD start\tCDD end\tCDD ID\tSuperfamily" > protein_dict.tsv

for protein_id in "${protein_ids[@]}"; do
    #Get Protein ID + tab
    line="$protein_id\t"
    #Save efetch output for repeated use; because it is multiline, must run echo "${efetch_res}" as opposed to echo ${efetch_res}
    efetch_res=`efetch -db protein -id $protein_id -format ft`
    #Add protein end position + tab
    line+=`echo "${efetch_res}" | grep "Protein$" | awk -F'\t' '{print $2"\t"}'`
    #Add gene name 
    line+=`echo "${efetch_res}" | grep "^[[:space:]]*gene" | cut -f5 | tr '\n' '|'`
    #Add protein name + tab
    line+=`echo "${efetch_res}" | grep "^[[:space:]]*product" | cut -f5 | tr '\n' '\t'`
    #Get region location
    cdd_loc=`echo "${efetch_res}" | grep "Region$" | awk -F'\t' '{print $1"\n"$2}'`
    #get region ID
    cdd_id=`echo "${efetch_res}" | grep -A2 "^[[:space:]]*region"| grep "^[[:space:]]*db_xref" | cut -f5 | awk -F':' '{print $NF}'`
    i=0
    cdd_id_array=(${cdd_id})
    cdd_loc_array=(${cdd_loc})
    for id in "${cdd_id_array[@]}"; do
        db_id=`esearch -db cdd -query $id | efetch -format docsum | grep '<Accession>' | tr '>' '\t'| tr '<' '\t' | cut -f3`
        super=`curl -s "https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=$id" | grep -o 'cl[0-9]\{5,\}' | head -n 1`
        #Add region location and ID
        j=$((i + 1))
        cdd_start="${cdd_loc_array[$i]}"
        echo "$cdd_start"
        cdd_end="${cdd_loc_array[$j]}"
        line+="${cdd_start}\t${cdd_end}\t${db_id}\t${super}\t"
        i=$((i + 2))
    done
    echo -e "${line}" >> protein_dict.tsv
done