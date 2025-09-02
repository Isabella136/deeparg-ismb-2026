amr_mapping_difference_files = [
    "deeparg_results/arg_alignment_identity_30/AMR_category_difference.tsv",
    "deeparg_results/arg_alignment_identity_50/AMR_category_difference.tsv",
    "deeparg_results/arg_alignment_identity_80/AMR_category_difference.tsv"
    ]

long_sequence_difference_summary_files = [
    "deeparg_results/arg_alignment_identity_30/difference_summary.tsv",
    "deeparg_results/arg_alignment_identity_50/difference_summary.tsv",
    "deeparg_results/arg_alignment_identity_80/difference_summary.tsv"
    ]

orf_info = dict()

for line in open("orf.gbk"):
    stripped = line.strip()
    if "/note" in stripped:
        stripped = stripped[7:]
        fields = stripped.split(';')
        id = fields[0].split('=')[1]
        conf = fields[6].split('=')[1]
        orf_info.update({id:conf})

for index in range(3):
    with open(long_sequence_difference_summary_files[index], 'w') as output:
        output.write("ORF ID\tORF conf\tARG best-hit\tAMR best-hit\tARG pred\tAMR pred\n")
        for i, line in enumerate(open(amr_mapping_difference_files[index])):
            fields = line.strip().split('\t')
            if fields[0] == 'line#':
                continue
            elif len(fields) == 1:
                output.write(line)
                continue
            orf_full_name = fields[5].split('_')
            orf_id = orf_full_name[1] + '_' + orf_full_name[-1]
            orf_conf = orf_info[orf_id]

            best_arg = fields[1].split('|')[-1].upper()
            best_amr = fields[1].split('|')[-2]
            
            pred_arg = fields[2]
            pred_amr = fields[6]
        
            output.write('\t'.join((orf_id, orf_conf, best_arg, best_amr, pred_arg, pred_amr)) + '\n')