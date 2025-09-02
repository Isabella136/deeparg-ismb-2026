amr_mapping_difference_files = [
    "deeparg_results/arg_alignment_identity_30/AMR_category_difference.tsv",
    "deeparg_results/arg_alignment_identity_50/AMR_category_difference.tsv",
    "deeparg_results/arg_alignment_identity_80/AMR_category_difference.tsv"
    ]

short_sequence_difference_summary_files = [
    "deeparg_results/arg_alignment_identity_30/difference_summary.tsv",
    "deeparg_results/arg_alignment_identity_50/difference_summary.tsv",
    "deeparg_results/arg_alignment_identity_80/difference_summary.tsv"
    ]

for index in range(3):
    mapping_frequencies = dict()
    footnote = ""
    for i, line in enumerate(open(amr_mapping_difference_files[index])):
        fields = line.strip().split('\t')
        if fields[0] == 'line#':
            continue
        elif len(fields) == 1:
            footnote = footnote + line
            continue
        best_arg = fields[1].split('|')[-1].upper()
        best_amr = fields[1].split('|')[-2]
        
        pred_arg = fields[2]
        pred_amr = fields[6]
    
        if (best_arg, best_amr, pred_arg, pred_amr) in mapping_frequencies:
            mapping_frequencies[(best_arg, best_amr, pred_arg, pred_amr)] += 1
        else:
            mapping_frequencies.update({(best_arg, best_amr, pred_arg, pred_amr):1})

    with open(short_sequence_difference_summary_files[index], 'w') as output:
        output.write("ARG best-hit\tAMR best-hit\tARG pred\tAMR pred\tFrequency\n")
        for keys, value in mapping_frequencies.items():
            output.write('\t'.join(keys) + "\t{}\n".format(value))
        output.write(footnote)