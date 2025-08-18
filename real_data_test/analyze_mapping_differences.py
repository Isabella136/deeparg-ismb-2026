from Bio.SeqIO import FastaIO

biosamples = [line.strip() for line in open("real_samples.txt", "r")]

arg_specific_cdd = dict()
arg_superfamily_cdd = dict()
for part in range(1, 26):
    for line in open("../CDD_features/Part{}_hitdata.txt".format(part), "r"):
        fields = line.strip().split('\t')
        if fields[0] == "Query":
            continue
        fields[0] = fields[0].split('>')[1]
        if (fields[0]) in arg_specific_cdd:
            arg_specific_cdd[fields[0]].append(fields[7])
            if fields[1] == "specific":
                arg_superfamily_cdd[fields[0]].append(fields[10])
            else:
                arg_superfamily_cdd[fields[0]].append(fields[7])
        else:
            arg_specific_cdd.update({fields[0]:[fields[7]]})
            if fields[1] == "specific":
                arg_superfamily_cdd.update({fields[0]:[fields[10]]})
            else:
                arg_superfamily_cdd.update({fields[0]:[fields[7]]})

specific_cdd_amr_ref_count = dict()
superfamily_cdd_amr_ref_count = dict()
arg_amr_ref_count = dict()
arg_list = list()
with open("../data/database/v2/features.fasta") as features:
    for values in FastaIO.SimpleFastaParser(features):
        arg = values[0].split('|')[-1].upper()
        amr = values[0].split('|')[-2]
        if values[0] in arg_specific_cdd:
            cdd_list = arg_specific_cdd[values[0]]
            cdd_sup_list = arg_superfamily_cdd[values[0]]
            for cdd in cdd_list:
                if (cdd, amr) in specific_cdd_amr_ref_count:
                    specific_cdd_amr_ref_count[(cdd, amr)] += 1
                else:
                    specific_cdd_amr_ref_count.update({(cdd, amr):1})
                if ((arg, cdd), amr) in arg_amr_ref_count:
                    arg_amr_ref_count[((arg, cdd), amr)] += 1
                else:
                    arg_amr_ref_count.update({((arg, cdd), amr):1})
            for cdd_sup in cdd_sup_list:
                if (cdd_sup, amr) in superfamily_cdd_amr_ref_count:
                    superfamily_cdd_amr_ref_count[(cdd_sup, amr)] += 1
                else:
                    superfamily_cdd_amr_ref_count.update({(cdd_sup, amr):1})
            if not (cdd_sup_list, cdd_list, arg) in arg_list:
                arg_list.append((cdd_sup_list, cdd_list, arg))
            
        else:
            if (arg, amr) in specific_cdd_amr_ref_count:
                specific_cdd_amr_ref_count[(arg, amr)] += 1
                superfamily_cdd_amr_ref_count[(arg, amr)] += 1
            else:
                specific_cdd_amr_ref_count.update({(arg, amr):1})
                superfamily_cdd_amr_ref_count.update({(arg, amr):1})
            if not ([arg], [arg], arg) in arg_list:
                arg_list.append(([arg], [arg], arg))
            if ((arg, arg), amr) in arg_amr_ref_count:
                arg_amr_ref_count[((arg, arg), amr)] += 1
            else:
                arg_amr_ref_count.update({((arg, arg), amr):1})

for prefix in ["spades/deeparg_results", "deeparg_results"]:
    for perc in [30, 50, 80]:
        output_file = "_".join(prefix.split('/')) + "_{}.tsv".format(perc)
        arg_map_diff_freq = dict()
        specific_map_diff_freq = dict()
        superfamily_map_diff_freq = dict()
        arg_map_diff = list()
        for bio in biosamples:
            input = "{}/{}/arg_alignment_identity_{}/AMR_category_difference.tsv".format(
                bio, prefix, perc)
            for line in open(input, "r"):
                fields = line.strip().split('\t')
                if fields[0] == "line#" or len(fields) == 1:
                    continue
                best_cdd_list = arg_specific_cdd[fields[1]]
                best_cdd_sup_list = arg_superfamily_cdd[fields[1]]
                best_arg = fields[1].split('|')[-1].upper()
                best_index = arg_list.index((best_cdd_sup_list, best_cdd_list, best_arg))
                best_amr = fields[1].split('|')[-2]
                pred_cdd_list = arg_specific_cdd[fields[7]]
                pred_cdd_sup_list = arg_superfamily_cdd[fields[7]]
                pred_arg = fields[2]
                pred_index = arg_list.index((pred_cdd_sup_list, pred_cdd_list, pred_arg))
                pred_amr = fields[6]
                for best_cdd, best_cdd_sup in zip(best_cdd_list, best_cdd_sup_list):
                    for pred_cdd, pred_cdd_sup in zip(pred_cdd_list, pred_cdd_sup_list):
                        if ((best_index, pred_index), (best_cdd_sup, best_cdd, best_arg, best_amr), 
                            (pred_cdd_sup, pred_cdd, pred_arg, pred_amr)) not in arg_map_diff:
                            arg_map_diff.append(
                                ((best_index, pred_index), (best_cdd_sup, best_cdd, best_arg, best_amr), 
                                (pred_cdd_sup, pred_cdd, pred_arg, pred_amr)))
                        if ((best_cdd_sup, best_cdd, best_arg, best_amr), 
                            (pred_cdd_sup, pred_cdd, pred_arg, pred_amr)) in arg_map_diff_freq:
                            arg_map_diff_freq[
                                ((best_cdd_sup, best_cdd, best_arg, best_amr), 
                                (pred_cdd_sup, pred_cdd, pred_arg, pred_amr))] += 1
                        else:
                            arg_map_diff_freq.update({
                                ((best_cdd_sup, best_cdd, best_arg, best_amr), 
                                (pred_cdd_sup, pred_cdd, pred_arg, pred_amr)):1})
                        if ((best_cdd_sup, best_cdd, best_amr), 
                            (pred_cdd_sup, pred_cdd, pred_amr)) in specific_map_diff_freq:
                            specific_map_diff_freq[
                                ((best_cdd_sup, best_cdd, best_amr), 
                                (pred_cdd_sup, pred_cdd, pred_amr))] += 1
                        else:
                            specific_map_diff_freq.update({
                                ((best_cdd_sup, best_cdd, best_amr), 
                                (pred_cdd_sup, pred_cdd, pred_amr)):1})
                        if ((best_cdd_sup, best_amr), 
                            (pred_cdd_sup, pred_amr)) in superfamily_map_diff_freq:
                            superfamily_map_diff_freq[
                                ((best_cdd_sup, best_amr), 
                                (pred_cdd_sup, pred_amr))] += 1
                        else:
                            superfamily_map_diff_freq.update({
                                ((best_cdd_sup, best_amr), 
                                (pred_cdd_sup, pred_amr)):1})

        with open(output_file, "w") as output:
            output.write("Superfamily best-hit\tCDD best-hit\tARG best-hit\tAMR best-hit\tARG best-hit index\t")
            output.write("Superfamily-AMR count best-hit\tCDD-AMR count best-hit\tARG-AMR count best-hit\t")
            output.write("Superfamily pred\tCDD pred\tARG pred\tAMR pred\tARG pred index\t")
            output.write("Superfamily-AMR count pred\tCDD-AMR count pred\tARG-AMR count pred\t")
            output.write("Superfamily Mapping Frequency\tCDD Mapping Frequency\tARG Mapping Frequency\t")
            output.write("Superfamily Vice-Versa Frequency\tCDD Vice-Versa Frequency\tARG Vice-Versa Frequency\n")

            for key in sorted(arg_map_diff):
                
                best_hit_super_count = superfamily_cdd_amr_ref_count[(key[1][0], key[1][3])]
                best_hit_cdd_count = specific_cdd_amr_ref_count[(key[1][1], key[1][3])]
                best_hit_arg_count = arg_amr_ref_count[((key[1][2], key[1][1]), key[1][3])]
                pred_super_count = superfamily_cdd_amr_ref_count[(key[2][0], key[2][3])]
                pred_cdd_count = specific_cdd_amr_ref_count[(key[2][1], key[2][3])]
                pred_arg_count = arg_amr_ref_count[((key[2][2], key[2][1]), key[2][3])]
                
                best_index = key[0][0]
                pred_index = key[0][1]

                arg_value = arg_map_diff_freq[(key[1], key[2])]
                if (key[2], key[1]) in arg_map_diff_freq:
                    arg_vice_versa = arg_map_diff_freq[(key[2], key[1])]
                else:
                    arg_vice_versa = 0
                if arg_vice_versa > arg_value:
                    continue

                cdd_key = (
                    (key[1][0], key[1][1], key[1][3]),
                    (key[2][0], key[2][1], key[2][3]))
                cdd_value = specific_map_diff_freq[cdd_key]
                if (cdd_key[1], cdd_key[0]) in specific_map_diff_freq:
                    cdd_vice_versa = specific_map_diff_freq[(cdd_key[1], cdd_key[0])]
                else:
                    cdd_vice_versa = 0

                super_key = (
                    (key[1][0], key[1][3]),
                    (key[2][0], key[2][3]))
                super_value = superfamily_map_diff_freq[super_key]
                if (super_key[1], super_key[0]) in superfamily_map_diff_freq:
                    super_vice_versa = superfamily_map_diff_freq[(super_key[1], super_key[0])]
                else:
                    super_vice_versa = 0

                output.write('\t'.join(key[1]) + '\t{}\t'.format(best_index))
                output.write("{}\t{}\t{}\t".format(
                    best_hit_super_count, best_hit_cdd_count, best_hit_arg_count))
                output.write('\t'.join(key[2]) + '\t{}\t'.format(pred_index))
                output.write("{}\t{}\t{}\t".format(
                    pred_super_count, pred_cdd_count, pred_arg_count))
                output.write(str(super_value) + '\t' + str(cdd_value) + '\t' + str(arg_value) + '\t')
                output.write(str(super_vice_versa) + '\t' + str(cdd_vice_versa) + '\t' + str(arg_vice_versa) + '\n')