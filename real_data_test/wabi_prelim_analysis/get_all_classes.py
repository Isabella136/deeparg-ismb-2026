from Bio.SeqIO import FastaIO

class_count = dict()

with open("../../data/database/v2/features.fasta") as features:
    for values in FastaIO.SimpleFastaParser(features):
        amr_class = values[0].split('|')[-2]
        if amr_class in class_count:
            class_count[amr_class] += 1
        else:
            class_count.update({amr_class:1})

print(class_count)