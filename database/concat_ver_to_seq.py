from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord

DB_V1 = "v1/features.fasta"
DB_V2 = "v2/features.fasta"
OUTPUT_V1 = "v1_features.fasta"
OUTPUT_V2 = "v2_features.fasta"

records: list[SeqRecord] = list()

with open(DB_V2) as handle:
    for record in FastaIO.FastaIterator(handle):
        records.append(SeqRecord(record.seq, id=record.id + "|v2"))

SeqIO.write(records, OUTPUT_V2, "fasta-2line")
records.clear()

with open(DB_V1) as handle:
    for record in FastaIO.FastaIterator(handle):
        records.append(SeqRecord(record.seq, id=record.id + "|v1"))

SeqIO.write(records, OUTPUT_V1, "fasta-2line")