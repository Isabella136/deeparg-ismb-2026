import pathlib

from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord

DB_V1 = "../data/database/v1/features.fasta"
DB_V2 = "../data/database/v2/features.fasta"
OUTPUT_V1 = "v1_features.fasta"
OUTPUT_V2 = "v2_features.fasta"

records: list[SeqRecord] = []

with pathlib.Path(DB_V2).open(encoding="utf-8") as handle:
    records.extend(
        SeqRecord(record.seq, id=record.id + "|v2")
        for record in FastaIO.FastaIterator(handle)
    )

SeqIO.write(records, OUTPUT_V2, "fasta-2line")
records.clear()

with pathlib.Path(DB_V1).open(encoding="utf-8") as handle:
    records.extend(
        SeqRecord(record.seq, id=record.id + "|v1")
        for record in FastaIO.FastaIterator(handle)
    )

SeqIO.write(records, OUTPUT_V1, "fasta-2line")
