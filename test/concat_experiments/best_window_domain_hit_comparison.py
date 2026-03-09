import os
import pathlib

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord

DEEPARG_ENV = "../../deeparg_env"
TARGET_WINDOW = "target_windows.fasta"
INDEXED_FEATURES = "../../data/database/v2/features.dmnd"
ALIGNMENT_OUTPUT_DIR = "alignment_output/"


def get_potential_frameshifting_target_windows(
        orf_windows: pd.DataFrame, orf_sequences: dict) -> pd.DataFrame:
    """Figure out which target windows can change the reading frame of the source sequence."""  # noqa: E501
    SeqIO.write(
        sequences=[
            SeqRecord(
                seq=seq[1][1][
                    orf_windows.at[seq[0], "window_start"] :
                    orf_windows.at[seq[0], "window_end"] + 1],
                id=seq[0],
                description=""
            )
            for seq in orf_sequences.items()
        ],
        handle=TARGET_WINDOW,
        format="fasta"
    )

    output_columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "qlen",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "score"
    ]

    og_align = f"{ALIGNMENT_OUTPUT_DIR}/OG_align.tsv"

    # Using a sequence identity of 5% and an e-value of 1000 to allow for
    # scenarios where the concatenated sequence identity and e-value are mostly
    # driven by the source sequence
    alignment_cmd = f"""
        eval "$(conda shell.bash hook)" && conda activate {DEEPARG_ENV} && \
        diamond blastx -q {TARGET_WINDOW} -d {INDEXED_FEATURES} -k 1000 --id \
        5 -e 1000 -o {og_align} --outfmt 6 {' '.join(output_columns)} \
        --sensitive"""

    if (not pathlib.Path(og_align).exists() and os.system(alignment_cmd) != 0):
        msg = "Diamond didn't run properly"
        raise RuntimeError(msg)

    alignment_results = pd.read_csv(
        filepath_or_buffer=og_align, delimiter="\t", header=None,
        index_col=[0, 1], names=output_columns[2:]
    )

    alignment_results = alignment_results.loc[
        (
            alignment_results["qlen"]
            - (np.maximum(
                alignment_results["qstart"], alignment_results["qend"]
            ))
        ) % 3 != 0
    ]

    # Default scoring matrix is BLOSUM62, gap open = 11 and gap extended = 1
    return alignment_results.loc[
        (
            alignment_results["qlen"]
            - (np.maximum(
                alignment_results["qstart"], alignment_results["qend"]
            ))
        ) + 11 < alignment_results["score"]
    ]


def get_from_nuc(index: str, from_codon: int) -> int:
    """Get the corresponding nucleotide index for the first codon of a domain."""  # noqa: E501
    name = index[0]
    frame = int(index[1].split("_")[0])  # noqa: PLC0207
    direction = index[1].split("_")[1]
    if direction == "FOR":
        return from_codon * 3 + frame - 1
    return orf_windows.at[name, "length"] - (from_codon * 3 + frame - 1)


def get_to_nuc(index: str, to_codon: int) -> int:
    """Get the corresponding nucleotide index for the last codon of a domain."""
    name = index[0]
    frame = int(index[1].split("_")[0])  # noqa: PLC0207
    direction = index[1].split("_")[1]
    if direction == "FOR":
        return to_codon * 3 + frame + 2
    return (orf_windows.at[name, "length"] - frame + 1) - to_codon * 3


orf_windows = pd.read_csv(
    filepath_or_buffer="best_windows.tsv", delimiter="\t", header=0, index_col=0
)

orf_sequences = {
    rec.id: (len(rec), rec.seq)
    for rec in FastaIO.FastaIterator(source="clean_ORFs.fa")
    if rec.id in orf_windows.index.to_list()
}

orf_windows["length"] = orf_windows.apply(
    lambda x: orf_sequences[x.name][0], axis=1
)

target_window = get_potential_frameshifting_target_windows(
    orf_windows, orf_sequences
)

all_hitdata: list[dict] = []
for direction in ["FOR", "REV"]:
    for frame in range(1, 4):
        hitdata = pd.read_csv(
            filepath_or_buffer=f"frame_{frame}_{direction}_hitdata.tsv",
            delimiter="\t",
            header=0,
        )
        hitdata["Query"] = hitdata["Query"].apply(
            lambda x: x.split(" ")[2][1:-2]
        )
        hitdata["Frame"] = f"{frame}_{direction}"
        all_hitdata.extend(hitdata.to_dict(orient="records"))

all_hitdata_df = pd.DataFrame(data=all_hitdata)
all_hitdata_df = all_hitdata_df.loc[
    all_hitdata_df["Query"].isin(orf_windows.index.to_list())
]
all_hitdata_df = all_hitdata_df.set_index(keys=["Query", "Frame"])

all_hitdata_df = all_hitdata_df.drop(
    columns=["Hit type", "PSSM-ID", "E-Value", "Bitscore", "Incomplete"]
)

all_hitdata_df["From Nuc"] = all_hitdata_df.apply(
    lambda x: get_from_nuc(x.name, x["From"]), axis=1
)
all_hitdata_df["To Nuc"] = all_hitdata_df.apply(
    lambda x: get_to_nuc(x.name, x["To"]), axis=1
)
all_hitdata_df["Window Start"] = all_hitdata_df.apply(
    lambda x: orf_windows.at[x.name[0], "window_start"], axis=1
)
all_hitdata_df["Window End"] = all_hitdata_df.apply(
    lambda x: orf_windows.at[x.name[0], "window_end"], axis=1
)

all_hitdata_df.to_csv("test.csv")
