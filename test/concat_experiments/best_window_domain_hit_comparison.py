import os
import pathlib
import warnings

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord

DEEPARG_ENV = "../../deeparg_env"
TARGET_WINDOW = "target_windows.fasta"
INDEXED_FEATURES = "../../data/database/v2/features.dmnd"
ALIGNMENT_OUTPUT_DIR = "alignment_output/"
OUTPUT_COLUMNS = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "qlen",
        "qframe",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "score"
    ]


def get_target_window_alignments(
        orf_windows: pd.DataFrame, orf_sequences: dict) -> pd.DataFrame:
    """Figure out alignments between target windows and database."""
    SeqIO.write(
        sequences=[
            SeqRecord(
                seq=seq[1][1][
                    orf_windows.at[seq[0], "window_start"] :
                    orf_windows.at[seq[0], "window_trunc_end"] + 1],
                id=seq[0],
                description=""
            )
            for seq in orf_sequences.items()
        ],
        handle=TARGET_WINDOW,
        format="fasta"
    )

    # Using a sequence identity of 5% and an e-value of 1000 to allow for
    # scenarios where the concatenated sequence identity and e-value are mostly
    # driven by the source sequence
    og_align = f"{ALIGNMENT_OUTPUT_DIR}target_align.tsv"
    alignment_cmd = f"""
        eval "$(conda shell.bash hook)" && conda activate {DEEPARG_ENV} && \
        diamond blastx -q {TARGET_WINDOW} -d {INDEXED_FEATURES} -k 1000 --id \
        5 -e 1000 -o {og_align} --outfmt 6 {' '.join(OUTPUT_COLUMNS)} \
        --sensitive"""

    if (not pathlib.Path(og_align).exists() and os.system(alignment_cmd) != 0):
        msg = "Diamond didn't run properly"
        raise RuntimeError(msg)

    alignment_results = pd.read_csv(
        filepath_or_buffer=og_align, delimiter="\t", header=None,
        index_col=[0, 1], names=OUTPUT_COLUMNS[2:]
    )

    # Default scoring matrix is BLOSUM62, gap open = 11 and gap extended = 1
    alignment_results["gap_score"] = alignment_results.apply(
        lambda x: (
            int((x["qlen"] - max(x["qstart"], x["qend"])) / 3)
            + (
                10
                if int((x["qlen"] - max(x["qstart"], x["qend"])) / 3) > 0
                else 0
            )
        ),
        axis=1,
    )

    return alignment_results

    # # Traceback ends if it reaches zero
    # return alignment_results.loc[
    #     alignment_results["gap_score"] < alignment_results["score"]
    # ]


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


def get_target_window_domains(orf_windows: pd.DataFrame) -> pd.DataFrame:
    """Figure out conserved domains that overlap with target windows."""
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
    all_hitdata_df["Window Trunc End"] = all_hitdata_df.apply(
        lambda x: orf_windows.at[x.name[0], "window_trunc_end"], axis=1
    )

    return all_hitdata_df


def get_og_source_alignments() -> pd.DataFrame:
    """Get original source sequence alignments along with e-value and raw scores."""  # noqa: E501
    og_align = f"{ALIGNMENT_OUTPUT_DIR}OG_source_align.tsv"
    alignment_cmd = f"""
        eval "$(conda shell.bash hook)" && conda activate {DEEPARG_ENV} && \
        diamond blastx -q ../ORFs.fa -d {INDEXED_FEATURES} -k 1000 \
        --id 10 -e 1 -o {og_align} --outfmt 6 {' '.join(OUTPUT_COLUMNS)} \
        --sensitive"""

    if (not pathlib.Path(og_align).exists() and os.system(alignment_cmd) != 0):
        msg = "Diamond didn't run properly"
        raise RuntimeError(msg)

    return pd.read_csv(
        filepath_or_buffer=og_align, delimiter="\t", header=None,
        index_col=[0, 1], names=OUTPUT_COLUMNS[2:]
    )


def get_frameshifted_source_alignments(
    frame: int,
    direction: str,
    og_source_alignments: pd.DataFrame,
    frame_dir: int,
) -> pd.DataFrame:
    """Get frameshifted source sequence alignments along with e-value and raw scores."""  # noqa: E501
    # Using a sequence identity of 5% and an e-value of 1000 to allow for
    # scenarios where the concatenated sequence identity and e-value are mostly
    # driven by the source sequence
    source_fasta = f"clean_frame_{frame}_{direction}_ORFs.fasta"
    source_align = (
        f"{ALIGNMENT_OUTPUT_DIR}frame_{frame}_{direction}_source_align.tsv"
    )
    alignment_cmd = f"""
        eval "$(conda shell.bash hook)" && conda activate {DEEPARG_ENV} && \
        diamond blastp -q {source_fasta} -d {INDEXED_FEATURES} -k 1000 --id 5 \
        -e 1000 -o {source_align} --outfmt 6 {' '.join(OUTPUT_COLUMNS)} \
        --sensitive"""

    if (
        not pathlib.Path(source_align).exists()
        and os.system(alignment_cmd) != 0
    ):
        msg = "Diamond didn't run properly"
        raise RuntimeError(msg)

    alignment_results = pd.read_csv(
        filepath_or_buffer=source_align, delimiter="\t", header=None,
        index_col=[0, 1], names=OUTPUT_COLUMNS[2:]
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # We only want frameshifted source-ref alignments if original contains
        # a corresponding source-ref alignment
        alignment_results = alignment_results.loc[
            alignment_results.apply(
                lambda x: (
                    False
                    if (x.name[0][:-2], x.name[1])
                    not in og_source_alignments.index
                    else frame_dir
                    not in og_source_alignments.loc[
                        (x.name[0][:-2], x.name[1]), "qframe"
                    ].tolist()
                ),
                axis=1,
            )
        ]

    # Default scoring matrix is BLOSUM62, gap open = 11 and gap extended = 1
    # However, this time we only care about gap extended since gap open is
    # already taken care of with the frameshifting_target_window
    if direction == "FOR":
        alignment_results["gap_score"] = alignment_results.apply(
            lambda x: (x["qstart"] - 1), axis=1,
        )
    else:
        alignment_results["gap_score"] = alignment_results.apply(
            lambda x: (x["qlen"] - x["qend"]), axis=1,
        )

    return alignment_results


def which_target_can_frameshift_source(
    frameshifted_source_alignments: pd.DataFrame,
    target_window_alignments: pd.DataFrame,
    n: int,
    ka_K: np.float64,  # noqa: N803
    ka_lambda: np.float64,
    frame_dir: int
) -> pd.DataFrame:
    """Figure out whether a target window and a source sequence can combine and cause alignment frameshifting.

    Arguments:
        frameshifted_source_alignments: alignments to source sequences that have been frameshifted from original
        target_window_alignments: alignments to target window sequences
        n: length of database
        ka_K: Karlin-Altschul K parameter
        ka_lambda: Karlin-Altschul lambda parameter
        frame_dir: reading frame and direction of alignment

    """  # noqa: E501
    # Let's start off checking if target windows and frameshift sources align
    # to same reference at same reading frame

    target_window_alignments = target_window_alignments.loc[
        target_window_alignments["qframe"] == frame_dir
    ]

    frameshifted_source_alignments["potential_targets"] = (
        frameshifted_source_alignments.apply(
            lambda x: (
                target_window_alignments.loc[pd.IndexSlice[:, x.name[1]], :]
                if x.name[1]
                in target_window_alignments.index.get_level_values(1).to_list()
                else ""
            ),
            axis=1,
        )
    )
    frameshifted_source_alignments = frameshifted_source_alignments.loc[
        ~frameshifted_source_alignments["potential_targets"].isin([""])
    ]
    frameshifted_source_alignments.to_csv(
        "test.csv", mode="a" if frame_dir != 1 else "w"
    )


# The current target windows are not all multiples of three, but what gets
# appended will be. I just don't know how Kamala made that work, but I'll assume
# that she removed the last one or two nucleotides.
orf_windows = pd.read_csv(
    filepath_or_buffer="best_windows.tsv", delimiter="\t", header=0, index_col=0
)
orf_windows["window_trunc_end"] = orf_windows.apply(
    lambda x: (
        int((x["window_end"] - x["window_start"] + 1) / 3) * 3
        + x["window_start"]
        - 1
    ),
    axis=1,
)

orf_sequences = {
    rec.id: (len(rec), rec.seq)
    for rec in FastaIO.FastaIterator(source="clean_ORFs.fa")
    if rec.id in orf_windows.index.to_list()
}
orf_windows["length"] = orf_windows.apply(
    lambda x: orf_sequences[x.name][0], axis=1
)
target_window_alignments = get_target_window_alignments(
    orf_windows, orf_sequences
)
target_window_domains = get_target_window_domains(orf_windows)
og_source_alignments = get_og_source_alignments()

# Estimate Karlin-Altschul parameters from og_source_alignments
# E-value is calculated by this formula: e = ka_K * m * n * exp(-ka_lambda * s)
# e             E-value
# m             Length of query sequence
# n             Length of database
# s             Raw score
# ka_K          Karlin-Altschul K parameter
# ka_lambda     Karlin-Altschul lambda parameter
e = og_source_alignments["evalue"].to_numpy()
m = np.floor_divide(og_source_alignments["qlen"].to_numpy(), 3)
n = sum(len(rec) for rec in FastaIO.FastaIterator(source="clean_ORFs.fa"))
s = og_source_alignments["score"].to_numpy()

# Can find lambda via ln((e[0] * m[1]) / (e[1] * m[0])) / (s[1] - s[0])
e0 = e[0]; m0 = m[0]; s0 = s[0]  # noqa: E702
e1 = e[1:]; m1 = m[1:]; s1 = s[1:]  # noqa: E702
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ka_lambda = np.nanmean(
        np.ma.masked_invalid(np.log((e0 * m1) / (e1 * m0)) / (s1 - s0))
    )
    ka_K = np.nanmean(e / (m * n * np.exp(-1 * ka_lambda * s)))  # noqa: N816

target_window_alignments.to_csv("target_window_alignments.csv")
target_window_domains.to_csv("target_window_domains.csv")
og_source_alignments.to_csv("og_source_alignments.csv")
for direction in ["FOR", "REV"]:
    for frame in range(1, 4):
        # In transeq, reverse reading frames are the reverse-complement of the
        # forward reading frame. In DIAMOND, reverse reading frames are the forward
        # reading frames of the reverse-complement of the original sequence.
        # Thus, -2 in transeq is -3 in DIAMOND, and -2 in DIAMOND is -3 in transeq
        frame_dir = (
            frame
            if direction == "FOR"
            else -1
            if frame == 1
            else -2
            if frame == 3  # noqa: PLR2004
            else -3
        )
        frameshifted_source_alignments = get_frameshifted_source_alignments(
            frame, direction, og_source_alignments, frame_dir
        )
        frameshifted_source_alignments.to_csv(
            f"frame_{frame}_{direction}_source_align.csv"
        )
        which_target_can_frameshift_source(
            frameshifted_source_alignments,
            target_window_alignments,
            n,
            ka_K,
            ka_lambda,
            frame_dir,
        )
