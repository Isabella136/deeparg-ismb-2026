from itertools import zip_longest

import numpy as np
import pandas as pd

DEEPARG_HIT_FILE = "X.mapping.ARG"
ALIGNMENT_FILE = "max_class_hit_per_query.csv"
PATH_TO_SS = "deeparg_results"
PATH_TO_LS = "spades/deeparg_results"


# Define helper functions
def compare_freq_with_deep(
    df: pd.DataFrame, idx_list: pd.MultiIndex, al_ids: int | list[int]
) -> list[bool]:
    """Compare most frequent class to DeepARG class for each query."""
    return [
        df.at[(idx[0], idx[1], idx[2], al_id), "Most Frequent Class"]
        == df.at[(idx[0], idx[1], idx[2], al_id), "DeepARG Class"]
        for idx, al_id in zip_longest(idx_list, al_ids, fillvalue=al_ids[0])
    ]


def get_freq_class_column(
    df: pd.DataFrame, idx_list: pd.MultiIndex, al_ids: int | list[int]
) -> list[str]:
    """Get most frequent class for each query."""
    return [
        df.at[(idx[0], idx[1], idx[2], al_id), "Most Frequent Class"]
        for idx, al_id in zip_longest(idx_list, al_ids, fillvalue=al_ids[0])
    ]


def get_deep_class_column(
    df: pd.DataFrame, idx_list: pd.MultiIndex, al_ids: int | list[int]
) -> list[str]:
    """Get DeepARG class for each query."""
    return [
        df.at[(idx[0], idx[1], idx[2], al_id), "DeepARG Class"]
        for idx, al_id in zip_longest(idx_list, al_ids, fillvalue=al_ids[0])
    ]


def get_other_threshold_info(df: pd.DataFrame, al_id: int) -> pd.DataFrame:
    """Add class information for al_id% threshold."""
    df[f"{al_id} pair-containing"] = False
    df[f"{al_id} Most frequent class"] = ""
    previous_pred_file = ""
    prediction_df = pd.DataFrame()
    prev_max_class_file = ""
    max_class_df = pd.DataFrame()
    for idx in df.index.to_list():
        sample, model, query = idx
        path_to_file = (
            f"../samples/{sample}/{PATH_TO_SS if model == 'SS' else PATH_TO_LS}"
        )
        file_to_open = (
            f"{path_to_file}/arg_alignment_identity_{al_id}/{DEEPARG_HIT_FILE}"
        )
        if previous_pred_file != file_to_open:
            previous_pred_file = file_to_open
            prediction_df = pd.read_csv(file_to_open, sep="\t", header=0, index_col=3)
        df.at[(sample, model, query), f"{al_id} predicted"] = (
            query in prediction_df.index.to_list()
        )
        df.at[(sample, model, query), f"{al_id} DeepARG class"] = (
            prediction_df.at[query, "predicted_ARG-class"]
            if df.at[(sample, model, query), f"{al_id} predicted"]
            else ""
        )
        file_to_open = (
            f"{path_to_file}/arg_alignment_identity_{al_id}/{ALIGNMENT_FILE}"
        )
        if prev_max_class_file != file_to_open:
            prev_max_class_file = file_to_open
            max_class_df = pd.read_csv(file_to_open, header=0, index_col=0)
        df.at[(sample, model, query), f"{al_id} aligned"] = (
            query in max_class_df.index.to_list()
        )
        df.at[(sample, model, query), f"{al_id} Most frequent class"] = (
            (
                frozenset(max_class_df.loc[query]["deeparg hit"].to_list())
                if np.sum(max_class_df.index.isin([query])) > 1
                else frozenset({max_class_df.at[query, "deeparg hit"]})
            )
            if df.at[(sample, model, query), f"{al_id} aligned"]
            else frozenset()
        )
    return df


# Get list of all pair-containing queries
label_count = pd.read_csv("label_counts.tsv", sep="\t", header=0)
pair_containing_queries = label_count[
    [
        "Sample ID",
        "Alignment Identity",
        "Model",
        "Most Frequent Class",
        "DeepARG Class",
        "Query",
    ]
].drop_duplicates()
del label_count
pair_containing_queries = pair_containing_queries.sort_values(
    by=["Query", "Alignment Identity"]
)
pair_containing_queries = pair_containing_queries.set_index(
    keys=["Query", "Alignment Identity"]
)

# Keep queries that class switched at least once
class_switch_queries = pair_containing_queries[
    pair_containing_queries.index.get_level_values(0).isin(
        pair_containing_queries[
            pair_containing_queries["Most Frequent Class"]
            != pair_containing_queries["DeepARG Class"]
        ]
        .index.get_level_values(0)
        .to_numpy(),
    )
]
del pair_containing_queries

# Calculate the number of class switches per alignment identity and model
class_switch_count = (
    class_switch_queries
    .loc[
        class_switch_queries["Most Frequent Class"]
        != class_switch_queries["DeepARG Class"]
    ]
    .groupby(["Alignment Identity", "Model"])
    .count()["Sample ID"]
)
class_switch_count.to_csv("class_switch_count.csv")
# Remove queries that were pair-containing for all three thresholds and
# that had same most frequent class for all three thresholds as well
thrice_pair_containing_queries = class_switch_queries[
    class_switch_queries.index.get_level_values(0).isin(
        class_switch_queries.index
        .get_level_values(0)
        .drop_duplicates()
        .to_numpy()[
            class_switch_queries.index
            .get_level_values(0)
            .value_counts(sort=False)
            .eq(3)
        ],
    )
]
thrice_pair_containing_queries_index = (
    thrice_pair_containing_queries.index
    .get_level_values(0)
    .drop_duplicates()
    .to_numpy()
)
thrice_class_switch_queries_index = thrice_pair_containing_queries_index[
    thrice_pair_containing_queries
    .reset_index(level=0)
    .drop_duplicates()
    .reset_index()
    .set_index(keys=["Query", "Alignment Identity"])
    .index.get_level_values(0)
    .value_counts(sort=False)
    .eq(1)
]
not_thrice_class_switch_queries = class_switch_queries[
    ~class_switch_queries.index.get_level_values(0).isin(
        thrice_class_switch_queries_index
    )
]
del (
    class_switch_queries,
    thrice_class_switch_queries_index,
    thrice_pair_containing_queries_index,
    thrice_pair_containing_queries,
)

# Reformat query names
not_thrice_class_switch_queries = not_thrice_class_switch_queries.reset_index()
not_thrice_class_switch_queries["Query"] = [
    query[13:] for query in not_thrice_class_switch_queries["Query"].to_list()
]
not_thrice_class_switch_queries = not_thrice_class_switch_queries.set_index(
    keys=["Sample ID", "Model", "Query"]
)

# Let's start with queries that were pair-containing for all three thresholds
# but didn't have the same most DeepARG or frequent class for all thresholds
thrice_pair_containing_not_same_class_switch_queries = not_thrice_class_switch_queries[
    not_thrice_class_switch_queries.index.isin(
        not_thrice_class_switch_queries.index.drop_duplicates().to_numpy()[
            not_thrice_class_switch_queries.index.value_counts(sort=False).eq(3)
        ]
    )
]
thrice_pair_containing_not_same_class_switch_queries = (
    thrice_pair_containing_not_same_class_switch_queries.set_index(
        keys="Alignment Identity", append=True
    )
)
thrice_pair_containing_not_same_class_switch_queries_info_index: pd.MultiIndex = (
    thrice_pair_containing_not_same_class_switch_queries.index.droplevel(
        3
    ).drop_duplicates()
)
thrice_pair_containing_not_same_class_switch_queries_info = pd.DataFrame(
    data={
        "30 same": compare_freq_with_deep(
            df=thrice_pair_containing_not_same_class_switch_queries,
            idx_list=thrice_pair_containing_not_same_class_switch_queries_info_index.to_list(),
            al_ids=[30],
        ),
        "30 pair-containing": True,
        "30 predicted": True,
        "30 aligned": True,
        "30 Most frequent class": get_freq_class_column(
            df=thrice_pair_containing_not_same_class_switch_queries,
            idx_list=thrice_pair_containing_not_same_class_switch_queries_info_index.to_list(),
            al_ids=[30],
        ),
        "30 DeepARG class": get_deep_class_column(
            df=thrice_pair_containing_not_same_class_switch_queries,
            idx_list=thrice_pair_containing_not_same_class_switch_queries_info_index.to_list(),
            al_ids=[30],
        ),
        "50 same": compare_freq_with_deep(
            df=thrice_pair_containing_not_same_class_switch_queries,
            idx_list=thrice_pair_containing_not_same_class_switch_queries_info_index.to_list(),
            al_ids=[50],
        ),
        "50 pair-containing": True,
        "50 predicted": True,
        "50 aligned": True,
        "50 Most frequent class": get_freq_class_column(
            df=thrice_pair_containing_not_same_class_switch_queries,
            idx_list=thrice_pair_containing_not_same_class_switch_queries_info_index.to_list(),
            al_ids=[50],
        ),
        "50 DeepARG class": get_deep_class_column(
            df=thrice_pair_containing_not_same_class_switch_queries,
            idx_list=thrice_pair_containing_not_same_class_switch_queries_info_index.to_list(),
            al_ids=[50],
        ),
        "80 same": compare_freq_with_deep(
            df=thrice_pair_containing_not_same_class_switch_queries,
            idx_list=thrice_pair_containing_not_same_class_switch_queries_info_index.to_list(),
            al_ids=[80],
        ),
        "80 pair-containing": True,
        "80 predicted": True,
        "80 aligned": True,
        "80 Most frequent class": get_freq_class_column(
            df=thrice_pair_containing_not_same_class_switch_queries,
            idx_list=thrice_pair_containing_not_same_class_switch_queries_info_index.to_list(),
            al_ids=[80],
        ),
        "80 DeepARG class": get_deep_class_column(
            df=thrice_pair_containing_not_same_class_switch_queries,
            idx_list=thrice_pair_containing_not_same_class_switch_queries_info_index.to_list(),
            al_ids=[80],
        ),
    },
    index=thrice_pair_containing_not_same_class_switch_queries_info_index,
)
# Let's narrow on queries where 30 and 50 are the same class switches
# but 80 is pair containing
same_class_switch_30_50_pair_containing_80_info = (
    thrice_pair_containing_not_same_class_switch_queries_info[
        ~thrice_pair_containing_not_same_class_switch_queries_info["30 same"]
        & ~thrice_pair_containing_not_same_class_switch_queries_info["50 same"]
        & thrice_pair_containing_not_same_class_switch_queries_info["80 same"]
    ]
)
same_class_switch_30_50_pair_containing_80_info = (
    same_class_switch_30_50_pair_containing_80_info.drop(
        columns=["30 same", "50 same", "80 same"]
    )
)
same_class_switch_30_50_pair_containing_80_info = (
    same_class_switch_30_50_pair_containing_80_info.reindex(
        columns=sorted(same_class_switch_30_50_pair_containing_80_info.columns)
    )
)

# Let's narrow on queries where 30 is a class switch
# but 50 and 80 are pair containing
class_switch_30_pair_containing_50_80_info = (
    thrice_pair_containing_not_same_class_switch_queries_info[
        ~thrice_pair_containing_not_same_class_switch_queries_info["30 same"]
        & thrice_pair_containing_not_same_class_switch_queries_info["50 same"]
        & thrice_pair_containing_not_same_class_switch_queries_info["80 same"]
    ]
)
class_switch_30_pair_containing_50_80_info = (
    class_switch_30_pair_containing_50_80_info.drop(
        columns=["30 same", "50 same", "80 same"]
    )
)
class_switch_30_pair_containing_50_80_info = (
    class_switch_30_pair_containing_50_80_info.reindex(
        columns=sorted(class_switch_30_pair_containing_50_80_info.columns)
    )
)

# Let's narrow on queries where 30 and 50 are pair containing
# but 80 is a class switch
pair_containing_30_50_class_switch_80_info = (
    thrice_pair_containing_not_same_class_switch_queries_info[
        thrice_pair_containing_not_same_class_switch_queries_info["30 same"]
        & thrice_pair_containing_not_same_class_switch_queries_info["50 same"]
        & ~thrice_pair_containing_not_same_class_switch_queries_info["80 same"]
    ]
)
pair_containing_30_50_class_switch_80_info = (
    pair_containing_30_50_class_switch_80_info.drop(
        columns=["30 same", "50 same", "80 same"]
    )
)
pair_containing_30_50_class_switch_80_info = (
    pair_containing_30_50_class_switch_80_info.reindex(
        columns=sorted(pair_containing_30_50_class_switch_80_info.columns)
    )
)

# Let's narrow on queries where 30 and 80 are pair containing
# but 50 is a class switch
pair_containing_30_class_switch_50_pair_containing_80_info = (
    thrice_pair_containing_not_same_class_switch_queries_info[
        thrice_pair_containing_not_same_class_switch_queries_info["30 same"]
        & ~thrice_pair_containing_not_same_class_switch_queries_info["50 same"]
        & thrice_pair_containing_not_same_class_switch_queries_info["80 same"]
    ]
)
pair_containing_30_class_switch_50_pair_containing_80_info = (
    pair_containing_30_class_switch_50_pair_containing_80_info.drop(
        columns=["30 same", "50 same", "80 same"]
    )
)
pair_containing_30_class_switch_50_pair_containing_80_info = (
    pair_containing_30_class_switch_50_pair_containing_80_info.reindex(
        columns=sorted(
            pair_containing_30_class_switch_50_pair_containing_80_info.columns
        )
    )
)

# Let's filter out queries that were pair-containing for all three thresholds
# but didn't have the same most DeepARG or frequent class for all thresholds
not_thrice_pair_containing_queries = not_thrice_class_switch_queries[
    ~not_thrice_class_switch_queries.index.isin(
        thrice_pair_containing_not_same_class_switch_queries_info_index
    )
]
del (
    not_thrice_class_switch_queries,
    thrice_pair_containing_not_same_class_switch_queries_info_index,
    thrice_pair_containing_not_same_class_switch_queries_info,
    thrice_pair_containing_not_same_class_switch_queries,
)

# Now, let's do queries that were pair-containing for two thresholds
# but didn't have the same most DeepARG or frequent class for those thresholds
twice_pair_containing_queries = not_thrice_pair_containing_queries[
    not_thrice_pair_containing_queries.index.isin(
        not_thrice_pair_containing_queries.index.drop_duplicates().to_numpy()[
            not_thrice_pair_containing_queries
            .reset_index()[
                [
                    "Sample ID",
                    "Model",
                    "Query",
                    "Most Frequent Class",
                    "DeepARG Class",
                ]
            ]
            .drop_duplicates()
            .set_index(keys=["Sample ID", "Model", "Query"])
            .index.value_counts(sort=False)
            .eq(2)
        ]
    )
]
twice_pair_containing_queries_first_al_id = twice_pair_containing_queries[
    ~twice_pair_containing_queries.index.duplicated(keep="first")
]["Alignment Identity"].to_list()
twice_pair_containing_queries_last_al_id = twice_pair_containing_queries[
    ~twice_pair_containing_queries.index.duplicated(keep="last")
]["Alignment Identity"].to_list()
twice_pair_containing_queries = twice_pair_containing_queries.set_index(
    keys="Alignment Identity", append=True
)
twice_pair_containing_queries_info_index: pd.MultiIndex = (
    twice_pair_containing_queries.index.droplevel(3).drop_duplicates()
)
twice_pair_containing_queries_info = pd.DataFrame(
    data={
        "first identity": twice_pair_containing_queries_first_al_id,
        "first identity same": compare_freq_with_deep(
            df=twice_pair_containing_queries,
            idx_list=twice_pair_containing_queries_info_index.to_list(),
            al_ids=twice_pair_containing_queries_first_al_id,
        ),
        "first identity Most frequent class": get_freq_class_column(
            df=twice_pair_containing_queries,
            idx_list=twice_pair_containing_queries_info_index.to_list(),
            al_ids=twice_pair_containing_queries_first_al_id,
        ),
        "first identity DeepARG class": get_deep_class_column(
            df=twice_pair_containing_queries,
            idx_list=twice_pair_containing_queries_info_index.to_list(),
            al_ids=twice_pair_containing_queries_first_al_id,
        ),
        "last identity": twice_pair_containing_queries_last_al_id,
        "last identity same": compare_freq_with_deep(
            df=twice_pair_containing_queries,
            idx_list=twice_pair_containing_queries_info_index.to_list(),
            al_ids=twice_pair_containing_queries_last_al_id,
        ),
        "last identity Most frequent class": get_freq_class_column(
            df=twice_pair_containing_queries,
            idx_list=twice_pair_containing_queries_info_index.to_list(),
            al_ids=twice_pair_containing_queries_last_al_id,
        ),
        "last identity DeepARG class": get_deep_class_column(
            df=twice_pair_containing_queries,
            idx_list=twice_pair_containing_queries_info_index.to_list(),
            al_ids=twice_pair_containing_queries_last_al_id,
        ),
    },
    index=twice_pair_containing_queries_info_index,
)

# Let's specifically look at queries that were pair-containing at 50% and
# class switched at 80%.
pair_containing_50_class_switch_80_info = twice_pair_containing_queries_info.loc[
    twice_pair_containing_queries_info["last identity"] == 80  # noqa: PLR2004
]
pair_containing_50_class_switch_80_info = (
    pair_containing_50_class_switch_80_info.rename(
        columns={
            "first identity Most frequent class": "50 Most frequent class",
            "first identity DeepARG class": "50 DeepARG class",
            "last identity Most frequent class": "80 Most frequent class",
            "last identity DeepARG class": "80 DeepARG class",
        }
    )
)
pair_containing_50_class_switch_80_info["50 pair-containing"] = True
pair_containing_50_class_switch_80_info["80 pair-containing"] = True
pair_containing_50_class_switch_80_info["50 predicted"] = True
pair_containing_50_class_switch_80_info["80 predicted"] = True
pair_containing_50_class_switch_80_info["50 aligned"] = True
pair_containing_50_class_switch_80_info["80 aligned"] = True
pair_containing_50_class_switch_80_info = pair_containing_50_class_switch_80_info.drop(
    columns=[
        "first identity",
        "last identity",
        "first identity same",
        "last identity same",
    ]
)
pair_containing_50_class_switch_80_info = get_other_threshold_info(
    pair_containing_50_class_switch_80_info, 30
)
pair_containing_50_class_switch_80_info = (
    pair_containing_50_class_switch_80_info.reindex(
        columns=sorted(pair_containing_50_class_switch_80_info.columns)
    )
)

# Now, let's specifically look at queries that were pair-containing at
# 30% and 50% but didn't have the same class switch
pair_containing_not_same_class_switch_30_50_info = (
    twice_pair_containing_queries_info.loc[
        twice_pair_containing_queries_info["last identity"] != 80  # noqa: PLR2004
    ]
)
del (
    twice_pair_containing_queries_info,
    twice_pair_containing_queries,
    twice_pair_containing_queries_first_al_id,
    twice_pair_containing_queries_last_al_id,
)
pair_containing_not_same_class_switch_30_50_info = (
    pair_containing_not_same_class_switch_30_50_info.rename(
        columns={
            "first identity same": "30 same",
            "first identity Most frequent class": "30 Most frequent class",
            "first identity DeepARG class": "30 DeepARG class",
            "last identity same": "50 same",
            "last identity Most frequent class": "50 Most frequent class",
            "last identity DeepARG class": "50 DeepARG class",
        }
    )
)
pair_containing_not_same_class_switch_30_50_info["30 pair-containing"] = True
pair_containing_not_same_class_switch_30_50_info["50 pair-containing"] = True
pair_containing_not_same_class_switch_30_50_info["30 predicted"] = True
pair_containing_not_same_class_switch_30_50_info["50 predicted"] = True
pair_containing_not_same_class_switch_30_50_info["30 aligned"] = True
pair_containing_not_same_class_switch_30_50_info["50 aligned"] = True
pair_containing_not_same_class_switch_30_50_info = (
    pair_containing_not_same_class_switch_30_50_info.drop(
        columns=["first identity", "last identity"]
    )
)
pair_containing_not_same_class_switch_30_50_info = get_other_threshold_info(
    pair_containing_not_same_class_switch_30_50_info, 80
)

# Narrow on queries that were pair-containing at 30% and class switched at 50%
pair_containing_30_class_switch_50_info = (
    pair_containing_not_same_class_switch_30_50_info[
        pair_containing_not_same_class_switch_30_50_info["30 same"]
    ]
)
pair_containing_30_class_switch_50_info = pair_containing_30_class_switch_50_info.drop(
    columns=["30 same", "50 same"]
)
pair_containing_30_class_switch_50_info = (
    pair_containing_30_class_switch_50_info.reindex(
        columns=sorted(pair_containing_30_class_switch_50_info.columns)
    )
)

# Narrow on queries that class switched at 30% and were pair-containing at 50%
class_switch_30_pair_containing_50_info = (
    pair_containing_not_same_class_switch_30_50_info[
        pair_containing_not_same_class_switch_30_50_info["50 same"]
    ]
)
class_switch_30_pair_containing_50_info = class_switch_30_pair_containing_50_info.drop(
    columns=["30 same", "50 same"]
)
class_switch_30_pair_containing_50_info = (
    class_switch_30_pair_containing_50_info.reindex(
        columns=sorted(class_switch_30_pair_containing_50_info.columns)
    )
)

# Narrow on queries that class switched at 30% and 50%
# but class switches weren't the same
not_same_class_switch_30_50_info = pair_containing_not_same_class_switch_30_50_info[
    pair_containing_not_same_class_switch_30_50_info["30 same"]
    == pair_containing_not_same_class_switch_30_50_info["50 same"]
]
del pair_containing_not_same_class_switch_30_50_info
not_same_class_switch_30_50_info = not_same_class_switch_30_50_info.drop(
    columns=["30 same", "50 same"]
)
not_same_class_switch_30_50_info = not_same_class_switch_30_50_info.reindex(
    columns=sorted(not_same_class_switch_30_50_info.columns)
)

# Let's filter out queries that were pair-containing for two thresholds
# but didn't have the same most DeepARG or frequent class for those thresholds
once_or_twice_same_class_switch_queries = not_thrice_pair_containing_queries[
    ~not_thrice_pair_containing_queries.index.isin(
        twice_pair_containing_queries_info_index
    )
]
del twice_pair_containing_queries_info_index, not_thrice_pair_containing_queries

# Now, let's do queries that had the same class switch for two thresholds
twice_same_class_switch_queries = once_or_twice_same_class_switch_queries[
    once_or_twice_same_class_switch_queries.index.isin(
        once_or_twice_same_class_switch_queries.index.drop_duplicates().to_numpy()[
            once_or_twice_same_class_switch_queries.index.value_counts(sort=False).eq(
                2
            )
        ]
    )
]
twice_same_class_switch_queries_first_al_id = (
    twice_same_class_switch_queries[
        ~twice_same_class_switch_queries.index.duplicated(keep="first")
    ]["Alignment Identity"]
).to_list()
twice_same_class_switch_queries_last_al_id = (
    twice_same_class_switch_queries[
        ~twice_same_class_switch_queries.index.duplicated(keep="last")
    ]["Alignment Identity"]
).to_list()
twice_same_class_switch_queries = twice_same_class_switch_queries.set_index(
    keys="Alignment Identity", append=True
)
twice_same_class_switch_queries_info_index: pd.MultiIndex = (
    twice_same_class_switch_queries.index.droplevel(3).drop_duplicates()
)
twice_same_class_switch_queries_info = pd.DataFrame(
    data={
        "first identity": twice_same_class_switch_queries_first_al_id,
        "last identity": twice_same_class_switch_queries_last_al_id,
        "Most frequent class": get_freq_class_column(
            df=twice_same_class_switch_queries,
            idx_list=twice_same_class_switch_queries_info_index.to_list(),
            al_ids=twice_same_class_switch_queries_first_al_id,
        ),
        "DeepARG class": get_deep_class_column(
            df=twice_same_class_switch_queries,
            idx_list=twice_same_class_switch_queries_info_index.to_list(),
            al_ids=twice_same_class_switch_queries_first_al_id,
        ),
    },
    index=twice_same_class_switch_queries_info_index,
)

# All queries that have the same class switch for two thresholds
# have it for 30 and 50%
same_class_switch_30_50_info = twice_same_class_switch_queries_info.loc[
    twice_same_class_switch_queries_info["last identity"] != 80  # noqa: PLR2004
]
del (
    twice_same_class_switch_queries_info,
    twice_same_class_switch_queries,
    twice_same_class_switch_queries_first_al_id,
    twice_same_class_switch_queries_last_al_id,
)
same_class_switch_30_50_info = same_class_switch_30_50_info.rename(
    columns={
        "Most frequent class": "30 Most frequent class",
        "DeepARG class": "30 DeepARG class",
    }
)
same_class_switch_30_50_info["30 pair-containing"] = True
same_class_switch_30_50_info["50 pair-containing"] = True
same_class_switch_30_50_info["30 predicted"] = True
same_class_switch_30_50_info["50 predicted"] = True
same_class_switch_30_50_info["30 aligned"] = True
same_class_switch_30_50_info["50 aligned"] = True
same_class_switch_30_50_info["50 Most frequent class"] = same_class_switch_30_50_info[
    "30 Most frequent class"
]
same_class_switch_30_50_info["50 DeepARG class"] = same_class_switch_30_50_info[
    "30 DeepARG class"
]
same_class_switch_30_50_info = same_class_switch_30_50_info.drop(
    columns=["first identity", "last identity"]
)
same_class_switch_30_50_info = get_other_threshold_info(
    same_class_switch_30_50_info, 80
)
same_class_switch_30_50_info = same_class_switch_30_50_info.reindex(
    columns=sorted(same_class_switch_30_50_info.columns)
)

# Let's only keep queries that class switched for one threshold
once_class_switch_queries = once_or_twice_same_class_switch_queries[
    ~once_or_twice_same_class_switch_queries.index.isin(
        twice_same_class_switch_queries_info_index
    )
]
once_class_switch_queries["Most Frequent Class"] = once_class_switch_queries[
    "Most Frequent Class"
].map(lambda most_freq_class: most_freq_class)
del (
    twice_same_class_switch_queries_info_index,
    once_or_twice_same_class_switch_queries,
)

# Let's look at queries that class switched at 30
class_switch_30_info = once_class_switch_queries.loc[
    once_class_switch_queries["Alignment Identity"] == 30  # noqa: PLR2004
]
class_switch_30_info = class_switch_30_info.rename(
    columns={
        "Most Frequent Class": "30 Most frequent class",
        "DeepARG Class": "30 DeepARG class",
    }
)
class_switch_30_info["30 pair-containing"] = True
class_switch_30_info["30 predicted"] = True
class_switch_30_info["30 aligned"] = True
class_switch_30_info = class_switch_30_info.drop(columns=["Alignment Identity"])
class_switch_30_info = get_other_threshold_info(class_switch_30_info, 50)
class_switch_30_info = get_other_threshold_info(class_switch_30_info, 80)
class_switch_30_info = class_switch_30_info.reindex(
    columns=sorted(class_switch_30_info.columns)
)

# Let's look at queries that class switched at 50
class_switch_50_info = once_class_switch_queries.loc[
    once_class_switch_queries["Alignment Identity"] == 50  # noqa: PLR2004
]
class_switch_50_info = class_switch_50_info.rename(
    columns={
        "Most Frequent Class": "50 Most frequent class",
        "DeepARG Class": "50 DeepARG class",
    }
)
class_switch_50_info["50 pair-containing"] = True
class_switch_50_info["50 predicted"] = True
class_switch_50_info["50 aligned"] = True
class_switch_50_info = class_switch_50_info.drop(columns=["Alignment Identity"])
class_switch_50_info = get_other_threshold_info(class_switch_50_info, 30)
class_switch_50_info = get_other_threshold_info(class_switch_50_info, 80)
class_switch_50_info = class_switch_50_info.reindex(
    columns=sorted(class_switch_50_info.columns)
)

# Let's look at queries that class switched at 80
class_switch_80_info = once_class_switch_queries.loc[
    once_class_switch_queries["Alignment Identity"] == 80  # noqa: PLR2004
]
del once_class_switch_queries
class_switch_80_info = class_switch_80_info.rename(
    columns={
        "Most Frequent Class": "80 Most frequent class",
        "DeepARG Class": "80 DeepARG class",
    }
)
class_switch_80_info["80 pair-containing"] = True
class_switch_80_info["80 predicted"] = True
class_switch_80_info["80 aligned"] = True
class_switch_80_info = class_switch_80_info.drop(columns=["Alignment Identity"])
class_switch_80_info = get_other_threshold_info(class_switch_80_info, 30)
class_switch_80_info = get_other_threshold_info(class_switch_80_info, 50)
class_switch_80_info = class_switch_80_info.reindex(
    columns=sorted(class_switch_80_info.columns)
)

# Concat all queries that only class switched at 30
class_switch_30_concat_info = pd.concat(
    objs=[
        class_switch_30_info,
        class_switch_30_pair_containing_50_info,
        class_switch_30_pair_containing_50_80_info,
    ]
)

is_class_in_list = np.vectorize(
    lambda deep_class, frequent_class_list: deep_class in frequent_class_list
)
do_lists_overlap = np.vectorize(
    lambda list_a, list_b: len(set(list_a).intersection(set(list_b))) > 0)

class_switch_30_concat_info["50 Most frequent same as 30 Most frequent"] = (
    do_lists_overlap(
        class_switch_30_concat_info["50 Most frequent class"],
        class_switch_30_concat_info["30 Most frequent class"],
    )
)
class_switch_30_concat_info["80 Most frequent same as 30 Most frequent"] = (
    do_lists_overlap(
        class_switch_30_concat_info["80 Most frequent class"],
        class_switch_30_concat_info["30 Most frequent class"],
    )
)
class_switch_30_concat_info["50 Most frequent same as 30 DeepARG"] = is_class_in_list(
    class_switch_30_concat_info["30 DeepARG class"],
    class_switch_30_concat_info["50 Most frequent class"]
)
class_switch_30_concat_info["80 Most frequent same as 30 DeepARG"] = is_class_in_list(
    class_switch_30_concat_info["30 DeepARG class"],
    class_switch_30_concat_info["80 Most frequent class"]
)
class_switch_30_concat_info = class_switch_30_concat_info.reset_index()
class_switch_30_concat_summary = class_switch_30_concat_info.groupby(
    by=[
        "Model",
        "30 Most frequent class",
        "30 DeepARG class",
        "50 aligned",
        "50 predicted",
        "50 Most frequent same as 30 Most frequent",
        "50 Most frequent same as 30 DeepARG",
        "80 aligned",
        "80 predicted",
        "80 Most frequent same as 30 Most frequent",
        "80 Most frequent same as 30 DeepARG",
    ]
).count()[["Query"]]
class_switch_30_concat_summary["Percentage"] = class_switch_30_concat_summary[
    "Query"
] / (
    class_switch_30_concat_summary.index.get_level_values(0).map(
        class_switch_count.loc[30]
    )
)
class_switch_30_concat_summary.to_csv("class_switch_30_concat_summary.csv")

# Concat all queries that only class switched at 50
class_switch_50_concat_info = pd.concat(
    objs=[
        class_switch_50_info,
        pair_containing_30_class_switch_50_info,
        pair_containing_30_class_switch_50_pair_containing_80_info
    ]
)

class_switch_50_concat_info["30 Most frequent same as 50 Most frequent"] = (
    do_lists_overlap(
        class_switch_50_concat_info["30 Most frequent class"],
        class_switch_50_concat_info["50 Most frequent class"],
    )
)
class_switch_50_concat_info["80 Most frequent same as 50 Most frequent"] = (
    do_lists_overlap(
        class_switch_50_concat_info["80 Most frequent class"],
        class_switch_50_concat_info["50 Most frequent class"],
    )
)
class_switch_50_concat_info["30 Most frequent same as 50 DeepARG"] = is_class_in_list(
    class_switch_50_concat_info["50 DeepARG class"],
    class_switch_50_concat_info["30 Most frequent class"]
)
class_switch_50_concat_info["80 Most frequent same as 50 DeepARG"] = is_class_in_list(
    class_switch_50_concat_info["50 DeepARG class"],
    class_switch_50_concat_info["80 Most frequent class"]
)
class_switch_50_concat_info = class_switch_50_concat_info.reset_index()
class_switch_50_concat_summary = class_switch_50_concat_info.groupby(
    by=[
        "Model",
        "50 Most frequent class",
        "50 DeepARG class",
        "30 aligned",
        "30 predicted",
        "30 Most frequent same as 50 Most frequent",
        "30 Most frequent same as 50 DeepARG",
        "80 aligned",
        "80 predicted",
        "80 Most frequent same as 50 Most frequent",
        "80 Most frequent same as 50 DeepARG",
    ]
).count()[["Query"]]
class_switch_50_concat_summary["Percentage"] = class_switch_50_concat_summary[
    "Query"
] / (
    class_switch_50_concat_summary.index.get_level_values(0).map(
        class_switch_count.loc[50]
    )
)
class_switch_50_concat_summary.to_csv("class_switch_50_concat_summary.csv")

# Concat all queries that only class switched at 30 and 50
# and had the same class switches
class_switch_30_50_concat_info = pd.concat(
    objs=[
        same_class_switch_30_50_info,
        same_class_switch_30_50_pair_containing_80_info
    ]
)

class_switch_30_50_concat_info["80 Most frequent same as 30/50 Most frequent"] = (
    do_lists_overlap(
        class_switch_30_50_concat_info["80 Most frequent class"],
        class_switch_30_50_concat_info["50 Most frequent class"],
    )
)
class_switch_30_50_concat_info["80 Most frequent same as 30/50 DeepARG"] = (
    is_class_in_list(
        class_switch_30_50_concat_info["50 DeepARG class"],
        class_switch_30_50_concat_info["80 Most frequent class"],
    )
)
class_switch_30_50_concat_info = class_switch_30_50_concat_info.reset_index()
class_switch_30_50_concat_summary = class_switch_30_50_concat_info.groupby(
    by=[
        "Model",
        "30 Most frequent class",
        "30 DeepARG class",
        "80 aligned",
        "80 predicted",
        "80 Most frequent same as 30/50 Most frequent",
        "80 Most frequent same as 30/50 DeepARG",
    ]
).count()[["Query"]]
class_switch_30_50_concat_summary["Percentage of 30% class switches"] = (
    class_switch_30_50_concat_summary["Query"]
    / (
        class_switch_30_50_concat_summary.index.get_level_values(0).map(
            class_switch_count.loc[30]
        )
    )
)
class_switch_30_50_concat_summary["Percentage of 50% class switches"] = (
    class_switch_30_50_concat_summary["Query"]
    / (
        class_switch_30_50_concat_summary.index.get_level_values(0).map(
            class_switch_count.loc[50]
        )
    )
)
class_switch_30_50_concat_summary.to_csv("class_switch_30_50_concat_summary.csv")

# Concat all queries that only class switched at 30 and 50
# and had different class switches
not_same_class_switch_30_50_info["80 Most frequent same as 30 Most frequent"] = (
    do_lists_overlap(
        not_same_class_switch_30_50_info["80 Most frequent class"],
        not_same_class_switch_30_50_info["30 Most frequent class"],
    )
)
not_same_class_switch_30_50_info["80 Most frequent same as 50 Most frequent"] = (
    do_lists_overlap(
        not_same_class_switch_30_50_info["80 Most frequent class"],
        not_same_class_switch_30_50_info["50 Most frequent class"],
    )
)
not_same_class_switch_30_50_info["30 Most frequent same as 50 Most frequent"] = (
    do_lists_overlap(
        not_same_class_switch_30_50_info["30 Most frequent class"],
        not_same_class_switch_30_50_info["50 Most frequent class"],
    )
)

not_same_class_switch_30_50_info["80 Most frequent same as 30 DeepARG"] = (
    is_class_in_list(
        not_same_class_switch_30_50_info["30 DeepARG class"],
        not_same_class_switch_30_50_info["80 Most frequent class"],
    )
)
not_same_class_switch_30_50_info["80 Most frequent same as 50 DeepARG"] = (
    is_class_in_list(
        not_same_class_switch_30_50_info["50 DeepARG class"],
        not_same_class_switch_30_50_info["80 Most frequent class"],
    )
)
not_same_class_switch_30_50_info["50 Most frequent same as 30 DeepARG"] = (
    is_class_in_list(
        not_same_class_switch_30_50_info["30 DeepARG class"],
        not_same_class_switch_30_50_info["50 Most frequent class"],
    )
)
not_same_class_switch_30_50_info["30 Most frequent same as 50 DeepARG"] = (
    is_class_in_list(
        not_same_class_switch_30_50_info["50 DeepARG class"],
        not_same_class_switch_30_50_info["30 Most frequent class"],
    )
)
not_same_class_switch_30_50_info = not_same_class_switch_30_50_info.reset_index()
not_same_class_switch_30_50_summary = not_same_class_switch_30_50_info.groupby(
    by=[
        "Model",
        "30 Most frequent class",
        "30 DeepARG class",
        "50 Most frequent class",
        "50 DeepARG class",
        "80 aligned",
        "80 predicted",
        "80 Most frequent same as 30 Most frequent",
        "80 Most frequent same as 30 DeepARG",
        "80 Most frequent same as 50 Most frequent",
        "80 Most frequent same as 50 DeepARG",
        "30 Most frequent same as 50 Most frequent",
        "30 Most frequent same as 50 DeepARG",
        "50 Most frequent same as 30 DeepARG"
    ]
).count()[["Query"]]
not_same_class_switch_30_50_summary["Percentage of 30% class switches"] = (
    not_same_class_switch_30_50_summary["Query"]
    / (
        not_same_class_switch_30_50_summary.index.get_level_values(0).map(
            class_switch_count.loc[30]
        )
    )
)
not_same_class_switch_30_50_summary["Percentage of 50% class switches"] = (
    not_same_class_switch_30_50_summary["Query"]
    / (
        not_same_class_switch_30_50_summary.index.get_level_values(0).map(
            class_switch_count.loc[50]
        )
    )
)
not_same_class_switch_30_50_summary.to_csv("not_same_class_switch_30_50_summary.csv")

# Concat all queries that only class switched at 80
class_switch_80_concat_info = pd.concat(
    objs=[
        class_switch_80_info,
        pair_containing_30_50_class_switch_80_info,
        pair_containing_50_class_switch_80_info
    ]
)

class_switch_80_concat_info["30 Most frequent same as 80 Most frequent"] = (
    do_lists_overlap(
        class_switch_80_concat_info["30 Most frequent class"],
        class_switch_80_concat_info["80 Most frequent class"],
    )
)
class_switch_80_concat_info["50 Most frequent same as 80 Most frequent"] = (
    do_lists_overlap(
        class_switch_80_concat_info["50 Most frequent class"],
        class_switch_80_concat_info["80 Most frequent class"],
    )
)
class_switch_80_concat_info["30 Most frequent same as 80 DeepARG"] = is_class_in_list(
    class_switch_80_concat_info["80 DeepARG class"],
    class_switch_80_concat_info["30 Most frequent class"]
)
class_switch_80_concat_info["50 Most frequent same as 80 DeepARG"] = is_class_in_list(
    class_switch_80_concat_info["80 DeepARG class"],
    class_switch_80_concat_info["50 Most frequent class"]
)
class_switch_80_concat_info = class_switch_80_concat_info.reset_index()
class_switch_80_concat_summary = class_switch_80_concat_info.groupby(
    by=[
        "Model",
        "80 Most frequent class",
        "80 DeepARG class",
        "30 aligned",
        "30 predicted",
        "30 Most frequent same as 80 Most frequent",
        "30 Most frequent same as 80 DeepARG",
        "50 aligned",
        "50 predicted",
        "50 Most frequent same as 80 Most frequent",
        "50 Most frequent same as 80 DeepARG",
    ]
).count()[["Query"]]
class_switch_80_concat_summary["Percentage"] = class_switch_80_concat_summary[
    "Query"
] / (
    class_switch_80_concat_summary.index.get_level_values(0).map(
        class_switch_count.loc[80]
    )
)
class_switch_80_concat_summary.to_csv("class_switch_80_concat_summary.csv")
