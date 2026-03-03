from Bio.SeqIO import FastaIO
import pandas as pd

def get_from_nuc(index, from_codon):
    name = index[0]
    frame = int(index[1].split('_')[0])
    direction = index[1].split('_')[1]
    if direction == 'FOR':
        return from_codon * 3 + frame - 1
    else:
        return orf_windows.at[name, 'length'] - (from_codon * 3 + frame - 1)
    
def get_to_nuc(index, to_codon):
    name = index[0]
    frame = int(index[1].split('_')[0])
    direction = index[1].split('_')[1]
    if direction == 'FOR':
        return to_codon * 3 + frame + 2
    else:
        return (orf_windows.at[name, 'length'] - frame + 1) - to_codon * 3

orf_windows = pd.read_csv(
    filepath_or_buffer="best_windows.tsv",
    delimiter='\t',
    header=0,
    index_col=0)

orf_sequences = {
    rec.id: len(rec)
    for rec in FastaIO.FastaIterator(source="clean_ORFs.fa")}

orf_windows["length"] = orf_windows.apply(
    lambda x: orf_sequences[x.name], axis=1)

all_hitdata: list[dict] = list()
for direction in ["FOR", "REV"]:
    for frame in range(1,4):
        hitdata = pd.read_csv(
            filepath_or_buffer=f"frame_{frame}_{direction}_hitdata.tsv",
            delimiter='\t',
            header=0)
        hitdata['Query'] = hitdata['Query'].apply(
            lambda x: x.split(' ')[2][1:-2])
        hitdata['Frame'] = f"{frame}_{direction}"
        all_hitdata.extend(hitdata.to_dict(orient='records'))

all_hitdata_df = pd.DataFrame(data=all_hitdata)
all_hitdata_df = all_hitdata_df.loc[all_hitdata_df['Query'].isin(orf_windows.index.to_list())]
all_hitdata_df.set_index(
    keys=['Query', 'Frame'],
    inplace=True)

all_hitdata_df.drop(
    columns=["Hit type","PSSM-ID","E-Value","Bitscore","Incomplete"],
    inplace=True)

all_hitdata_df["From Nuc"] = all_hitdata_df.apply(
    lambda x: get_from_nuc(x.name, x['From']), axis=1)
all_hitdata_df["To Nuc"] = all_hitdata_df.apply(
    lambda x: get_to_nuc(x.name, x['To']), axis=1)
all_hitdata_df["Window Start"] = all_hitdata_df.apply(
    lambda x: orf_windows.at[x.name[0], 'window_start'], axis=1)
all_hitdata_df["Window End"] = all_hitdata_df.apply(
    lambda x: orf_windows.at[x.name[0], 'window_end'], axis=1)

all_hitdata_df.to_csv("test.csv")