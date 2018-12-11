import pandas as pd

from pyrle import PyRles

import numpy as np

f1 = "examples/srf_huds_Gm12878_rep1.bed"
f2 = "examples/backgr_huds_Gm12878_rep1.bed"

names = "Chromosome Start End Strand".split()

usecols = [0, 1, 2, 5]

def extend_df(df, read_width):

    df.Start += 1

    width = df.End - df.Start
    gaps = np.array(np.floor((read_width - width) / 2), dtype=np.long)
    df.Start = df.Start - gaps - 1
    df.End = df.Start + width + (2 * gaps)

    return df




df1 = pd.read_table(f1, sep="\t", usecols=usecols, names=names)
df1_prime = extend_df(df1, 100)
rle1 = PyRles(df1_prime, stranded=True)
print(rle1)
df2 = pd.read_table(f2, sep="\t", usecols=usecols, names=names)
df2_prime = extend_df(df2, 100)
rle2 = PyRles(df2_prime, stranded=True)
print(rle2)
