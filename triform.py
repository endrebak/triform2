from itertools import product
import pandas as pd

from pyrle import PyRles, Rle

import numpy as np

f1 = "examples/srf_huds_Gm12878_rep1.bed"
f2 = "examples/backgr_huds_Gm12878_rep1.bed"

chip_files = [f1]
input_files = [f2]
names = "Chromosome Start End Strand".split()

usecols = [0, 1, 2, 5]

read_width = 100
flank_distance = 150

def extend_df(df, read_width):

    df.Start += 1

    width = df.End - df.Start
    gaps = np.array(np.floor((read_width - width) / 2), dtype=np.long)
    df.Start = df.Start - gaps - 1
    df.End = df.Start + width + (2 * gaps)

    return df

def create_coverage(f, read_width):

    df1 = pd.read_table(f1, sep="\t", usecols=usecols, names=names)
    df1_prime = extend_df(df1, 100)
    return PyRles(df1_prime, stranded=True)


def init_chip(coverage, flank_distance):

    """Creating left by chopping off rightmost fd, shifting remaining flank_distance to the right and vice-versa
    """

    replicates = list(coverage)

    outdict = {}
    for replicate, direction, location in product(
            replicates, ["+", "-"], ["left", "right", "center"]):

        cvg = coverage[replicate]

        chromosomes = [k[0] for k in cvg.keys()]
        for c in chromosomes:
            rle = cvg[c, direction]

            if location == "left":
                length = np.sum(rle.runs) - flank_distance
                newrle = rle[:length]
                newruns = np.concatenate([np.array([flank_distance], dtype=np.long), newrle.runs])
                newvals = np.concatenate([np.array([0], dtype=np.double), newrle.values])
                newrle = Rle(newruns, newvals)
                outdict[c, direction, replicate, location] = newrle
            elif location == "right":
                length = np.sum(rle.runs)
                newrle = Rle(rle.runs, rle.values)
                newrle = newrle[flank_distance:]
                # could be optimized by ending getitem after found start...
                newruns = np.concatenate([newrle.runs, np.array([flank_distance], dtype=np.long)])
                newvals = np.concatenate([newrle.values, np.array([0], dtype=np.double)])
                newrle = Rle(newruns, newvals)
                outdict[c, direction, replicate, location] = newrle
            else:
                outdict[c, direction, replicate, location] = rle


    return outdict


def init_background(coverage, flank_distance):

    """Creating left by chopping off rightmost fd, shifting remaining flank_distance to the right and vice-versa
    """

    replicates = list(coverage)

    outdict = {}
    for replicate, direction, location in product(
            replicates, ["+", "-"], ["center"]):

        cvg = coverage[replicate]

        chromosomes = [k[0] for k in cvg.keys()]
        for c in chromosomes:
            rle = cvg[c, direction]

            outdict[c, direction, replicate, location] = rle


    return outdict


chip_coverage = {f: create_coverage(f, read_width) for f in chip_files}
input_coverage = {f: create_coverage(f, read_width) for f in input_files}

# print(chip_coverage)

# chip = init_chip(chip_coverage, flank_distance)
# for k, v in chip.items():
#     print(k)
#     print(v)


background = init_background(input_coverage, flank_distance)
for k, v in background.items():
    print(k)
    print(v)
