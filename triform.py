from itertools import product
from collections import defaultdict
from functools import reduce
import pandas as pd

from pyrle import PyRles, Rle

import numpy as np

c1 = "examples/srf_huds_Gm12878_rep1.bed"
c2 = "examples/srf_huds_Gm12878_rep2.bed"
b1 = "examples/backgr_huds_Gm12878_rep1.bed"
b2 = "examples/backgr_huds_Gm12878_rep2.bed"

chip_files = [c1]
input_files = [b1]

# chip_files = [c1, c2]
# input_files = [b1, b2]

names = "Chromosome Start End Strand".split()

usecols = [0, 1, 2, 5]

read_width = 100
flank_distance = 150
max_p = 0.1

def extend_df(df, read_width):

    df.Start += 1

    width = df.End - df.Start
    gaps = np.array(np.floor((read_width - width) / 2), dtype=np.long)
    df.Start = df.Start - gaps - 1
    df.End = df.Start + width + (2 * gaps)

    return df

def create_df(f, read_width):

    df1 = pd.read_table(f, sep="\t", usecols=usecols, names=names)
    df1 = extend_df(df1, read_width)

    return df1

def create_coverage(df):

    return PyRles(df, stranded=True)

# def get_sizes(f, df):

#     returnlen(df)

#     return sizes



def init_chip(coverage, flank_distance):

    """Creating left by chopping off rightmost fd, shifting remaining flank_distance to the right and vice-versa
    """

    replicates = list(coverage)

    outdict = {}
    for replicate, location in product(
            replicates, ["left", "right", "center"]):

        cvg = coverage[replicate]
        _cvg = {}
        for direction in ["+", "-"]:

            chromosomes = [k[0] for k in cvg.keys()]
            for c in chromosomes:

                try:
                    rle = cvg[c, direction]
                except:
                    continue

                if location == "left":
                    length = np.sum(rle.runs) - flank_distance
                    newrle = rle[:length]
                    newruns = np.concatenate([np.array([flank_distance], dtype=np.long), newrle.runs])
                    newvals = np.concatenate([np.array([0], dtype=np.double), newrle.values])
                    newrle = Rle(newruns, newvals)
                    _cvg[c, direction] = newrle
                elif location == "right":
                    length = np.sum(rle.runs)
                    newrle = Rle(rle.runs, rle.values)
                    newrle = newrle[flank_distance:]
                    # could be optimized by ending getitem after found start...
                    newruns = np.concatenate([newrle.runs, np.array([flank_distance], dtype=np.long)])
                    newvals = np.concatenate([newrle.values, np.array([0], dtype=np.double)])
                    newrle = Rle(newruns, newvals)
                    _cvg[c, direction] = newrle
                else:
                    _cvg[c, direction] = rle


        outdict[replicate, location] = PyRles(_cvg)


    return outdict


def init_background(coverage, flank_distance):

    """Creating left by chopping off rightmost fd, shifting remaining flank_distance to the right and vice-versa
    """

    replicates = list(coverage)

    outdict = {}
    for replicate, location in product(
            replicates, ["center"]):

        cvg = coverage[replicate]
        chromosomes = [k[0] for k in cvg.keys()]
        _cvg = {}
        for direction in ["+", "-"]:

            for c in chromosomes:

                try:
                    rle = cvg[c, direction]
                except:
                    continue

                rle = cvg[c, direction]
                _cvg[c, direction] = rle

        outdict[replicate, location] = PyRles(_cvg)


    return outdict


input_dfs = {f: create_df(f, read_width) for f in input_files}
input_coverage = {f: create_coverage(df) for f, df in input_dfs.items()}

# print(chip_coverage)



background = init_background(input_coverage, flank_distance)

# print(background)


def sum_background(background):


    first_file = list(background.keys())[0]
    first_key = list(background[first_file].keys())[0]
    background_sum = PyRles({first_key: Rle([1], [0])})

    for v in background.values():

        background_sum += v

    return background_sum

background_sum = sum_background(background)

# background_sum

# print(background)

# for k, v in background.items():
#     print(k)
#     print(v)

# def chromosome(background):

def get_locs(chip, where):

    results = {}
    for (f, loc), v in chip.items():

        if not loc == where:
            continue

        results[f] = v

    return results

def sum_chip(chip, where):


    locs = list(get_locs(chip, where).values())

    first_loc = locs[0]
    first_key = list(first_loc.keys())[0]
    chip_sum = PyRles({first_key: Rle([1], [0])})

    for v in locs:

        chip_sum += v

    return chip_sum


def sum_data(chip):

    center = sum_chip(chip, "center")
    left = sum_chip(chip, "left")
    right = sum_chip(chip, "right")

    return center, left, right



chip_dfs = {f: create_df(f, read_width) for f in chip_files}
chip_coverage = {f: create_coverage(df) for f, df in chip_dfs.items()}
chip = init_chip(chip_coverage, flank_distance)

cvg, left, right = sum_data(chip)
center = get_locs(chip, "center")

def qnorm(max_p):

    from scipy.stats import norm

    return norm.ppf(1 - max_p)

# def compute_peaks_and_zscores(cvg, center, left, right, chip, input, ratios, ratio, args):

def merge_runs(s):
    _new_rledict = {}
    for k, v in s.items():
        v.values[v.values < 0] = 0
        v.values[v.values > 0] = 1
        v = Rle(v.runs, v.values)
        _new_rledict[k] = v

    return PyRles(_new_rledict)

def compute_ok1(chip):

    sign_sum = None

    files = set([f[0] for f in chip.keys()])
    items = [[chip[f1, "center"], chip[f1, "left"], chip[f1, "right"]] for f1 in files]

    for c, l, r in items:
        s = merge_runs(((c * 2)  - l) - r)

        if sign_sum is None:
            sign_sum = s
        else:
            sign_sum *= s

    return sign_sum


def compute_ok23(chip, _type):
    # _type = "left" if _type == 2 else "right"
    sign_sum = None

    files = set([f[0] for f in chip.keys()])
    items = [[chip[f1, "center"], chip[f1, _type]] for f1 in files]

    for c, o in items:
        s = merge_runs(c - o)

        if sign_sum is None:
            sign_sum = s
        else:
            sign_sum *= s

    return sign_sum


def compute_ok4(ratios, center, background):

    # print(ratios.keys())
    # print(center.keys())
    # assert ratios.keys() == center.keys()

    sign_sum = {}
    for strand in ["+", "-"]:
        for (f, _strand), ratio in ratios.items():
            if _strand != strand:
                continue

            c = center[f]
            s = merge_runs((c[strand] * (1/ratio)) - background[strand])

            if not sign_sum.get(strand):
                sign_sum[strand] = s
            else:
                sign_sum[strand] *= s


    merge_strands = {}
    for d in sign_sum.values():
        for k, v in d.items():
            merge_strands[k] = v

    return PyRles(merge_strands)


def compute_peaks_and_zscores(cvg, center, left, right, chip, background, ratios, ratio, args):

    # center is the center coverage list for each chip_file
    # cvg is the summed coverage across files (same with right, left)

    min_z = qnorm(max_p)

    # left_right = left + right

    ok1 = compute_ok1(chip)
    # print(ok1["-"])

    ok2 = compute_ok23(chip, "left")
    # print("ok2" * 50)
    # print(ok2["-"])
    ok3 = compute_ok23(chip, "right")
    # print("ok3" * 50)
    # print(ok3["-"])
    ok4 = compute_ok4(ratios, center, background)
    print("ok4" * 50)
    print(ok4["-"])



def get_ratios(chip, background):

    # print("background", background)
    # print("chip", chip)

    # chip_sizes = {f: len(df) for f, df in chip_dfs.items()}
    chip_sizes = {}
    for f, df in chip_dfs.items():
        for s, sdf in df.groupby("Strand"):
            chip_sizes[f, s] = len(sdf)


    background_sizes = defaultdict(int) # sum(len(df) for df in background.values())
    for df in background.values():
        for s, sdf in df.groupby("Strand"):
            background_sizes[s] += len(sdf)

    per_file_ratios = {}
    ratios = {}
    if background is not None:
        for (f, strand), size in chip_sizes.items():
            # print("ratio compute" * 50)
            # print(f, len(df), background_sizes)
            per_file_ratios[f, strand] = size / background_sizes[strand]
        ratios[strand] = background_sizes[strand] / sum(size for (f, s), size in per_file_ratios.items() if strand == s)

    else:
        ratios = None
        ratios = {"+": 1, "-": 1}

    return per_file_ratios, ratios



ratios, ratio = get_ratios(chip_dfs, input_dfs)

args = {}
print(background_sum)
compute_peaks_and_zscores(cvg, center, left, right, chip, background_sum, ratios, ratio, args)
