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
min_width = 10
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
            s = merge_runs((c[strand] * (ratio)) - background[strand])

            if not sign_sum.get(strand):
                sign_sum[strand] = s
            else:
                sign_sum[strand] *= s


    merge_strands = {}
    for d in sign_sum.values():
        for k, v in d.items():
            merge_strands[k] = v

    return PyRles(merge_strands)


def _zscores(x, y, ratio=1):

    _zscores = r("""
function(x,y,r) {  # r = size.y/size.x
  dif <- (r*x-y)
  zs <- dif/sqrt(r*(x+y))
  zs[!dif] <- 0
  zs
}
""")
    return _zscores(x, y, ratio)


def zscores(x, y, ratio=1):

    _ratio = ratio

    # zs = ratio * (x - y)
    new_pyrle = {}
    print("ratio")
    print(ratio)
    for k, v in (x + y).items():

        if isinstance(ratio, dict):
            _ratio = ratio[k[1]]

        _ratio = _ratio
        difference = (_ratio * x[k]) - y[k]

        denominator = Rle(v.runs, np.nan_to_num(np.sqrt((_ratio * v).values)))
        zs = difference / denominator
        zs = zs.numbers_only()
        zs.values[zs.values < 0] = 0
        zs = zs.defragment()
        new_pyrle[k] = zs


    # denominator = PyRles(denominator)
    # print("denominator")
    # print(denominator)
    # zs = difference / denominator

    return PyRles(new_pyrle)

def slice_min_z(pyrle, min_z):

    new_pyrle = {}
    for k, v in pyrle.items():
        # print("k " * 50, k)
        v = v.copy()
        # print(v)
        v.values[v.values < min_z] = 0
        # print(v)
        new_pyrle[k] = v.defragment()
        # print(new_pyrle[k])

    return PyRles(new_pyrle)


def remove_too_short(pyrle, min_width):

    new_pyrle = {}
    for k, v in pyrle.items():
        v = v.copy()
        v.values[v.values != 0] = 1
        v = v.defragment()
        v.values[v.runs <= min_width] = 0
        new_pyrle[k] = v.defragment()

    return PyRles(new_pyrle)



def remove_empty(df, _):
    df = df[~(df.Score == 0)]

    return df

def add_1_to_start(df, _):

    df.Start += 1
    return df

# result = ranges1.apply(lambda df, _: df[~(df.Score == 0)]).cluster(strand=True)

# result = ranges1.apply(remove_empty).cluster(strand=True).apply(add_1_to_start)

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
    # print("ok4" * 50)
    # print(ok4["-"])

    # print("")
    # print(" cvg " * 50)
    # print(cvg)
    # print(left + right)
    zs1 = (ok1 * zscores(cvg, left + right, 2)).defragment(numbers_only=True)
    zs2 = (ok2 * zscores(cvg, left)).defragment(numbers_only=True)
    zs3 = (ok3 * zscores(cvg, right)).defragment(numbers_only=True)
    zs4 = (ok4 * zscores(cvg, background, ratio)).defragment(numbers_only=True)

    # print("zs1")
    # print(zs1["-"])
    peaks1 = slice_min_z(zs1, min_z)
    peaks2 = slice_min_z(zs2, min_z)
    peaks3 = slice_min_z(zs3, min_z)
    peaks4 = slice_min_z(zs4, min_z)

    subset1 = remove_too_short(peaks1, min_width)
    subset2 = remove_too_short(peaks2, min_width)
    subset3 = remove_too_short(peaks3, min_width)
    subset4 = remove_too_short(peaks4, min_width)

    # print("peaks1")
    # print(peaks1["-"])
    # print("subset1")
    # print(subset1["-"])
    peaks1 *= subset1
    peaks2 *= subset2
    peaks3 *= subset3
    peaks4 *= subset4

    peaks1 *= peaks4 # * subset1
    peaks2 *= peaks4 # * subset2
    peaks3 *= peaks4 # * subset3
    # peaks4 *= subset4
    # print("peaks1")
    # print(peaks1["-"])

    subset1 = remove_too_short(peaks1, min_width)
    subset2 = remove_too_short(peaks2, min_width)
    subset3 = remove_too_short(peaks3, min_width)

    peaks1 *= subset1
    peaks2 *= subset2
    peaks3 *= subset3

    # only calling remove_too_short here to set values to 1
    peaks1 = remove_too_short(peaks1, min_width)
    peaks2 = remove_too_short(peaks2, min_width)
    peaks3 = remove_too_short(peaks3, min_width)

    # print("peaks1")
    # print(peaks1["-"])
    # print("peaks2")
    # print(peaks2["-"])
    # print("peaks3")
    # print(peaks3["-"])

    _peaks = [peaks1, peaks2, peaks3]
    _peaks = [p.to_ranges().apply(remove_empty).cluster(strand=True).apply(add_1_to_start) for p in _peaks]
    _zscores = [zs1, zs2, zs3]

    return _peaks, _zscores


# ranges1 = peaks1.to_ranges()

def remove_empty(df, _):
    df = df[~(df.Score == 0)]

    return df

def add_1_to_start(df, _):

    df.Start += 1
    return df

# result = ranges1.apply(lambda df, _: df[~(df.Score == 0)]).cluster(strand=True)



    # return _peaks, _zscores

    # print("peaks4 " * 50)
    # print(peaks4["-"])

    # runs = peaks2["chrY", "-"].runs
    # values = peaks2["chrY", "-"].values
    # import pandas as pd
    # pd.DataFrame(runs, values).to_csv("chry_m.txt", sep=" ")


    # print("zs1" * 50)
    # print(zs1["-"])
    # print("zs2" * 50)
    # print(zs2["-"])
    # print("zs3" * 50)
    # print(zs3["-"])
    # print("zs4" * 50)
    # print(zs4["-"])




def get_ratios(chip, background):

    chip_sizes = {}
    for f, df in chip_dfs.items():
        for s, sdf in df.groupby("Strand"):
            chip_sizes[f, s] = len(sdf)


    background_sizes = defaultdict(int)
    for df in background.values():
        for s, sdf in df.groupby("Strand"):
            background_sizes[s] += len(sdf)

    per_file_ratios = {}
    ratios = defaultdict(int)
    if background is not None:
        for (f, strand), size in chip_sizes.items():
            per_file_ratios[f, strand] = size / background_sizes[strand]
            ratios[strand] += size

        ratios = {s: background_sizes[s]/v for s, v in ratios.items()}

    else:
        per_file_ratios = None
        ratios = {"+": 1, "-": 1}

    return per_file_ratios, ratios

ratios, ratio = get_ratios(chip_dfs, input_dfs)

args = {}
# print(background_sum)
all_peaks, zs = compute_peaks_and_zscores(cvg, center, left, right, chip, background_sum, ratios, ratio, args)

for peak_type, peaks in enumerate(all_peaks, 1):
    print(peak_type)
    print(peaks)
    # peaks1 = peaks[0]

# print("as coverage" * 50)
# print(peaks1)
# rle1 = peaks1["chrY", "+"]


# import pandas as pd
# pd.DataFrame(rle1.runs, rle1.values).to_csv("chry_f.txt", sep=" ")
# 431 + 434

# print(result)
# result.df.to_csv("pr1_2.txt", sep=" ")
# print(result)
# print(zs[0])
