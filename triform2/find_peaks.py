from itertools import product

import pandas as pd

from collections import defaultdict

from pyrle import PyRles, Rle

import numpy as np

from scipy.stats import norm

from pyranges import PyRanges

from natsort import natsorted


def qnorm(max_p):

    return norm.ppf(1 - max_p)


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



def zscores(x, y, ratio=1):

    _ratio = ratio

    new_pyrle = {}
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

    return PyRles(new_pyrle)

def slice_min_z(pyrle, min_z):

    new_pyrle = {}
    for k, v in pyrle.items():
        v = v.copy()
        v.values[v.values < min_z] = 0
        new_pyrle[k] = v.defragment()

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

def _compute_peaks_and_zscores(cvg, center, left, right, chip, background, ratios, ratio, args):

    # center is the center coverage list for each chip_file
    # cvg is the summed coverage across files (same with right, left)

    min_z = qnorm(args["max_p"])
    min_width = args["min_width"]

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
    # _3 = zs3["chrY", "-"].values
    # _3 = _3[_3 != 0]
    # print(_3)
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
    # print(zs1)
    zs1 *= peaks1
    # print(zs1)
    zs2 *= peaks2
    # _3 = zs3["chrY", "-"].values
    # _3 = _3[_3 != 0]
    # print(pd.Series(_3))
    # print(zs3)
    zs3 *= peaks3
    # print(zs3)
    # print(zs3["chrY", "-"].values)
    # _3 = zs3["chrY", "-"].values
    # _3 = _3[_3 != 0]
    # print(pd.Series(_3))
    # print("peaks1")
    # print(peaks1["-"])
    # print("peaks2")
    # print(peaks2["-"])
    # print("peaks3")
    # print(peaks3["-"])# #

    _peaks = [peaks1, peaks2, peaks3]
    _peaks = [p.to_ranges().apply(remove_empty).cluster(strand=True).apply(add_1_to_start) for p in _peaks]
    _zscores = [zs1, zs2, zs3]

    return _peaks, _zscores


def find_max(zscores):

    # TODO: cythonize me?

    max_zs_dict = defaultdict(list)
    for k, v in zscores.items():
        max_zs = []
        max_z = 0

        first_nonzero = 0
        for i in v.values:
            if i == 0:
                first_nonzero += 1
            else:
                break

        for i in v.values[first_nonzero:]:
            if i == 0:
                max_zs.append(max_z)
                max_z = 0
            else:
                max_z = max(max_z, i)

        if max_z != 0:
            max_zs.append(max_z)

        max_zs_dict[k] = np.array(max_zs)

    return max_zs_dict


def pnorm(max_z):

    np.seterr(divide="ignore")
    r = np.log(1 - norm.cdf(max_z))
    np.seterr(divide="warn")

    return r



def compute_peaks_and_zscores(cvg, center, left, right, chip, background_sum, ratios, ratio, args):

    all_peaks, zs = _compute_peaks_and_zscores(cvg, center, left, right, chip, background_sum, ratios, ratio, args)

    min_er = args["min_enrichment"]

    peaks_with_info = {}
    for peak_type, peaks in enumerate(all_peaks, 1):

        max_zs = find_max(zs[peak_type - 1])

        result = {k: -(pnorm(v)/np.log(10)) for k, v in max_zs.items()}
        peaks.NLP = np.around(np.concatenate([result[k] for k in natsorted(result)]), 3)

        peaks.Location = np.array(np.ceil((peaks.Start + peaks.End)/2), dtype=np.long)

        peaks.Type = peak_type

        peaks_loc = PyRanges(seqnames=peaks.Chromosome, starts=peaks.Location, ends=peaks.Location + 1, strands=peaks.Strand)
        loc_cvg = peaks_loc.coverage()

        chip_cvg = loc_cvg * cvg
        bg_cvg = loc_cvg * background_sum

        peak_enrich_cvg_f = 1 + (ratio["+"] * chip_cvg["+"])
        peak_enrich_cvg_r = 1 + (ratio["-"] * chip_cvg["-"])
        peak_enrich_cvg = PyRles({k: v for k, v in list(peak_enrich_cvg_r.items() + peak_enrich_cvg_f.items())})

        peak_enrich_ref = 1 + (bg_cvg)
        peak_enrich = peak_enrich_cvg / peak_enrich_ref

        vals_f = np.concatenate([peak_enrich[k].values for k in peak_enrich["+"].keys()])
        vals_r = np.concatenate([peak_enrich[k].values for k in peak_enrich["-"].keys()])
        vals_f = vals_f[np.isfinite(vals_f)]
        vals_r = vals_r[np.isfinite(vals_r)]

        vals_f = vals_f[vals_f > 1]
        vals_r = vals_r[vals_r > 1]

        if peak_type == 1:
            min_er_f = np.percentile(vals_f, min_er * 100)
            min_er_r = np.percentile(vals_r, min_er * 100)

        vals_f = vals_f > min_er_f
        vals_r = vals_r > min_er_r

        peaks["+"].Enrichment = vals_f
        peaks["-"].Enrichment = vals_r

        peaks_loc["+"].Enrichment = vals_f
        peaks_loc["-"].Enrichment = vals_r

        peaks = peaks.apply(lambda df, _: df[df.Enrichment].drop("Enrichment", axis=1))
        peaks_loc = peaks_loc.apply(lambda df, _: df[df.Enrichment].drop("Enrichment", axis=1))
        peaks_loc.Start += 1
        peaks_loc.End += 1

        chip_cvg = np.array(np.concatenate([cvg[k][peaks[k].Location] for k in cvg.keys()]), dtype=np.long)
        left_cvg = np.array(np.concatenate([left[k][peaks[k].Location] for k in left.keys()]), dtype=np.long)
        right_cvg = np.array(np.concatenate([right[k][peaks[k].Location] for k in right.keys()]), dtype=np.long)

        peaks.CVG = chip_cvg
        peaks.SURL = left_cvg
        peaks.SURR = right_cvg

        peaks_with_info[peak_type] = peaks

    return peaks_with_info
