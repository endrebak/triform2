from itertools import product

import logging

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


        np.seterr(all="ignore")
        denominator = Rle(v.runs, np.nan_to_num(np.sqrt((_ratio * v).values)))
        zs = difference / denominator
        np.seterr(all="warn")
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


def set_all_nonzero_to_one(pyrle):

    new_pyrle = {}
    for k, v in pyrle.items():
        v = v.copy()
        v.values[v.values != 0] = 1
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

    logging.info("ok1")
    ok1 = compute_ok1(chip)
    logging.info("ok2")
    ok2 = compute_ok23(chip, "left")
    logging.info("ok3")
    ok3 = compute_ok23(chip, "right")
    logging.info("ok4")
    ok4 = compute_ok4(ratios, center, background)

    logging.info("zs1")
    zs1 = (ok1 * zscores(cvg, left + right, 2)).defragment(numbers_only=True)
    logging.info("zs2")
    zs2 = (ok2 * zscores(cvg, left)).defragment(numbers_only=True)
    logging.info("zs3")
    zs3 = (ok3 * zscores(cvg, right)).defragment(numbers_only=True)
    logging.info("zs4")
    zs4 = (ok4 * zscores(cvg, background, ratio)).defragment(numbers_only=True)

    logging.info("slice min")
    peaks1 = set_all_nonzero_to_one(slice_min_z(zs1, min_z)).to_ranges()
    peaks2 = set_all_nonzero_to_one(slice_min_z(zs2, min_z)).to_ranges()
    peaks3 = set_all_nonzero_to_one(slice_min_z(zs3, min_z)).to_ranges()
    peaks4 = set_all_nonzero_to_one(slice_min_z(zs4, min_z)).to_ranges()
    logging.info("without zero?")
    print(peaks4)


    # logging.info("remove_too_short")
    # peaks1 = peaks1.apply(lambda df, _: df[df.Score > 0])
    # peaks2 = peaks2.apply(lambda df, _: df[df.Score > 0])
    # peaks3 = peaks3.apply(lambda df, _: df[df.Score > 0])
    # peaks4 = peaks4.apply(lambda df, _: df[df.Score > 0])


    logging.info("multiply subset")
    peaks1 = peaks1.intersect(peaks4)
    peaks2 = peaks2.intersect(peaks4)
    peaks3 = peaks3.intersect(peaks4)

    # logging.info("remove_too_short")
    peaks1 = peaks1.apply(lambda df, _: df[(df.End - df.Start) > min_width])
    peaks2 = peaks2.apply(lambda df, _: df[(df.End - df.Start) > min_width])
    peaks3 = peaks3.apply(lambda df, _: df[(df.End - df.Start) > min_width])

    _peaks = [peaks1, peaks2, peaks3]

    # print(len(zs1))
    zs1 = zs1[peaks1]
    zs2 = zs2[peaks2]
    zs3 = zs3[peaks3]
    # print(len(zs1))
    # print(zs1)
    # raise
    # print(peaks1)

    _zscores = [zs1, zs2, zs3]
    # for peaks in _peaks:
        # print("ooooo")
        # print(peaks.apply(lambda df, _: df[df.index.duplicated()]))

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

    print("peaks and zscores")
    all_peaks, zs = _compute_peaks_and_zscores(cvg, center, left, right, chip, background_sum, ratios, ratio, args)
    print("peaks and zscores done")

    min_er = args["min_enrichment"]

    peaks_with_info = {}
    for peak_type, peaks in enumerate(all_peaks, 1):

        # print("find max start")
        # print(list(len(v) for v in zs[peak_type - 1].values()))
        # print(peaks)
        # print(zs[peak_type - 1].values())
        # t1 = list(zs[peak_type - 1].values())[0]
        # print(t1)
        # print(max(t1[1]))
        max_zs = {}
        for k, v in zs[peak_type - 1].items():
            max_zs[k] = np.array([max(v2[1]) for v2 in v])

        # max_zs = np.array(max_zs)
        # print("find max end")

        # print("len max_zs:", sum(len(v) for v in max_zs.values()))

        result = {k: -(pnorm(v)/np.log(10)) for k, v in max_zs.items()}
        # print(len(peaks))
        # print(len(np.concatenate([result[k] for k in natsorted(result)])))
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

        # print(len(vals_f))
        vals_f = vals_f[vals_f > 1]
        vals_r = vals_r[vals_r > 1]

        if peak_type == 1:
            min_er_f = np.percentile(vals_f, min_er * 100)
            min_er_r = np.percentile(vals_r, min_er * 100)

        vals_f = vals_f > min_er_f
        vals_r = vals_r > min_er_r

        # print(np.sum(vals_f))
        # print(len(vals_f))
        # print(peaks["+"])

        peaks["+"].Enrichment = vals_f
        peaks["-"].Enrichment = vals_r

        peaks_loc["+"].Enrichment = vals_f
        peaks_loc["-"].Enrichment = vals_r

        peaks = peaks.apply(lambda df, _: df[df.Enrichment].drop("Enrichment", axis=1))
        peaks_loc = peaks_loc.apply(lambda df, _: df[df.Enrichment].drop("Enrichment", axis=1))
        peaks_loc.Start += 1
        peaks_loc.End += 1

        chip_cvg = np.array(np.concatenate([cvg[k][peaks[k].Location] for k in cvg.keys() if not peaks[k].empty()]), dtype=np.long)
        left_cvg = np.array(np.concatenate([left[k][peaks[k].Location] for k in left.keys() if not peaks[k].empty()]), dtype=np.long)
        right_cvg = np.array(np.concatenate([right[k][peaks[k].Location] for k in right.keys() if not peaks[k].empty()]), dtype=np.long)

        peaks.CVG = chip_cvg
        peaks.SURL = left_cvg
        peaks.SURR = right_cvg

        peaks.drop_empty()

        peaks_with_info[peak_type] = peaks

    return peaks_with_info
