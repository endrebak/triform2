from itertools import product

import pandas as pd

from collections import defaultdict

from pyrle import PyRles, Rle

import numpy as np

def extend_df(df, read_width):

    df.Start += 1

    width = df.End - df.Start
    gaps = np.array(np.floor((read_width - width) / 2), dtype=np.long)
    df.Start = df.Start - gaps - 1
    df.End = df.Start + width + (2 * gaps)

    return df

def create_df(f, read_width):

    usecols = [0, 1, 2, 5]
    names = "Chromosome Start End Strand".split()
    df = pd.read_table(f, sep="\t", usecols=usecols, names=names)
    df = extend_df(df, read_width)

    return df

def create_coverage(df):

    return PyRles(df, stranded=True)



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




def sum_background(background):


    first_file = list(background.keys())[0]
    first_key = list(background[first_file].keys())[0]
    background_sum = PyRles({first_key: Rle([1], [0])})

    for v in background.values():

        background_sum += v

    return background_sum




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



def get_ratios(chip_sizes, background_sizes):


    per_file_ratios = {}
    ratios = defaultdict(int)
    if background_sizes is not None:
        for (f, strand), size in chip_sizes.items():
            per_file_ratios[f, strand] = size / background_sizes[strand]
            ratios[strand] += size

        ratios = {s: background_sizes[s]/v for s, v in ratios.items()}

    else:
        per_file_ratios = None
        ratios = {"+": 1, "-": 1}

    return per_file_ratios, ratios


def init_files(treatment_files, input_files, args):

    # do chromosome by chromosome?

    flank_distance = args["flank_distance"]
    read_width = args["read_width"]

    chip_dfs = {f: create_df(f, read_width) for f in treatment_files}

    chip_sizes = {}
    for f, df in chip_dfs.items():
        for s, sdf in df.groupby("Strand"):
            chip_sizes[f, s] = len(sdf)

    chip_coverage = {f: create_coverage(df) for f, df in chip_dfs.items()}
    # del defs
    chip = init_chip(chip_coverage, flank_distance)
    # del chip_coverage

    cvg, left, right = sum_data(chip)
    center = get_locs(chip, "center")

    if input_files:
        input_dfs = {f: create_df(f, read_width) for f in input_files}

        background_sizes = defaultdict(int)

        for df in input_dfs.values():
            for s, sdf in df.groupby("Strand"):
                background_sizes[s] += len(sdf)

        input_coverage = {f: create_coverage(df) for f, df in input_dfs.items()}
        # del input_dfs

        background = init_background(input_coverage, flank_distance)
        background_sum = sum_background(background)
        # del background
    else:
        background_sizes = None
        background_sum = None

    ratios, ratio = get_ratios(chip_sizes, background_sizes)

    return cvg, center, left, right, chip, background_sum, ratios, ratio


def cython_init_files(treatment_files, input_files, args):

    from triform2.src.files_to_coverage import files_to_coverage

    chip, chip_sizes = files_to_coverage(treatment_files, "chip", args["lenient"], args["read_width"], args["drop_duplicates"])
    chip = init_chip(chip, args["flank_distance"])

    cvg, left, right = sum_data(chip)
    center = get_locs(chip, "center")

    if input_files:
        background_sum, background_sizes = files_to_coverage(input_files, "input", args["lenient"], args["read_width"], args["drop_duplicates"])

    ratios, ratio = get_ratios(chip_sizes, background_sizes)

    return cvg, center, left, right, chip, background_sum, ratios, ratio
