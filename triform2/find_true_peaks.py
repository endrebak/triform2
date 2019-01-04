import numpy as np
import pandas as pd

from pyranges import PyRanges

import scipy.signal as ss

from triform2.find_peaks import pnorm

def ccf(x, y, lag_max = 100):

    result = ss.correlate(y - np.mean(y), x - np.mean(x), method='direct') / (np.std(y) * np.std(x) * len(y))
    length = (len(result) - 1) // 2
    lo = length - lag_max
    hi = length + (lag_max + 1)

    return result[lo:hi]




def _zscores(x, y, r=1):

    diff = (r * x) - y
    zs = diff/np.sqrt(r * (x + y))

    return zs


def compute_lags(peaks, cvg):


    def compute_lag(n, p):
        res = ccf(n, p)
        return res.argmax() - 100

    lags = []
    for k in peaks.chromosomes:
        df = peaks[k].df
        df = df["Start End".split()]

        p = cvg[k, "+"][df]
        n = cvg[k, "-"][df]

        n = [np.repeat(r.values, r.runs) for r in n]
        p = [np.repeat(r.values, r.runs) for r in p]

        for _n, _p in zip(n, p):
            lags.append(compute_lag(_p, _n))

    return np.array(lags)



def remove_overlapping_peaks(peaks):

    peak_types = set(peaks)

    if len(peak_types) == 3:

        peaks1, peaks2, peaks3 = [peaks[k] for k in sorted(peaks)]

        peaks2 = peaks2.no_overlap(peaks1)
        peaks3 = peaks3.no_overlap(peaks1)

        peaks3 = peaks3.no_overlap(peaks2)

        peaks = [peaks1, peaks2, peaks3]

    elif peak_types == set([1, 2]) or peak_types == set([1, 3]):

        peaks1, peaks_other = [peaks[k] for k in sorted(peaks)]
        peaks_other = peaks_other.no_ovelap(peaks1)

    elif peak_types == set([2, 3]):

        peaks = [peaks[k] for k in sorted(peaks)]

    else:

        peaks = [ list(peaks.values())[0] ]

    df = pd.concat([p.df for p in peaks])

    return df


def find_peaks(possible_peaks, cvg, args):

    flank_distance = args["flank_distance"]

    min_shift = args["min_shift"]

    _new_peaks = {}
    for peak_type, peaks in possible_peaks.items():

        peaks_f = peaks["+"].overlap(peaks["-"], strandedness=False)
        peaks_r = peaks["-"].overlap(peaks["+"], strandedness=False)

        # TODO: try False instead of first
        peaks_f_bool = np.concatenate(list(peaks_f.apply(lambda df, _: ~df.index.duplicated(keep="first"), as_pyranges=False).values()))
        peaks_r_bool = np.concatenate(list(peaks_r.apply(lambda df, _: ~df.index.duplicated(keep="first"), as_pyranges=False).values()))
        peaks_bool = peaks_r_bool & peaks_f_bool
        peaks_f.Keep = peaks_bool
        peaks_r.Keep = peaks_bool
        peaks_f = peaks_f.apply(lambda df, _: df[df.Keep].drop("Keep", 1))
        peaks_r = peaks_r.apply(lambda df, _: df[df.Keep].drop("Keep", 1))

        assert np.all(peaks_f.Chromosome.values == peaks_r.Chromosome.values)

        new_peaks = PyRanges(seqnames = peaks_f.Chromosome, starts = np.minimum(peaks_f.Start, peaks_r.Start), ends = np.maximum(peaks_f.End, peaks_r.End))

        if len(new_peaks) == 0:
            continue

        if peak_type == 1:
            new_peaks.Start -= flank_distance
            new_peaks.End += flank_distance
        elif peak_type == 2:
            new_peaks.Start -= flank_distance
        elif peak_type == 3:
            new_peaks.End += flank_distance
        else:
            assert 0

        lags = compute_lags(new_peaks, cvg)

        peaks_f.Lag = lags
        peaks_r.Lag = lags

        peaks_f = peaks_f.apply(lambda df, _: df[df.Lag > min_shift])
        peaks_r = peaks_r.apply(lambda df, _: df[df.Lag > min_shift])

        if len(peaks_f) == 0:
            continue

        new_locs = np.array(np.round((peaks_f.Location.values + peaks_r.Location.values)/2), dtype=np.long)
        _cvg = (peaks_f.CVG.values + peaks_r.CVG.values)
        surl = (peaks_f.SURL.values + peaks_r.SURL.values)
        surr = (peaks_f.SURR.values + peaks_r.SURR.values)

        if peak_type == 1:
            zs = _zscores(_cvg, surl + surr, 2)
            max_zs = _zscores(_cvg + surl + surr, 0, 2)
        elif peak_type == 2:
            zs = _zscores(_cvg, surl)
            max_zs = _zscores(_cvg + surl, 0)
        elif peak_type == 3:
            zs = _zscores(_cvg, surr)
            max_zs = _zscores(_cvg + surr, 0)

        peak_nlp = -pnorm(zs)/np.log(10)
        max_nlp = -pnorm(max_zs)/np.log(10)

        new_starts = np.minimum(peaks_f.Start, peaks_r.Start)
        new_ends = np.maximum(peaks_f.End, peaks_r.End)
        new_peaks = PyRanges(seqnames=peaks_f.Chromosome, starts=new_starts, ends=new_ends)
        new_peaks.Location = new_locs
        new_peaks.CVG = _cvg
        new_peaks.SURL = surl
        new_peaks.SURR = surr
        new_peaks.Form = peak_type
        new_peaks.NLP = peak_nlp
        new_peaks.MAX_NLP = max_nlp

        _new_peaks[peak_type] = new_peaks

    return remove_overlapping_peaks(_new_peaks)
