from pyrle.src.coverage import _coverage
from pyrle import Rle, PyRles
import numpy as np
import pandas as pd

cimport cython

from libc.stdint cimport uint32_t
from collections import defaultdict

from cython.operator import dereference, postincrement
from libcpp.algorithm cimport sort as stdsort
from libcpp.map cimport map as cppmap
from libcpp.algorithm cimport unique
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.map cimport map

cimport cpp_read_files as cr

import sys
import logging

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef class Vector32:

    cdef vector[uint32_t] wrapped_vector

    cdef push_back(self, uint32_t num):
        self.wrapped_vector.push_back(num)

    def sort(self):
        stdsort(self.wrapped_vector.begin(), self.wrapped_vector.end())

    def unique(self):
        self.wrapped_vector.erase(unique(self.wrapped_vector.begin(), self.wrapped_vector.end()), self.wrapped_vector.end())

    def __str__(self):
        return "[" + ", ".join([str(i) for i in self.wrapped_vector]) + "]"

    def __repr__(self):
        return str(self)

    def __len__(self):
        return self.wrapped_vector.size()

    def __iter__(self):
        # slow, only implemented to ease testing
        return (v for v in self.wrapped_vector)


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cpdef files_to_coverage(files, datatype, bool lenient, uint32_t read_width, bool drop_duplicates):

    cdef:
        Vector32 v
        Vector32 v2

        long[::1] positions
        double[::1] values

        long maxlength = -1

        long[::1] bin_arr
        bytes py_bytes
        char* c_string
        cr.genome_map cpp_tags
        map[cr.key, cr.intvec].iterator it
        uint32_t i = 0


    logging.info("Parsing {} file(s):".format(datatype))
    sys.stderr.flush()
    coverage = dict()
    sizes = dict()
    # not always used, just init for compiler
    _strict_coverage = dict()
    _strict_sizes = defaultdict(int)

    for f in files:

        tags = dict()
        logging.info("  " + f)
        sys.stderr.flush()


        py_bytes = f.encode()
        c_string = py_bytes


        if f.endswith(".bed"):
            cpp_tags = cr.read_bed(c_string, drop_duplicates, read_width)
        elif f.endswith(".bedpe"):
            cpp_tags = cr.read_bedpe(c_string, drop_duplicates, read_width)
        # elif f.endswith(".bam") or f.endswith(".sam"):
        #     cpp_tags = read_bam(f, drop_duplicates)
        elif f.endswith(".bed.gz"):
            cpp_tags = cr.read_bed_gz(c_string, drop_duplicates, read_width)
        elif f.endswith(".bedpe.gz"):
            cpp_tags = cr.read_bedpe_gz(c_string, drop_duplicates, read_width)

        it = cpp_tags.begin();

        # turn C++ based map into py dict
        while (it != cpp_tags.end()):
            chromosome = dereference(it).first.first.decode()

            strand = chr(dereference(it).first.second)

            v = Vector32()
            v.wrapped_vector = dereference(it).second
            tags[chromosome, strand] = v
            maxlength = max(maxlength, len(v))

            postincrement(it)

        positions_arr = np.zeros(maxlength, dtype=np.long)
        values_arr = np.ones(maxlength, dtype=np.double)
        values_arr[1::2] = -values_arr[1::2]
        positions = positions_arr
        values = values_arr

        for (chromosome, strand), v in tags.items():

            # print("1" * 50, chromosome, strand)
            sys.stdout.flush()
            for i in range(len(v)):
                positions[i] = v.wrapped_vector[i]

            # print("2" * 50, chromosome, strand, len(v))
            sys.stdout.flush()

            _mypos = np.copy(positions_arr[:len(v)])
            _myvals = np.copy(values_arr[:len(v)])

            s = pd.Series(index=_mypos, data=_myvals, name="Value").sort_index()
            # print(s)
            # print(s.index.values)

            cvg = _coverage(s.index.values, s.values)

            # print("here", chromosome, strand)
            # print(cvg[0])
            # print(cvg[1])
            sys.stdout.flush()

            cvg = Rle(cvg[0], cvg[1])
            print("cvg " * 50)
            print(cvg)


            if datatype == "input" or lenient:

                if (chromosome, strand) not in coverage:
                    coverage[chromosome, strand] = cvg
                    sizes[strand] = len(v) / 2
                else:
                    coverage[chromosome, strand] = cvg + coverage[chromosome, strand]
                    sizes[strand] += len(v) / 2

            else: # chip and not lenient

                _strict_sizes[strand] += len(v) / 2
                if (chromosome, strand) not in coverage:
                    _strict_coverage[chromosome, strand] = cvg
                else:
                    _strict_coverage[chromosome, strand] = cvg + _strict_coverage[chromosome, strand]

        if datatype == "chip" and not lenient:
            coverage[f] = PyRles(_strict_coverage)
            sizes[f] = _strict_sizes

            _strict_coverage = dict()
            _strict_sizes = defaultdict(int)


    if datatype == "input":
        coverage = PyRles(coverage)

    if datatype == "chip" and lenient:
        # so that strict/lenient effortlessly can be processed the same way
        coverage["mock_file"] = PyRles(coverage)

    return coverage, sizes
