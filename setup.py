
import os
import sys
from setuptools import setup, find_packages

from setuptools import find_packages, Extension, Command
from Cython.Build import cythonize

from subprocess import check_output
macros = []

# install_requires = ["coveralls"]
install_requires = ["pyranges", "pyrle"]

compile_options = ["-Ofast", "-Wall", "-std=c++11"] #, "-frename-registers", "-funroll-loops"] # , "-lgzstream", "-lz"

conda_path = check_output("which conda", shell=True).decode().strip()

conda_include = []
conda_lib = []
if conda_path:
    conda_base = conda_path.replace("bin/conda", "")
    conda_include.append(conda_base + "include/")
    conda_lib.append(conda_base + "lib/")

extensions = [Extension("triform2.src.files_to_coverage",
                        ["triform2/src/files_to_coverage.pyx", "triform2/src/gzstream.cpp"], language="c++",
                        include_dirs=conda_include,
                        library_dirs=conda_lib,
                        extra_compile_args=compile_options,
                        libraries=["z"])]

setup(
    name="triform2",
    packages=find_packages(),
    # package_dir=find_packages(),
    ext_modules = cythonize(extensions, annotate=True),
    scripts=["bin/triform2"],
    version="0.0.1",
    description=
    "Improved sensitivity, specificity and control of false discovery rates in ChIP-Seq peak finding.",
    author="Endre Bakken Stovner",
    author_email="endrebak85@gmail.com",
    url="http://github.com/endrebak/triform2",
    keywords=["ChIP-Seq"],
    license=["MIT"],
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Other Environment", "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"
    ],
    long_description=
    "Improved sensitivity, specificity and control of false discovery rates in ChIP-Seq peak finding.")
