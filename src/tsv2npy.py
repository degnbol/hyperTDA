#!/usr/bin/env python3
# USAGE: ./tsv2npy.py FILENAME*.tsv
# writes to FILENAME.npy
import numpy as np
import sys
from os.path import splitext

FILENAMES = sys.argv[1:]

for FILENAME in FILENAMES:
    name = splitext(FILENAME)[0]
    np.save(name + ".npy", np.loadtxt(FILENAME))

