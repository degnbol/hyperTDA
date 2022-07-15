#!/usr/bin/env python3
# USAGE: go to folder with .npy and it will print.
import os, sys
import numpy as np

def anynan(x, axis=0):
    return np.any(np.isnan(x), axis=axis)

def get_mat_dists(mat):
    return np.sqrt(np.sum(np.diff(mat, axis=0) ** 2, axis=1))

mats = [np.load(f) for f in os.listdir()]
mats = [mat[~anynan(mat, axis=1), :] for mat in mats]

mindists = [np.min(get_mat_dists(mat)) for mat in mats]

print(np.mean(mindists))
print(np.median(mindists))

