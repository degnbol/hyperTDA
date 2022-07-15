#!/usr/bin/env python3
# coding: utf-8
import numpy as np
from sys import argv
from os.path import basename, splitext
import glob
import pandas as pd

def anynan(x, axis=0):
    return np.any(np.isnan(x), axis=axis)

def nNaN(curve):
    """
    Get number of points in a curve that contain NaN.
    """
    return np.sum(anynan(curve, axis=1))

def interpolate_nan(curve):
    """
    Replace points containing NaN with linear interpolation of neighbors.
    curve: matrix with xyz
    """
    if anynan(curve[0]) or anynan(curve[-1]):
        raise NotImplementedError
    
    curve_interp = []
    i = 0
    while i < len(curve):
        if not anynan(curve[i]):
            curve_interp.append(curve[i])
            i += 1
        else:
            # n == number of nans in a row to convert.
            n = 1
            while anynan(curve[i+n]): n += 1
            for j in range(1, n+1):
                curve_interp.append((1-j/(n+1))*curve[i-1] + (j/(n+1))*curve[i+n])
            i += n
    
    return np.asarray(curve_interp)

def interpolate(curve, dist=0.02):
    """
    Add interpolated points between xyzs of curve.
    dist: place points approximately this far apart
    """
    curve_interp = []
    interp_idx = []
    
    for i in range(len(curve)-1):
        # add first real point of a segment
        curve_interp.append(curve[i])
        interp_idx.append(False)
        
        d = np.linalg.norm(curve[i] - curve[i+1])
        # number of interpolation points propertional to distance between two real points.
        n = int(np.ceil(d / dist))
        # starting from 1 means we don't add the first real point twice
        for j in range(1, n):
            curve_interp.append((1-j/n)*curve[i] + (j/n)*curve[i+1])
            interp_idx.append(True)
        
    # add terminal real point
    curve_interp.append(curve[-1])
    interp_idx.append(False)
    
    return np.asarray(curve_interp), np.asarray(interp_idx)

def main():   
    # USAGE: ./interpolate.py INDIR/ OUTDIR/ DIST MAXNAN
    # DIST: approx distance between interpolated points
    # MAXNAN: number of points with NaNs allowed.
    # Outfiles will be given the same name as infiles but placed in OUTDIR
    INDIR, OUTDIR, DIST, MAXNAN = argv[1], argv[2], float(argv[3]), int(argv[4])

    for fname in glob.glob(f"{INDIR}/*.npy"):
        curve = np.load(fname)
        name = splitext(basename(fname))[0]
        
        if nNaN(curve) > 0:
            if nNaN(curve) > MAXNAN:
                print(f"{nNaN(curve)} > {MAXNAN} NaNs. SKIP {fname}")
                continue
            else:
                print(f"{nNaN(curve)} <= {MAXNAN} NaNs. Allowing {fname}")
        
            try: curve = interpolate_nan(curve)
            except NotImplementedError:
                print("Not implemented first or last NaN point.")
                continue
        
        curve_interp, interp_idx = interpolate(curve, DIST)
        assert nNaN(curve_interp) == 0, nNaN(curve_interp)
        
        df = pd.DataFrame(curve_interp, columns=list("xyz"))
        df["interp"] = interp_idx.astype(int)
        df.to_csv(f"{OUTDIR}/{name}.tsv", sep='\t', index=False)
        print(f"{name}:\t{len(curve)}\t->\t{len(curve_interp)}")

if __name__ == "__main__":
    main()

