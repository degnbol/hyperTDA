#!/usr/bin/env python3
import numpy as np
import json

def read_hyperedges(file):
    """
    Read hyperedges from file.
    Either read from file with integers separated by comma for a hyperedge on each line,
    OR read json. If json, then also accept field "representatives" as we have been using from PH.
    """
    with open(file) as fh:
        if file.endswith(".json"):
            d = json.load(fh)
            if "representatives" in d:
                d = d["representatives"]
            
            # flatten each element to a hyperedge, i.e. set of vertices
            return [list(np.asarray(he).flatten()) for he in d]
        else:
            hes = []
            for line in fh:
                l = line.strip().split(',')
                hes.append([int(v) for v in l])
            return hes


                

