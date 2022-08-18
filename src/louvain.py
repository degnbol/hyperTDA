#!/usr/bin/env python3
import os, sys
import numpy as np
import pandas as pd
import networkx as nx
from community import best_partition as louvain
import json

"""
Partition nodes into communities using persistent homology representives. 
A graph is constructed by densely connecting nodes in each representative, 
weighing them by persistence and then applying Louvain partitioning.
USE: louvain.py PH/ pointClouds/ communities.json
It is assumed that files in PH/ are named the same as in pointClouds for each 
curve, except the file extension.
Point clouds are only used to get the total number of points.
Prints curve length and number of communities for each curve.
"""

def load_PH(filename):
    with open(filename) as fh: d = json.load(fh)
    barcode = np.array(d['barcode']).T
    representatives = d['representatives']
    return barcode, representatives

def create_graph_from_PH(barcodes, representatives, nPoints):
    G = nx.Graph()
    G.add_nodes_from(range(nPoints))
    
    for b, r in zip(barcodes, representatives):
        persistence = abs(b[1] - b[0])
        
        vx = []
        for el in r:
            vx.append(el[0])
            vx.append(el[1])
        vx = list(set(vx))
        
        # add edge for all unique pairs of nodes in this representative
        for k in range(len(vx)):
            for j in range(k+1, len(vx)):
                # -1 due to zero indexing
                G.add_edge(vx[k]-1, vx[j]-1, weight=persistence)
    
    return G

def communities(barcodes, representatives, nPoints):
    G = create_graph_from_PH(barcodes, representatives, nPoints)
    comms = louvain(G, resolution=1)
    # louvain returns dict mapping from node id to partition id.
    # Node ids are simply the range from 0 to nPoints so we can throw keys away 
    # by calling .values().
    assertMsg = "If assert fails then node ids aren't simply a range that can be discarded"
    assert np.all(np.asarray(list(comms.keys())) == np.arange(len(comms))), assertMsg
    return list(comms.values())

if __name__ == "__main__":
    
    PH_dir, PC_dir, outfile = sys.argv[1:]
    
    partitions = {}
    
    print("name\tnPoints\tnCommunities")
    
    for filename in sorted(os.listdir(PH_dir)):
        filename = os.path.join(PH_dir, filename)
        name, ext = os.path.splitext(filename)
        if ext != ".json": continue
        if os.path.getsize(filename) == 0:
            sys.stderr.write(f"Empty file: {filename}\n")
            continue
        name = os.path.basename(name)
        filename_PC = os.path.join(PC_dir, name)
        if os.path.isfile(filename_PC + ".npy"):
            n = len(np.load(filename_PC))
        elif os.path.isfile(filename_PC + ".tsv"):
            n = len(pd.read_table(filename_PC + ".tsv"))
        else:
            sys.stderr.write("File not found: "+filename_PC+".{npy,tsv}. nPoints can't be determined.\n")
            continue
        
        b, r = load_PH(filename)
        partitions[name] = communities(b, r, n)
        # also works to show which curves analysis succeeded for
        print(name, n, len(set(partitions[name])), sep='\t')
    
    with open(outfile, 'w') as fh:
        json.dump(partitions, fh, indent=True)

