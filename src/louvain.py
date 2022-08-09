#!/usr/bin/env python3
import os, sys
import numpy as np
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

def communities(barcodes, representatives, nPoints) -> dict:
    G = create_graph_from_PH(barcodes, representatives, nPoints)
    return louvain(G, resolution=1)

if __name__ == "__main__":
    
    PH_dir, PC_dir, outfile = sys.argv[1:]
    
    partitions = {}
    
    print("name\tnPoints\tnCommunities")
    
    for filename in sorted(os.listdir(PH_dir)):
        filename = os.path.join(PH_dir, filename)
        name, ext = os.path.splitext(filename)
        if ext != ".json": continue
        name = os.path.basename(name)
        filename_PC = os.path.join(PC_dir, name + ".npy")
        if not os.path.isfile(filename_PC): continue
        
        b, r = load_PH(filename)
        n = len(np.load(filename_PC))
        partitions[name] = communities(b, r, n)
        # also works to show which curves analysis succeeded for
        print(name, n, len(set(partitions[name].values())), sep='\t')
    
    with open(outfile, 'w') as fh:
        json.dump(partitions, fh, indent=True)

