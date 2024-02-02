#!/usr/bin/env python3
import leidenalg
import igraph as ig

def leiden(mat):
    G = ig.Graph.Weighted_Adjacency(mat)
    parts = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition)
    comm = np.zeros(mat.shape[0])
    for i, part in enumerate(parts):
        comm[part] = i+1
    return comm

