import pandas as pd
import numpy as np
import yaml
import copy
import matplotlib
import matplotlib.pyplot as plt
import warnings
import random
import collections
import networkx as nx
import time
import joblib
from joblib import Parallel, delayed
import itertools
import pickle, gzip, pickletools
import os
# Input File (with disruption)
nodes=pd.read_csv("nodes_new_s.csv")
edges=pd.read_csv("edges_new_s.csv") 

# Original Network
G = nx.Graph()
G.add_nodes_from(list(nodes.Node))
G.add_weighted_edges_from([(edges.iloc[x,2], edges.iloc[x,3], edges.iloc[x,-1]) for x in range(edges.shape[0])])

# Giant connected component
def calculate_gc_size(G):
    """
    gc = giant component
    """
    gc_size = [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)][0]
    return gc_size

# to get graph of specific data
def get_Graph_t (G, t, nodes, edges):
    if t == 0:
        return G
    Gt = G.copy()
    Gt.remove_edges_from([(np.array(edges.loc[ edges.iloc[:, 15+t] == 1][['Source', 'Target']])[i,0], np.array(edges.loc[ edges.iloc[:, 15+t] == 1][['Source', 'Target']])[i,1]) for i in range(len(np.array(edges.loc[ edges.iloc[:, 15+t] == 1][['Source', 'Target']])))])
    Gt.remove_nodes_from(list(nodes.loc[ nodes.iloc[:, 3+t] == 1 ]['Node']))
    return Gt

# to generate timeseries result in seperate folders
def PickleObj(obj, filepath):
    with gzip.open(filepath, "wb") as f:
        pickled = pickle.dumps(obj)
        optimized_pickle = pickletools.optimize(pickled)
        f.write(optimized_pickle)

def UnpickleObj(filepath):
    with gzip.open(filepath, 'rb') as f:
        p = pickle.Unpickler(f)
        obj = p.load()
    return obj

def tic():
    return time.time()

def tac(t_start):
    val_secs = round(time.time() - t_start)
    (t_min, t_sec) = divmod(val_secs, 60)
    (t_hour, t_min) = divmod(t_min, 60)
    return '{} hour: {} mins: {} secs'.format(t_hour, t_min, t_sec), val_secs

def my_func(G, t, nodes, edges):
    t_start = tic()
    G_t = get_Graph_t(G, t, nodes, edges)
    GCC = (calculate_gc_size(G_t))/57032 #(size of Giant connected component of original network)
    G_LCC_t = [G_t.subgraph(x).copy() for x in nx.connected_components(G_t)][0]
    
    G_LCC_t_dia = nx.diameter(G_LCC_t)
    # ‘dijkstra’, ‘bellman-ford’, ‘floyd-warshall’ and ‘floyd-warshall-numpy’
    shortest_path=nx.average_shortest_path_length(G_LCC_t, weight='weight', method='dijkstra')
    time_run, val_secs = tac(t_start)
    filepath = "SL_Data_Sim/{x}.pkl".format(x = t)
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    print([len(list(G.nodes)),t, GCC, len(list(G_LCC_t.nodes)), G_LCC_t_dia, shortest_path, time_run])
    PickleObj([t, GCC, G_LCC_t_dia, shortest_path, val_secs], filepath)
    print(filepath)
#'''
with Parallel(n_jobs = 32, verbose = 10) as parallel:
    parallel(delayed(my_func)(G, t, nodes, edges) for t in range(1,62+1))
#'''


