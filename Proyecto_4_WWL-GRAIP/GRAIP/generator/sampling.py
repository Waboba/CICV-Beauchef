import networkx as nx
import numpy as np
from random import random


def sample(G, N, n_gl, count_func, savepath=None):
    """
    Compute the mean and standard deviation on the number of nodes and edges, the degree histogram and the graphlet frequencies by sampling a probabilistic graph G.
    
    Parameters
    ----------
    G : NetworkX graph, of which each edge has an attribute 'probability' (float) between 0 and 1.
    N : int
        Number of samples.
    n_gl : int, optional
        Number of graphlets to be counted. Should be compatible with the count function:
        3-node graphlets: 2
        4-node graphlets: 8
        5-node graphlets: 29
    count_func : function
        The function used for counting graphlets. It should take a graph as input and return an array of graphlet counts.
    savepath : str, optional
        The results are saved to the given path in .npz file format. The default is None.

    Returns
    -------
    E_nodes : float
        Mean number of nodes.
    std_nodes : float
        Standard deviation on number of nodes.
    E_edges : float
        Mean number of edges.
    std_edges : float
        Standard deviation on number of edges.
    E_graphlets : array
        Mean graphlet frequencies.
    std_graphlets : array
        Standard deviation on graphlet frequencies.
    E_degrees : array
        Mean degree histogram.
    std_degrees : array
        Standard deviation on the degree histogram.
    """
    
    print("Sampling started.")
    
    if savepath is None:
        print("Warning: No save path given, results will not be stored.")
    
    # Initialize all variables
    
    E_nodes = 0
    E_nodes_sq = 0      # squared, store for std
    E_edges = 0
    E_edges_sq = 0
    E_graphlets = np.zeros(n_gl)
    E_graphlets_sq = np.zeros(n_gl)
    E_degrees = np.zeros(len(nx.degree_histogram(G)))
    E_degrees_sq = np.zeros(len(E_degrees))
    
    p = 1
    for i in range(N):
        
        # Construct sample graph
        
        Gi = nx.create_empty_copy(G)
        for edge in G.edges():
            r = random()
            if r <= G.edges[edge]["probability"]:
                Gi.add_edge(edge[0], edge[1])
        
        Gi = Gi.subgraph(max(nx.connected_components(Gi), key=len))
        
        # Compute properties
        
        E_nodes += Gi.number_of_nodes()
        E_nodes_sq += Gi.number_of_nodes()**2
        
        E_edges += Gi.number_of_edges()
        E_edges_sq += Gi.number_of_edges()**2
        
        counts = count_func(Gi).astype('float64')
        E_graphlets += counts
        E_graphlets_sq += counts**2
        
        deg_hist = np.array(nx.degree_histogram(Gi))
        deg_hist.resize(len(E_degrees))
        E_degrees += deg_hist
        E_degrees_sq += deg_hist**2
        
        # Logging
        
        progress = i/N*100
        if progress >= p*10:
            p += 1
            print("%s/%s samples done" %(i,N))
    
    E_nodes /= N
    std_nodes = np.sqrt(E_nodes_sq/N - E_nodes**2)
    
    E_edges /= N
    std_edges = np.sqrt(E_edges_sq/N - E_edges**2)
    
    E_graphlets /= N
    std_graphlets = np.sqrt(E_graphlets_sq/N - E_graphlets**2)
    
    E_degrees /= N
    std_degrees = np.sqrt(E_degrees_sq/N - E_degrees**2)
    
    print("Sampling finished.")
    
    if savepath is not None:
        np.savez(savepath, samples=N, E_nodes=E_nodes, std_nodes=std_nodes, E_edges=E_edges, std_edges=std_edges, E_graphlets=E_graphlets, std_graphlets=std_graphlets, E_degrees=E_degrees, std_degrees=std_degrees)
            

    return E_nodes, std_nodes, E_edges, std_edges, E_graphlets, std_graphlets, E_degrees, std_degrees