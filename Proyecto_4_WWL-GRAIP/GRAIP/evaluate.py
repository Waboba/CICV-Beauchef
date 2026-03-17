import networkx as nx
import numpy as np
from random import random
from sklearn.metrics.pairwise import pairwise_kernels
from functools import partial
from generator.twwl import Twwl

from generator.graphlet_counts import three_counts, four_counts, five_counts



##### ---------- Quality metrics ---------- #####


def gaussian(x, y, sigma=1.0):
    dist = np.linalg.norm(x - y, 2)
    return np.exp(-1*dist**2/(2*sigma**2))

def kernel_compute(metric, X, Y=None, is_hist=True):

    def preprocess(X, max_len, is_hist):
        X_p = np.zeros((len(X), max_len))
        for i in range(len(X)):
            X_p[i, :len(X[i])] = X[i]

        if is_hist:
            row_sum = np.sum(X_p, axis=1)
            X_p = X_p/row_sum[:, None]

        return X_p
    
    
    max_len = max([len(x) for x in X])
    if Y is not None:
        max_len = max(max_len, max([len(y) for y in Y]))
    X = preprocess(X, max_len, is_hist)

    if Y is not None:
        Y = preprocess(Y, max_len, is_hist)

    return pairwise_kernels(X, Y, metric=metric)

def compute_mmd(sampled, generated, metric, is_hist=True):

    X = kernel_compute(metric, sampled, is_hist=is_hist)
    Y = kernel_compute(metric, generated, is_hist=is_hist)
    Z = kernel_compute(metric, sampled, Y=generated, is_hist=is_hist)

    return np.mean(X) + np.mean(Y) - 2*np.mean(Z)


def mmd_degree(generated_graphs, sample_graphs=None, target_network=None, batches=1):
    """
    Compute the squared MMD score of the degree distribution.

    Parameters
    ----------
    generated_graphs : list
        A list of generated NetworkX graphs.
    sampled_graphs : list, optional
        A list of sample graphs. Should contain at least as many graphs as generated_graphs. If not provided, len(generated_graphs) samples are first obtained from the target network. The default is None.
    target_network : NetworkX graph, optional
        The probabilistc target network. Is only required if no list of sample graphs is provided. The default is None.
    batches : int, optional
        Graphs are divided in batches and the final score is the mean over the batches. Should be a factor of len(generated_graphs). The default is 1.

    Returns
    -------
    MMD_deg : float
        The MMD squared score, averaged over the batches.
    """
    
    if sample_graphs is None:
        if target_network is None:
            raise ValueError("Provide either a list of sample graphs or the probabilistic target network.")
        else:
            
            # Construct sample graphs
            
            sample_graphs = []
            for _ in range(len(generated_graphs)):
                Gi = nx.create_empty_copy(target_network)
                for edge in target_network.edges():
                    r = random()
                    if r <= target_network.edges[edge]["probability"]:
                        Gi.add_edge(edge[0], edge[1])
                Gi = Gi.subgraph(max(nx.connected_components(Gi), key=len))
                sample_graphs.append(Gi)
    
    # Compute mmd
    
    mmd_deg = []
    for b in range(batches):
        hist_generated = []
        hist_sampled = []
        for n in range(len(generated_graphs)//batches):
            hist_generated.append(nx.degree_histogram(generated_graphs[b*batches+n]))
            hist_sampled.append(nx.degree_histogram(sample_graphs[b*batches+n]))
        
        if b == 0:
            sigma_ = []
            for i in range(len(hist_sampled)-1):
                for j in range(i+1, len(hist_sampled)):
                    hi = np.array(hist_sampled[i])
                    hj = np.array(hist_sampled[j])
                    max_len = max(len(hi), len(hj))

                    # Pad both arrays to have the same length
                    hi = np.pad(hi, (0, max_len-len(hi)), 'constant', constant_values=0)
                    hj = np.pad(hj, (0, max_len-len(hj)), 'constant', constant_values=0)
                    sigma_.append(np.linalg.norm(hi-hj, 2))
            sigma = np.median(sigma_)
        
        mmd_deg.append(compute_mmd(hist_sampled, hist_generated, metric=partial(gaussian, sigma=sigma), is_hist=True))
        
    MMD_deg = np.mean(mmd_deg)
    
    print("MMD score on degree distribution: %s" %MMD_deg)
    
    return MMD_deg


def mmd_graphlets(generated_graphs, max_gl_size, sample_graphs=None, target_network=None, batches=1):
    """
    Compute the squared MMD score of the graphlet frequencies.

    Parameters
    ----------
    generated_graphs : list
        A list of generated NetworkX graphs.
    max_gl_size : int
        The maximum order of graphlets taken into account. This can be 3, 4 or 5.
    sampled_graphs : list, optional
        A list of sample graphs. Should contain at least as many graphs as generated_graphs. If not provided, len(generated_graphs) samples are first obtained from the target network. The default is None.
    target_network : NetworkX graph, optional
        The probabilistc target network. Is only required if no list of sample graphs is provided. The default is None.
    batches : int, optional
        Graphs are divided in batches and the final score is the mean over the batches. Should be a factor of len(generated_graphs). The default is 1.

    Returns
    -------
    MMD_deg : float
        The MMD squared score, averaged over the batches.
    """
    
    if max_gl_size == 3:
        count_func = three_counts
    elif max_gl_size == 4:
        count_func = four_counts
    elif max_gl_size == 5:
        count_func = five_counts
    else:
        raise ValueError("max_gl_size can only be 3, 4 or 5. Larger graphlet sizes not implemented.")
    
    if sample_graphs is None:
        if target_network is None:
            raise ValueError("Provide either a list of sample graphs or the probabilistic target network.")
        else:
            
            # Construct sample graphs
            
            sample_graphs = []
            for _ in range(len(generated_graphs)):
                Gi = nx.create_empty_copy(target_network)
                for edge in target_network.edges():
                    r = random()
                    if r <= target_network.edges[edge]["probability"]:
                        Gi.add_edge(edge[0], edge[1])
                Gi = Gi.subgraph(max(nx.connected_components(Gi), key=len))
                sample_graphs.append(Gi)
    
    # Compute mmd
    
    mmd_graphlets = []
    
    for b in range(batches):
        counts_generated = []
        counts_sampled = []
        for n in range(len(generated_graphs)//batches):
            counts_generated.append(count_func(generated_graphs[b*batches+n]))
            counts_sampled.append(count_func(sample_graphs[b*batches+n]))
        
        if b == 0:
            sigma_ = []
            for i in range(len(counts_sampled)-1):
                for j in range(i+1, len(counts_sampled)):
                    sigma_.append(np.linalg.norm(counts_sampled[i]-counts_sampled[j], 2))
            sigma = np.median(sigma_)
        
        mmd_graphlets.append(compute_mmd(counts_sampled, counts_generated, metric=partial(gaussian, sigma=sigma), is_hist=False))
        
    MMD_graphlets = np.mean(mmd_graphlets)
    
    print("MMD score on graphlet frequencies: %s" %MMD_graphlets)
    
    return MMD_graphlets




##### ---------- Randomness metrics ---------- #####


def spread_diameter(generated_graphs, sample_graphs=None, target_network=None, samples=10000):
    """
    Compute the relative spread on the graph diameter.

    Parameters
    ----------
    generated_graphs : list
        A list of generated NetworkX graphs.
    sampled_graphs : list, optional
        A list of sample graphs. Should contain at least as many graphs as generated_graphs. If not provided, len(generated_graphs) samples are first obtained from the target network. The default is None.
    target_network : NetworkX graph, optional
        The probabilistc target network. Is only required if no list of sample graphs is provided. The default is None.
    samples : int, optional
        Defines how many sample graphs are constructed if no list of sample graphs is provided. Otherwise, this argument is ignored. The default is 10000.
    
    Returns
    -------
    score : float
        The relative spread score.
    """
    
    if sample_graphs is None:
        if target_network is None:
            raise ValueError("Provide either a list of sample graphs or the probabilistic target network.")
        else:
            
            # Construct sample graphs
            
            sample_graphs = []
            for _ in range(samples):
                Gi = nx.create_empty_copy(target_network)
                for edge in target_network.edges():
                    r = random()
                    if r <= target_network.edges[edge]["probability"]:
                        Gi.add_edge(edge[0], edge[1])
                Gi = Gi.subgraph(max(nx.connected_components(Gi), key=len))
                sample_graphs.append(Gi)
    
    diam_generated = []
    diam_sampled = []
    
    for i in range(len(generated_graphs)):
        diam_generated.append(nx.diameter(generated_graphs[i]))
        diam_sampled.append(nx.diameter(sample_graphs[i]))
    
    score = (np.percentile(diam_generated, 95) - np.percentile(diam_generated, 5)) / (np.percentile(diam_sampled, 95) - np.percentile(diam_sampled, 5))
    
    print("Relative spread on graph diameters: %s" %score)
    
    return score


def spread_cc(generated_graphs, sample_graphs=None, target_network=None, samples=10000):
    """
    Compute the relative spread on the average local clusterin coefficient.

    Parameters
    ----------
    generated_graphs : list
        A list of generated NetworkX graphs.
    sampled_graphs : list, optional
        A list of sample graphs. Should contain at least as many graphs as generated_graphs. If not provided, len(generated_graphs) samples are first obtained from the target network. The default is None.
    target_network : NetworkX graph, optional
        The probabilistc target network. Is only required if no list of sample graphs is provided. The default is None.
    samples : int, optional
        Defines how many sample graphs are constructed if no list of sample graphs is provided. Otherwise, this argument is ignored. The default is 10000.
    
    Returns
    -------
    score : float
        The relative spread score.
    """
    
    if sample_graphs is None:
        if target_network is None:
            raise ValueError("Provide either a list of sample graphs or the probabilistic target network.")
        else:
            
            # Construct sample graphs
            
            sample_graphs = []
            for _ in range(samples):
                Gi = nx.create_empty_copy(target_network)
                for edge in target_network.edges():
                    r = random()
                    if r <= target_network.edges[edge]["probability"]:
                        Gi.add_edge(edge[0], edge[1])
                Gi = Gi.subgraph(max(nx.connected_components(Gi), key=len))
                sample_graphs.append(Gi)
    
    cc_generated = []
    cc_sampled = []
    
    for i in range(len(generated_graphs)):
        cc_generated.append(nx.average_clustering(generated_graphs[i]))
        cc_sampled.append(nx.average_clustering(sample_graphs[i]))
    
    score = (np.percentile(cc_generated, 95) - np.percentile(cc_generated, 5)) / (np.percentile(cc_sampled, 95) - np.percentile(cc_sampled, 5))
    
    print("Relative spread on average local cluster coefficient: %s" %score)
    
    return score

# Modificación que considera la nueva métrica de twwl

def mmd_twwl(generated_graphs, sample_graphs=None, target_network=None, H=5, lambda_=1.0, batches=1):
    """
    Compute the squared MMD score using the TWWL kernel.

    Parameters
    ----------
    generated_graphs : list
        A list of generated NetworkX graphs.
    sample_graphs : list, optional
        A list of sample graphs. Should contain at least as many graphs as generated_graphs.
        If not provided, len(generated_graphs) samples are first obtained from the target network.
    target_network : NetworkX graph, optional
        The probabilistic target network. Is only required if no list of sample graphs is provided.
    H : int, optional
        Number of WL iterations.
    lambda_ : float, optional
        Parameter for the exponential kernel exp(-lambda * distance).
    batches : int, optional
        Graphs are divided in batches and the final score is the mean over the batches.
        Should be a factor of len(generated_graphs).

    Returns
    -------
    mmd : float
        The MMD squared score, averaged over the batches.
    """
    if sample_graphs is None:
        if target_network is None:
            raise ValueError("Provide either a list of sample graphs or the probabilistic target network.")
        # Construir muestras (igual que en mmd_degree)
        sample_graphs = []
        for _ in range(len(generated_graphs)):
            Gi = nx.create_empty_copy(target_network)
            for edge in target_network.edges():
                r = random()
                if r <= target_network.edges[edge]["probability"]:
                    Gi.add_edge(edge[0], edge[1])
            Gi = Gi.subgraph(max(nx.connected_components(Gi), key=len))
            sample_graphs.append(Gi)

    # Asegurar que todos los nodos tengan atributo 'label' (si no, asignar grado)
    for G in generated_graphs + sample_graphs:
        for node in G.nodes:
            if 'label' not in G.nodes[node]:
                G.nodes[node]['label'] = str(G.degree(node))  # asignar grado como etiqueta

    # Inicializar y ajustar TWWL con todos los grafos (generados + muestras)
    all_graphs = generated_graphs + sample_graphs
    n_gen = len(generated_graphs)
    n_sam = len(sample_graphs)

    twwl = Twwl(H=H, label_attr='label', measure_attr=None)  # measure se generará automáticamente
    twwl.fit(all_graphs, update_measure=True, verbose=False)

    # Calcular matriz de distancias
    dist = twwl.tree_metric(verbose=False)

    # Dividir en batches si es necesario
    mmd_scores = []
    batch_size = n_gen // batches
    for b in range(batches):
        start = b * batch_size
        end = (b+1) * batch_size
        # Índices de los generados y muestras en la matriz global
        gen_idx = list(range(start, end))
        sam_idx = list(range(n_gen + start, n_gen + end))

        # Extraer submatrices
        K_gg = np.exp(-lambda_ * dist[np.ix_(gen_idx, gen_idx)])
        K_ss = np.exp(-lambda_ * dist[np.ix_(sam_idx, sam_idx)])
        K_gs = np.exp(-lambda_ * dist[np.ix_(gen_idx, sam_idx)])

        mmd = np.mean(K_gg) + np.mean(K_ss) - 2 * np.mean(K_gs)
        mmd_scores.append(mmd)

    mmd_avg = np.mean(mmd_scores)
    print("MMD score on TWWL kernel: %s" % mmd_avg)
    return mmd_avg
