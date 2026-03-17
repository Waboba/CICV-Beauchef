import networkx as nx
import numpy as np
from random import random, choice, choices
from copy import deepcopy
from math import floor, ceil
import matplotlib.pyplot as plt

from generator.graphlet_counts import update_counts_edge, update_counts_node


##### ---------- GRAIP ---------- #####


### First some utility functions.

def bin_degrees(deg, std, samples):
    """
    Bin the degree histogram such that no bin contains less than one node.

    Parameters
    ----------
    deg : array
        Mean degree histogram.
    std : array
        Standard deviation on degree histogram.
    samples : int
        Number of samples used to obtain the degree statistics.

    Returns
    -------
    bins : array
        List of bin edges. The right edge is not included in the bin, except for the last one.
    weights : array
        List of weights of the degree histgoram.
    bounds : array
        Bounds on the degree histogram of the generated graph.
    """
    
    bins = [len(deg)]
    hist = []
    
    one = 1-0.1/samples             # Avoid issues with rounding
    while bins[-1] > 1:
        nodes = 0
        bin_edge = bins[-1]
        zeros = False
        while deg[bin_edge-1] == 0:
            zeros = True
            bin_edge -= 1
        
        if not zeros:
            while nodes < one and bin_edge > 1:
                bin_edge -= 1
                nodes += deg[bin_edge]
        
        bins.append(bin_edge)
        if abs(nodes-1) < 0.1/samples:
            hist.append(1)
        else:
            hist.append(nodes)
    
    bins = np.array(bins[::-1])
    hist = np.array(hist[::-1])
    
    # Bin widths currently depend on number of samples, because peaks will be broadened at larger sample size.
    # Therefore, exclude parts that together add up to less than 1% of sample size.
    
    limit = 0.01
    for i in range(len(bins)-1):
        if bins[i+1]-bins[i] > 1 and hist[i] != 0:
            # Fix lower bound
            lower_part = deg[bins[i]] + deg[bins[i]+1]
            while lower_part < limit:
                bins[i] += 1
                lower_part += deg[bins[i]+1]
            # Fix upper bound
            upper_part = deg[bins[i+1]-1] + deg[bins[i+1]-2]
            while upper_part < limit:
                bins[i+1] -= 1
                upper_part += deg[bins[i+1]-2]
    
    # Standard deviation is not used to define bounds for bins with width > 1.
    
    bounds = np.zeros_like(hist)
    for i in range(len(bins)-1):
        width = bins[i+1] - bins[i]
        if width == 1:
            bounds[i] = std[i+1]*2
        else:
            bounds[i] = max(hist[i]-floor(hist[i]), ceil(hist[i])-hist[i])*1.1      # Slightly larger because number of nodes might be different
    
    return bins, hist, bounds


def custom_degree_histogram(deg_hist, bins):
    custom_deg_hist, _ = np.histogram(np.arange(len(deg_hist)), bins=bins, weights=deg_hist)
    return custom_deg_hist


### Get the number of edges between neighbors of a given node.

def neighbor_edges(G, n):
    neighbors = list(G.neighbors(n))
    if len(neighbors) < 2:
        return 0
    edges = 0
    for i, n1 in enumerate(neighbors):
        for n2 in neighbors[i+1:]:
            if n1 in G._adj[n2]:
                edges += 1
    return edges



def GRAIP(properties, samples, graphlet_config, max_steps=None, node_step=5, w=2/3, max_rej=None, savepath=None):
    """
    GRAphlet-based Incremental generator for Probabilistic networks.

    Parameters
    ----------
    properties : list
        List of graph properties as returned by the sample function.
    samples : int
        Number of samples used.
    graphlet_config : list
        Information about the graphlet types to be taken into account. This contains:
            - List of graphlet names
            - Dictionary with valid bit-strings per graphlet, with the above names as keys
            - The corresponding count function
    max_steps : int, optional
        Maximum number of iterations. If None, the maximum is set to 100 times the mean number of edges.
    node_step : int, optional
        A node is added or removed every node_step steps. The default is 5.
    w : float, optional
        Weight factor for cost function. The default is 2/3.
    max_rej : int, optional
        Number of rejected steps in a row before a step is automatically accepted. If None, this is set to 2% of the mean number of edges.
    savepath : str, optional
        If given, the generated graph is saved to the path in edgelist format. The default is None. 

    Returns
    -------
    H : NetworkX graph
        The newly generated graph.
    """
    
    # The cost function
    
    def Cost(P, gl, P_target, P_bounds, E_gl, bounds_gl, w):
        
        P_sum = np.cumsum(P[::-1])
        Pt_sum = np.cumsum(P_target[::-1])
        cost_deg = np.sum(np.abs(P_sum-Pt_sum)/Pt_sum)/len(Pt_sum)
        
        cost_gl = 0
        for i in range(len(gl)):
            c = gl[i]
            E_c = E_gl[i]
            b_c = bounds_gl[i]
            
            if E_c == 0:
                continue
            p = b_c/E_c
            if c == 0 and E_c > b_c:                      # Count is 0, take 0.1 to avoid issues with log
                cost_gl += np.log(0.1/E_c)/np.log(1-p)
            elif c < E_c - b_c:
                cost_gl += np.log(c/E_c)/np.log(1-p)
            elif c > E_c + b_c:
                cost_gl += np.log(c/E_c)/np.log(1+p)
        cost_gl /= len(gl)
        
        if np.any(gl[E_gl==0]!=0):
            cost_gl *= 10                 # Penalty for having a graphlet that never appears in the target graph
        
        cost = w*cost_deg + (1-w)*cost_gl
        
        return cost
    
    
    # Preprocessing
    
    E_n, std_n, E_e, std_e, E_gl, std_gl, E_deg, std_deg = properties
    
    E_n = float(E_n)
    std_n = float(std_n)
    E_e = float(E_e)
    std_e = float(std_e)
    
    E_deg = np.trim_zeros(E_deg, 'b')
    std_deg = np.trim_zeros(std_deg, 'b')
    
    bins, hist, bounds = bin_degrees(E_deg, std_deg, samples)
    P_target = hist/E_n
    P_bounds = bounds/E_n
    
    bounds_gl = std_gl*2
    
    
    Cg = 3*E_gl[1]/(E_gl[0]+3*E_gl[1])      # Global clustering coefficient
    
    if max_steps is None:
        max_steps = round(E_e*100)
        print('Maximum number of iterations set to %s (100 x E_e) by default.' %max_steps)
    
    graphlets, valid, count_func = graphlet_config
    
    if max_rej == None:
        max_rej = round(0.02*E_e)
    
    
    # Construct seed graph.
    # A really bad seed can very rarely slow down graph generation a lot, so generate 10 and keep the best one to be safe.
    
    N = round(E_n*0.2)
    lowest_cost = 1000
    for _ in range(10):
        H_test = nx.barabasi_albert_graph(N, round(E_e/E_n))
        gl_test = count_func(H_test)
        P_test = custom_degree_histogram(nx.degree_histogram(H_test), bins)/N
        cost = Cost(P_test, gl_test, P_target, P_bounds, E_gl, bounds_gl, w)
        if cost < lowest_cost:
            lowest_cost = cost
            H = deepcopy(H_test)
    gl = count_func(H)
    deg_hist = np.array(nx.degree_histogram(H))
    H_cost = lowest_cost
    
    
    # Incremental graph generation
    
    steps = 0
    counter = 0
    node_step_counter = 0
    new_node = N        # For convenience, the index of newly added nodes is not decreased when a node is removed.
    
    while True:
        
        if node_step_counter == node_step:
            node_step_counter = 0
            
            node_deviation = (N - E_n)/std_n
            r = random()
            node_addition = r > 1/(1+np.exp(-node_deviation))
            
            if node_addition:                     # Add a node
                
                n = new_node
                weights = np.array([degree for node, degree in H.degree()])
                first_neighbor = choices(list(H.nodes()), weights=weights/np.sum(weights))[0]
                
                largest_clique = max(nx.find_cliques(H, nodes=[first_neighbor]), key=len) 
                
                if len(largest_clique) >= 4:    # First neighbor is part of a clique size >= 4
                    for nb in largest_clique:
                        H.add_edge(n, nb)
                
                else:                           # No clique, match global CC
                    ne = neighbor_edges(H, first_neighbor)
                    k0 = H.degree(first_neighbor)
                    prob = 0.5*(k0+1)*Cg - ne/k0
                    for nb in H.neighbors(first_neighbor):
                        if random() < prob:
                            H.add_edge(n, nb)
                    H.add_edge(n, first_neighbor)
                
                deg_hist_temp = np.array(nx.degree_histogram(H))
            
            
            else:                               # Remove a node
                n = choice(list(H.nodes()))
                deg_hist_temp = np.copy(deg_hist)
                deg_hist_temp[H.degree(n)] -= 1
                for nb in H.neighbors(n):
                    deg_hist_temp[H.degree(nb)] -= 1
                    deg_hist_temp[H.degree(nb)-1] += 1
                
            gl_temp = gl.copy()
            gl_delta = update_counts_node(H, n, graphlets, valid)
            if node_addition:
                gl_temp += gl_delta
            else:
                gl_temp -= gl_delta
        
        else:
            node_step_counter += 1
        
            edge_deviation = (H.number_of_edges()*E_n/N - E_e)/std_e
            r = random()
            edge_addition = r > 1/(1+np.exp(-edge_deviation))
            
            if edge_addition:               # Add an edge
                while True:
                    n1, n2 = np.random.choice(H.nodes(), 2, replace=False)
                    if n1 not in H._adj[n2]:
                        break
        
            else:                           # Remove an edge
                while True:
                    n1 = np.random.choice(H.nodes())
                    if H.degree(n1) >= 1:
                        break
                n2 = np.random.choice(list(H.neighbors(n1)))
        
            deg_hist_temp = np.array(nx.degree_histogram(H)+[0])        # Extra element in case max degree is increased
            
            gl_temp = gl.copy()
            gl_temp += update_counts_edge(H, (n1,n2), graphlets, valid)
            deg_hist_temp[H.degree(n1)] -= 1
            deg_hist_temp[H.degree(n2)] -= 1
            if edge_addition:
                deg_hist_temp[H.degree(n1)+1] += 1
                deg_hist_temp[H.degree(n2)+1] += 1
            else:
                deg_hist_temp[H.degree(n1)-1] += 1
                deg_hist_temp[H.degree(n2)-1] += 1
        
        P_temp = custom_degree_histogram(deg_hist_temp, bins)/H.number_of_nodes()
        temp_cost = Cost(P_temp, gl_temp, P_target, P_bounds, E_gl, bounds_gl, w)
        
        
        if H_cost > temp_cost or counter == max_rej:
            counter = 0
            H_cost = temp_cost
            if node_step_counter == 0:
                if node_addition:
                    N += 1
                    new_node += 1
                else:
                    H.remove_node(n)
                    N -= 1
            else:
                if edge_addition:
                    H.add_edge(n1,n2)
                else:
                    H.remove_edge(n1,n2)
            deg_hist = deg_hist_temp.copy()
            gl = gl_temp.copy()
        else:
            if node_step_counter == 0 and node_addition:
                H.remove_node(n)
            counter += 1
        
        # Every once in a while, remove parts disconnected from the main component (these will be small).
        
        if steps%round(E_e) == 0:
            if nx.number_connected_components(H) > 1:
                H = H.subgraph(max(nx.connected_components(H), key=len)).copy()
                N = H.number_of_nodes()
                gl = count_func(H)
                P = custom_degree_histogram(nx.degree_histogram(H), bins)/N
                cost = Cost(P, gl, P_target, P_bounds, E_gl, bounds_gl, w)
        
        steps += 1
        
        
        P = custom_degree_histogram(deg_hist, bins)/N
        
        if steps >= max_steps or (np.all(P>=P_target-P_bounds) and np.all(P<=P_target+P_bounds) and np.all(gl>=E_gl-bounds_gl) and np.all(gl<=E_gl+bounds_gl)):
            
            H = H.subgraph(max(nx.connected_components(H), key=len)).copy()
            
            if steps >= max_steps:
                print("Maximum number of iterations reached.")
            else:
                print("Graph generation finished after %s iterations." %steps)
            
            if savepath is not None:
                nx.write_edgelist(H, savepath, data=False)
                
            return H



##### ---------- SwapCon ---------- #####


def SwapCon(properties, count_func, temperature=0.01, cooling=0.99, threshold=0.05, max_reject=None, savepath=None):
    """
    The swapping model extended to uncertain networks.

    Parameters
    ----------
    properties : list
        List of graph properties as returned by the sample function.
    count_func : funcion
        The graphlet counting function
    temperature : float, optional
        Initial temperature for the simulated annealing algorithm. The default is 0.01.
    cooling : float, optional
        Cooling factor for the simulated annealing algorithm.. The default is 0.99.
    threshold : float, optional
        Energy threshold for convergence of the simulated annealing algorithm. The default is 0.05.
    max_reject : int, optional
        Alternative convergence criterion: the algorithm is stopped after max_reject swaps in a row have been rejected. If None, this is set to the mean number of edges.
    savepath : str, optional
        If given, the generated graph is saved to the path in edgelist format. The default is None.

    Returns
    -------
    H : NetworkX graph
        The newly generated graph.
    """
    
    # The energy function.
    
    def Energy(gl, E_gl):
        N = len(gl)
        energy = 0
        for i in range(N):
            if gl[i] != 0 or E_gl[i] != 0:
                energy += abs(gl[i]-E_gl[i])/(gl[i]+E_gl[i])
        return energy/N
    
    # Function to generate a degree sequence.
    
    def generate_stubs(deg, E_n, P):
        while True:
            stubs = np.random.choice(deg, size=E_n, p=P)
            if nx.is_graphical(stubs):
                return stubs
    
    
    # Preprocessing
    
    E_n, std_n, E_e, std_e, E_gl, std_gl, E_deg, std_deg = properties
    
    E_n = round(float(E_n))
    
    E_deg = np.trim_zeros(E_deg, 'b')
    P = E_deg/np.sum(E_deg)
    deg = np.arange(len(P))
    
    if max_reject is None:
        max_reject = round(float(E_e))
    
    # Generate initial graph with configuration model.
    
    # First extract a degree sequance.
    
    stubs = generate_stubs(deg, E_n, P)
    
    # Now randomly connect stubs.
    
    H = nx.Graph()
    H.add_nodes_from(range(E_n))
    stubs_full = stubs.copy()
    new_seq = 0
    while np.any(stubs != 0):
        t = 0
        while True:
            if np.count_nonzero(stubs) < 2 or t == 100:     # Algorithm got stuck, restart from empty graph.
                if new_seq == 100:      # Restarted 100 times in a row, extract a new degree sequence.
                    new_seq = 0
                    stubs_full = generate_stubs(deg, E_n, P)
                H = nx.Graph()
                H.add_nodes_from(range(E_n))
                stubs = stubs_full.copy()
                new_seq += 1
            n1, n2 = np.random.choice(np.arange(E_n), size=2, replace=False, p=stubs/sum(stubs))
            if n1 not in H._adj[n2]:
                break
            else:
                t += 1
        
        stubs[n1] -= 1
        stubs[n2] -= 1
        H.add_edge(n1,n2)
    
    gl = count_func(H)
    energy_H = Energy(gl, E_gl)
    
    
    # Swap edges randomly.
    
    steps = 0
    reject = 0
    while True:
        
        while True:
            n1, n2 = choice(list(H.edges()))
            n3, n4 = choice(list(H.edges()))
            if (n1,n2) != (n3,n4) and n1 != n4 and n2 != n3 and n1 not in H._adj[n4] and n2 not in H._adj[n3]:
                break
        
        T = deepcopy(H)
        
        T.remove_edge(n1,n2)
        T.remove_edge(n3,n4)
        T.add_edge(n1,n4)
        T.add_edge(n2,n3)
        
        gl_T = count_func(T)
        energy_T = Energy(gl_T, E_gl)
        
        r = random()
        if energy_T < energy_H or (energy_T != energy_H and r < np.exp((energy_H-energy_T)/temperature)):
            reject = 0
            H = deepcopy(T)
            gl = gl_T.copy()
            energy_H = energy_T
        else:
            reject += 1
        
        steps += 1
        temperature *= cooling
        
        if reject >= max_reject or energy_H <= threshold:
            
            H = H.subgraph(max(nx.connected_components(H), key=len)).copy()
            
            if reject >= max_reject:
                print("Maximum number of rejected graphs in a row reached after %s swaps." %steps)
            else:
                print("Energy threshold reached after %s swaps." %steps)
            
            if savepath is not None:
                nx.write_edgelist(H, savepath, data=False)
                
            return H



##### ---------- (dual) BA model ---------- #####


def BA_graph(E_n, std_n, E_e, savepath=None):
    """
    The BA model used for comparison in the paper.
    """
    
    n = np.random.normal(E_n, std_n)
    e = n*E_e/E_n
    
    n_ER = round(0.1*n)
    e_ER = round(0.1*e)
    seed = nx.gnm_random_graph(n_ER, e_ER)      # Start from a small ER graph (NetworkX default is star graph)
    n -= n_ER
    e -= e_ER
    
    m = int(e/n)
    n1 = (m+1)*n - e
    
    # Add n1 nodes with m edges and (n-n1) nodes with m+1 edges.
    
    H = nx.dual_barabasi_albert_graph(n+n_ER, m, m+1, n1/n, initial_graph=seed)
    H = H.subgraph(max(nx.connected_components(H), key=len)).copy()
    
    if savepath is not None:
        nx.write_edgelist(H, savepath, data=False)
    
    return H

