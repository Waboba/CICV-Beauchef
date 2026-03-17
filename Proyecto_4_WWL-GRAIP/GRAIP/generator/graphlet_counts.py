import numpy as np
import networkx as nx
from itertools import combinations, product


#### This file contains all functions related to counting graphlets.

#### The combinatorial counting functions are our implementation of the ESCAPE algorithm, see:
#### Pinar, Ali, Comandur Seshadhri, and Vaidyanathan Vishal. "Escape: Efficiently counting all 5-vertex subgraphs." Proceedings of the 26th international conference on world wide web. 2017.


### First some simple utility functions.

def ordered(n1,n2):
    if n1 > n2:
        return (n2, n1)
    return (n1, n2)


class ZeroDict(dict):
    def __missing__(self, key):
        return 0

def topological_ordering(G):
    """
    Compute the topological ordering of an undirected graph G using the minimum vertex removal algorithm.

    Parameters
    ----------
    G : NetworkX graph
        An undirected graph.

    Returns
    -------
    DAG : NetworkX digraph
        The DAG version of G, in which each edge of G is directed (ni -> nj for i < j)
    """

    DAG = nx.DiGraph()
    DAG.add_nodes_from(G.nodes())
    
    neighbors = {n:list(G.neighbors(n)) for n in G.nodes()}
    degrees = dict(G.degree())
    deg_list = [[] for _ in range(max(degrees.values())+1)]
    
    for node in G.nodes():
        deg_list[degrees[node]].append(node)        # deg_list[i] is list of vertices with degree i
    min_degree = min(degrees.values())
    
    for _ in range(G.number_of_nodes()):
        
        while len(deg_list[min_degree]) == 0:       # deg_list[min_degree] may be empty due to removal of vertices
            min_degree += 1
        
        source = deg_list[min_degree].pop()         # Get vertex of minimal degree still in G
        for node in neighbors[source]:        # Loop over neighbors of source node
        
            # Update some lists to 'remove' source node from G
            
            deg = degrees[node]
            deg_list[deg].remove(node)
            deg_list[deg-1].append(node)
            if deg-1 < min_degree:
                    min_degree -= 1
            degrees[node] -= 1
            neighbors[node].remove(source)
            
            DAG.add_edge(source, node)          # Add directed edge (source -> node) to DAG
        
        del neighbors[source]
    
    return DAG


### Now the real counting functions.


##### ---------- 3-node graphlets ---------- #####


def three_counts(G):
    """
    Count all three-node (induced) graphlets.

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    counts : array
        Counts of the 2 three-node graphlets.
    """
    
    wedge = 0
    triangle = 0
    
    for n1 in G.nodes():                                    # Loop over all nodes n1
        for n2, n3 in combinations(G.neighbors(n1), 2):     # Loop over neighbors of n1, in pairs of 2
            if n3 in G.neighbors(n2):
                triangle += 1             # If neighbors are connected, (n1,n2,n3) froms a triangle
            else:
                wedge += 1                # Otherwise, this is a wedge
    return np.array([wedge, triangle//3])


### The following functions are used when counting higher order graphlets.

### Count the number of (non-induced) directed wedges between n1 and n2.
### There are three types:  outout: n1 <- n3 -> n2
###                         inin:   n1 -> n3 <- n2
###                         inout:  n1 -> n3 -> n2

def directed_wedges(DG):
    outout = ZeroDict()
    inin = ZeroDict()
    inout = ZeroDict()
    
    for node in DG.nodes():
        for (n1, n2) in combinations(DG.successors(node), 2):                 # Look for outout wedges
            pair = ordered(n1,n2)
            if pair in outout:
                outout[pair] += 1
            else:
                outout[pair] = 1
        
        for (n1, n2) in combinations(DG.predecessors(node), 2):               # Look for inin wedges
            pair = ordered(n1,n2)
            if pair in inin:
                inin[pair] += 1
            else:
                inin[pair] = 1
        
        for (n1, n2) in product(DG.predecessors(node), DG.successors(node)):  # Look for inout wedges
            pair = ordered(n1,n2)
            if pair in inout:
                inout[pair] += 1
            else:
                inout[pair] = 1
    
    return outout, inin, inout


### Count the number of (non-induced) wedges between n1 and n2.

def all_wedges(G):
    wedges = ZeroDict()
    
    for node in G.nodes():
        for (n1, n2) in combinations(G.neighbors(node), 2):
            pair = ordered(n1,n2)
            if pair in wedges:
                wedges[pair] += 1
            else:
                wedges[pair] = 1
    return wedges


### Count the number of triangles adjacent to each vertex.
### This function takes the DAG version of G as input, to avoid double counting.

def triangle_info(DG):
    tri_vertex = {n:0 for n in DG.nodes()}

    for n1 in DG.nodes():
        for n2, n3 in combinations(DG.successors(n1), 2):             # All triangles contain an outout wedge
            if n2 in DG._adj[n3] or n3 in DG._adj[n2]:
                tri_vertex[n1] += 1
                tri_vertex[n2] += 1
                tri_vertex[n3] += 1
    
    return tri_vertex



##### ---------- 4-node graphlets ---------- #####


def four_counts(G):
    """
    Count all three- and four-node (induced) graphlets.

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    counts : array
        Counts of the 2 three-node graphlets and 6 four-node graphlets.
    """
    
    star = 0
    path = 0
    tailed_tri = 0
    cycle = 0
    diamond = 0
    clique = 0
    
    # Compute everything that is needed later
    
    DG = topological_ordering(G)
    wedges = all_wedges(G)
    
    # 3-node graphlets
    
    W = sum(wedges.values())
    T = 0
    
    # Cut is a vertex
    
    for n1 in G.nodes():
        deg1 = G.degree(n1)
        
        star += deg1*(deg1-1)*(deg1-2)//6           # Pick 3 neighbors of n1 to form a star
        
        # Cut is an edge
        
        for n2 in DG.successors(n1):                # Only loop over successors to avoid double counting
            deg2 = G.degree(n2)
            w12 = wedges[ordered(n1,n2)]            # Number of triangles adjacent to this edge
            
            T += w12
            path += (deg1-1)*(deg2-1)               # Pick a different neighbor of both n1 and n2
            tailed_tri += w12*(deg1 + deg2 - 4)     # Pick a triangle adjacent to the edge and a neighbor of n1 or n2
            diamond += w12*(w12-1)//2               # Pick 2 triangles adjacent to the edge
        
            # Cliques:
            # There is only 1 directed clique (up to isomorphisms), so use its specific structure to avoid double counting.
            
            for n3 in set(DG.successors(n1)) & set(DG.successors(n2)):
                clique += len(set(DG.predecessors(n1)) & set(DG.predecessors(n2)) & set(DG.predecessors(n3)))
    
    # Cycles
    
    for w in wedges.values():
        cycle += w*(w-1)//2                         # Pick 2 wedges adjacent to a pair of nodes
    
    # Corrections:
    # Some graphlets were counted more than once and sometimes a lower order graphlets was counted as well.
    # E.g., when counting 4-paths by counting the number of neighbors of each endpoint of an edge,
    # a triangle (where the neighbors of both endpoints refer to the same node) is counted 3 times as well.
    
    T //= 3
    path -= 3*T
    tailed_tri //= 2
    cycle //= 2
    
    # Collect counts and transform into induced.
    
    three_induced_counts = [W-3*T, T]
    four_induced_counts = [star-tailed_tri+2*diamond-4*clique, path-2*tailed_tri-4*cycle+6*diamond-12*clique, tailed_tri-4*diamond+12*clique, cycle-diamond+3*clique, diamond-6*clique, clique]
    
    return np.concatenate((three_induced_counts, four_induced_counts))




##### ---------- 5-node graphlets ---------- #####


def five_counts(G):
    """
    Count all three-, four- and five-node (induced) graphlets.

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    counts : array
        Counts of the 2 three-node graphlets, 6 four-node graphlets and 21 five-node graphlets.
    """
    
    star = 0
    prong = 0
    path = 0
    fork_tailed_tri = 0
    long_tailed_tri = 0
    double_tailed_tri = 0
    tailed_cycle = 0
    hourglass = 0
    cycle = 0
    cobra = 0
    stingray = 0
    hatted_cycle = 0
    three_wedge = 0
    three_tri = 0
    tailed_clique = 0
    triangle_strip = 0
    diamond_wedge = 0
    wheel = 0
    hatted_clique = 0
    bipyramid = 0
    five_clique = 0
    
    # Compute everything that is needed later
    
    DG = topological_ordering(G)
    outout, inin, inout = directed_wedges(DG)
    tri_vertex = triangle_info(DG)
    
    # 3- and 4-node graphlets
    
    W = 0
    T = sum(tri_vertex.values())//3
    S4, P4, TT, C4, D, K4 = 0, 0, 0, 0, 0, 0
    
    # Cut is a vertex
    
    for n1 in G.nodes():
        deg1 = G.degree(n1)
        tri_v = tri_vertex[n1]
        
        S4 += deg1*(deg1-1)*(deg1-2)//6     # Pick 3 neighbors of n1
        TT += tri_v*(deg1-2)                # Pick a triangle adjacent to n2 and a neighbor not part of that triangle
        
        star += deg1*(deg1-1)*(deg1-2)*(deg1-3)//24         # Pick 4 neighbors of n1
        fork_tailed_tri += tri_v*(deg1-2)*(deg1-3)//2       # Pick triangle adjacent to n1 and 2 neighbors not part of that triangle
        hourglass += tri_v*(tri_v-1)//2                     # Pick 2 triangles adjacent to n1

        # Cut is an edge
        
        for n2 in DG.predecessors(n1):           # Only loop over predecessors to avoid double counting
            deg2 = G.degree(n2)
            pair12 = ordered(n1,n2)
            w12 = outout[pair12] + inin[pair12] + inout[pair12]
            
            P4 += (deg1-1)*(deg2-1)         # Pick a neighbor of n1 (except n2) and a neighbor of n2 (except n1)
            
            prong += (deg2-1)*(deg1-1)*(deg1-2)//2 + (deg1-1)*(deg2-1)*(deg2-2)//2  # Pick 1 neighbor of n1 (n2) and 2 neighbors of n2 (n1)
            double_tailed_tri += w12*(deg1-2)*(deg2-2)          # Pick a triangle adjacent to e12 and a neighbor of n1/n2 not part of that triangle
            stingray += w12*(w12-1)//2*(deg1-3 + deg2-3)        # Pick 2 triangles adjacent to e12 and a neighbor of n1 or n2 not part of either triangle
            three_tri += w12*(w12-1)*(w12-2)//6                 # Pick 3 triangles adjacent to e12

            # Cycle-related counts
            
            four_cycles = 0          # First count number of 4cycles adjacent to this edge
            
            for n3 in G.neighbors(n2):
                if n1 == n3:
                    continue
                pair13 = ordered(n1,n3)
                four_cycles += outout[pair13] + inin[pair13] + inout[pair13] - 1
            
            C4 += four_cycles
            
            tailed_cycle += four_cycles*(deg1-2 + deg2-2)       # Pick a cycle adjacent to e12 and a neighbor of n1 or n2 not part of that cycle
            hatted_cycle += w12*four_cycles                     # Pick a triangle and a cycle adjacent to e12
    
    # Cut is a wedge
    
    for (n1,n2) in outout.keys() | inin.keys() | inout.keys():
        count = outout[(n1,n2)] + inin[(n1,n2)] + inout[(n1,n2)]
        deg1 = G.degree(n1)
        deg2 = G.degree(n2)
        tri1 = tri_vertex[n1]
        tri2 = tri_vertex[n2]
        
        W += count
        path += count*(deg1-1)*(deg2-1)                         # Pick a neighbor, not part of the wedge, of each endpoint
        long_tailed_tri += count*(tri1 + tri2)                  # Pick a triangle adjacent to an endpoint
        three_wedge += count*(count-1)*(count-2)//6             # Pick three wedges
    
    # Diamond-related counts
    
    for n1 in G.nodes():
        for (n2,n3) in combinations(G.neighbors(n1), 2):
            dias = len(set(G.neighbors(n1)) & set(G.neighbors(n2)) & set(G.neighbors(n3)))
            if dias == 0:
                continue
            deg2 = G.degree(n2)
            deg3 = G.degree(n3)
            pair12 = ordered(n1,n2)
            pair13 = ordered(n1,n3)
            pair23 = ordered(n2,n3)
            w12 = outout[pair12] + inin[pair12] + inout[pair12]
            w13 = outout[pair13] + inin[pair13] + inout[pair13]
            w23 = outout[pair23] + inin[pair23] + inout[pair23]
            
            D += dias
            
            cobra += dias*(deg2-2 + deg3-2)             # Pick a diamond and a node adjacent to one of the unconnected nodes
            triangle_strip += dias*(w12-1 + w13-1)      # Pick a diamond and a triangle adjacent to one of the edges
            diamond_wedge += dias*(w23-2)               # Pick a diamond and an edge not part of the diamond
            wheel += dias*(dias-1)//2                   # Pick 2 diamonds adjacent to the triplet

    # Clique-related counts
    
    # Start from triangles.
    # There is only one type of DAG triangle, so use its specific structure to avoid double-counting.
    
    for n1 in G.nodes():
        deg1 = G.degree(n1)
        
        for n2 in DG.successors(n1):
            deg2 = G.degree(n2)
            pair12 = ordered(n1,n2)
            w12 = outout[pair12] + inin[pair12] + inout[pair12]
            
            for n3 in (set(DG.successors(n1)) & set(DG.successors(n2))):
                
                # A biparymid consists of two 4 cliques that share a triangle.
                # But there is no guarantee the fourth node of both 4cliques will be in intersect (see below).
                # We need to look for all common neighbors of n1, n2 and n3
                
                cliques = len(set(G.neighbors(n1)) & set(G.neighbors(n2)) & set(G.neighbors(n3)))
                if cliques == 0:
                    continue
                
                bipyramid += cliques*(cliques-1)//2
                
                # Every DAG 4-clique contains exactly one node for which the other 3 nodes are successors.
                
                intersect = set(DG.predecessors(n1)) & set(DG.predecessors(n2)) & set(DG.predecessors(n3))
                
                deg3 = G.degree(n3)
                pair13 = ordered(n1,n3)
                pair23 = ordered(n2,n3)
                w13 = outout[pair13] + inin[pair13] + inout[pair13]
                w23 = outout[pair23] + inin[pair23] + inout[pair23]

                for n4 in intersect:
                    deg4 = G.degree(n4)
                    pair14 = ordered(n1,n4)
                    pair24 = ordered(n2,n4)
                    pair34 = ordered(n3,n4)
                    w14 = outout[pair14] + inin[pair14] + inout[pair14]
                    w24 = outout[pair24] + inin[pair24] + inout[pair24]
                    w34 = outout[pair34] + inin[pair34] + inout[pair34]
                    
                    K4 += 1
                    tailed_clique += deg1 + deg2 + deg3 + deg4 - 12             # Pick a node adjacent to one of the 4clique nodes
                    hatted_clique += w12 + w13 + w14 + w23 + w24 + w34 - 12     # Pick a triangle adjacent to one of the 4clique edges
                
                for (n4,n5) in combinations(intersect, 2):
                    if n4 in G._adj[n5]:
                        five_clique += 1            # Pick two different 4cliques and check if the nodes are connected
    
    # Count of five-cycles based on DAG. This is the most annoying part.
    # Every five-cycle contains a n1 <- n3 <- n4 -> n2 pattern, with an additional wedge between n1 and n2.
    # This wedge can be inout (both directions) or outout (inin is already counted in inout).
    # We start from this wedge.
    
    dir_TT = 0          # Certain directed tailed triangles will also be counted and have to be corrected for
    
    for (n1,n2), count in inout.items():
        
        if n2 in DG._adj[n1]:
            dir_TT += count*(DG.out_degree(n1)-2)
        elif n1 in DG._adj[n2]:
            dir_TT += count*(DG.out_degree(n2)-2)
        
        # Now count possible n1 <- n3 <- n4 -> n2 paths.
        
        for n3 in DG.predecessors(n2):
            if n3 != n1:
                cycle += count*outout[ordered(n1,n3)]
            
        for n3 in DG.predecessors(n1):
            if n3 != n2:
                cycle += count*outout[ordered(n2,n3)]     
    
    for (n1,n2), count in outout.items():
        
        if n1 in DG._adj[n2]:
            dir_TT += count*DG.out_degree(n1)
            dir_TT += count*(DG.out_degree(n2)-1)
        if n2 in DG._adj[n1]:
            dir_TT += count*DG.out_degree(n2)
            dir_TT += count*(DG.out_degree(n1)-1)
        
        # Now count possible n1 <- n3 <- n4 -> n2 paths.
        
        for n3 in DG.predecessors(n2):
            if n3 != n1:
                cycle += count*outout[ordered(n1,n3)]
        
        for n3 in DG.predecessors(n1):
            if n3 != n2:
                cycle += count*outout[ordered(n2,n3)]
    
    
    # Corrections:
    # Some graphlets were counted more than once and sometimes a lower order graphlets was counted as well.
    # E.g., when counting 4-paths by counting the number of neighbors of each endpoint of an edge,
    # a triangle (where the neighbors of both endpoints refer to the same node) is counted 3 times as well.
    
    P4 -= 3*T
    C4 //= 4
    D //= 2
    
    prong -= 2*TT
    path -= (4*C4 + 2*TT + 3*T)
    long_tailed_tri -= (2*TT + 4*D +6*T)
    double_tailed_tri -= 2*D
    tailed_cycle = tailed_cycle//2 - 2*D
    cycle -= dir_TT
    hourglass -= 2*D
    cobra = cobra//2 - 12*K4
    diamond_wedge //= 2
    hatted_cycle -= 4*D
    triangle_strip = triangle_strip//2 - 12*K4
    wheel //= 2
    
    non_induced_counts = [star, prong, path, fork_tailed_tri, long_tailed_tri, double_tailed_tri, tailed_cycle, cycle, 
                          hourglass, cobra, stingray, hatted_cycle, three_wedge, three_tri, tailed_clique, triangle_strip, 
                          diamond_wedge, wheel, hatted_clique, bipyramid, five_clique]
    
    # Matrix transforming non-induced into induced counts
    
    transform = np.array([[1,0,0,-1,0,0,0,0,1,0,1,0,0,-2,-1,-1,0,1,2,-3,5],
                          [0,1,0,-2,-1,-2,-2,0,4,4,5,4,6,-12,-9,-10,-10,20,20,-36,60],
                          [0,0,1,0,-2,-1,-2,-5,4,4,2,7,6,-6,-6,-10,-14,24,18,-36,60],
                          [0,0,0,1,0,0,0,0,-2,0,-2,0,0,6,3,3,0,-4,-8,15,-30],
                          [0,0,0,0,1,0,0,0,-4,-2,0,-2,0,0,3,6,6,-16,-12,30,-60],
                          [0,0,0,0,0,1,0,0,0,-2,-2,-1,0,6,6,5,4,-12,-14,30,-60],
                          [0,0,0,0,0,0,1,0,0,-1,-1,-2,-6,6,3,4,8,-16,-12,30,-60],
                          [0,0,0,0,0,0,0,1,0,0,0,-1,0,0,0,1,2,-4,-2,6,-12],
                          [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,-1,0,2,2,-6,15],
                          [0,0,0,0,0,0,0,0,0,1,0,0,0,0,-3,-2,-2,8,8,-24,60],
                          [0,0,0,0,0,0,0,0,0,0,1,0,0,-6,-3,-2,0,4,10,-24,60],
                          [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-2,-4,12,6,-24,60],
                          [0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,-1,2,1,-4,10],
                          [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1,3,-10],
                          [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-2,6,-20],
                          [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-4,-4,18,-60],
                          [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-4,-1,9,-30],
                          [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-3,15],
                          [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-6,30],
                          [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-10],
                          [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]])
    
    induced_counts = transform @ non_induced_counts
    induced_counts = [int(count) for count in induced_counts]
    
    # Collect three- and four-node counts and transform into induced
    
    three_induced_counts = [W-3*T, T]
    four_induced_counts = [S4-TT+2*D-4*K4, P4-2*TT-4*C4+6*D-12*K4, TT-4*D+12*K4, C4-D+3*K4, D-6*K4, K4]
    
    return np.concatenate((three_induced_counts, four_induced_counts, induced_counts))



##### ---------- Update counts ---------- #####


### The following functions can be used to update the graphlet counts when a single node or edge is added or removed from a graph.


### First some auxiliary functions.

### This function converts a graphlet into a bit string that encodes which edges are present.
### Nodes are ordered at random. The least significant bit is 1 if there is an edge between node 1 and 2, 
### and 0 otherwise. The next bits refer to edges between nodes 1 and 3, 2 and 3, 1 and 4, ...
### Example: 4path with edges (1,2), (1,3), (2,4) -> 010|01|1 = 19

def get_code(G, nodes):
    N = len(nodes)
    code = 0
    for i in range(N-1):
        for j in range(i+1, N):
            if nodes[i] in G._adj[nodes[j]]:
                code |= 1<<(j*(j-1)//2+i)
    return code

### Update the code of an order k graphlet when a new node is added to form an order k+1 graphlet.
### The node ordering should be the same as for the original code, and the new node should appear last in 'nodes'.

def update_code(code, G, nodes):
    N = len(nodes)-1
    new_node = nodes[-1]
    for i in range(N):
        if nodes[i] in G._adj[new_node]:
            code |= 1<<(N*(N-1)//2+i)
    return code

### Given a dictionary of possible codes for each graphlet and a code, this function returns the name of the corresponding graphlet.
### If the code refers to an disconnected graphlet, the function will return None.

def code_to_graphlet(valid_codes, code):
    for name, val in valid_codes.items():
        if code in val:
            return name


### Now the update functions.


def update_counts_node(G, n1, graphlets, valid):
    """
    Compute the change in graphlet counts when a single node is added to or removed from the graph (along with edges linked to it).

    Parameters
    ----------
    G : NetworkX graph
        The graph AFTER the node is added or BEFORE the node is removed (i.e., it should contain the node).
    n1 : node
        The modified node.
    graphlets : list
        List of names of graphlets to be counted.
    valid : dictionary
        Precomputed valid bit-strings for each graphlet in 'graphlets'.

    Returns
    -------
    delta : array
        Change in graphlet counts (same size as 'graphlets'), such that new_counts = old_counts + delta
    """
    
    N = len(graphlets)

    delta = [0 for _ in range(N)]      # Record change in counts: new_counts = old_counts + delta
    
    # Keep track of visited nodes to avoid double counting.
    
    blacklist3 = set()
    blacklist4 = set()
    blacklist5 = set()
    
    for n2 in G.neighbors(n1):
        for s1 in (n1,n2):                      # Pick a source
            for n3 in G.neighbors(s1):          # Iterate over neighbors of source
                if n3 == n1 or n3 == n2:
                    continue
                pair = tuple(sorted((n2,n3)))
                if pair in blacklist3:
                    continue
                
                new3 = get_code(G, [n1,n2,n3])
                delta[graphlets.index(code_to_graphlet(valid, new3))] += 1
                
                blacklist3.add(pair)
                
                if N == 2:
                    continue
                
                for s2 in (n1,n2,n3):
                    for n4 in G.neighbors(s2):
                        if n4 == n1 or n4 == n2 or n4 == n3:
                            continue
                        triplet = tuple(sorted((n2,n3,n4)))
                        if triplet in blacklist4:
                            continue
                        
                        new4 = update_code(new3, G, [n1,n2,n3,n4])
                        delta[graphlets.index(code_to_graphlet(valid, new4))] += 1
                        
                        blacklist4.add(triplet)
                        
                        if N == 8:
                            continue
                        
                        for s3 in (n1,n2,n3,n4):
                            for n5 in G.neighbors(s3):
                                if n5 == n1 or n5 == n2 or n5 == n3 or n5 == n4:
                                    continue
                                quad = tuple(sorted((n2,n3,n4,n5)))
                                if quad in blacklist5:
                                    continue
                                
                                new5 = update_code(new4, G, [n1,n2,n3,n4,n5])
                                delta[graphlets.index(code_to_graphlet(valid, new5))] += 1
                                
                                blacklist5.add(quad)
    
    return np.array(delta, dtype=int)


### Function for edge updates. G should be the graph BEFORE the update.

def update_counts_edge(G, e, graphlets, valid):
    """
    Compute the change in graphlet counts when a single edge is added to or removed from the graph.

    Parameters
    ----------
    G : NetworkX graph
        The graph BEFORE the edge is added or removed.
    e : edge
        The modified edge.
    graphlets : list
        List of names of graphlets to be counted.
    valid : dictionary
        Precomputed valid bit-strings for each graphlet in 'graphlets'.

    Returns
    -------
    delta : array
        Change in graphlet counts (same size as 'graphlets'), such that new_counts = old_counts + delta
    """
    
    N = len(graphlets)
    n1, n2 = e
    
    delta = [0 for _ in range(N)]      # Record change in counts: new_counts = old_counts + delta
    
    # Keep track of visited nodes to avoid double counting.
    
    blacklist3 = set()
    blacklist4 = set()
    blacklist5 = set()
    
    for s1 in (n1,n2):                      # Pick a source
        for n3 in G.neighbors(s1):          # Iterate over neighbors of source
            if n3 == n1 or n3 == n2:
                continue
            if n3 in blacklist3:
                continue
            
            old3 = get_code(G, [n1,n2,n3])
            old_name = code_to_graphlet(valid, old3)
            if old_name != None:
                delta[graphlets.index(old_name)] -= 1
            
            new = old3^1
            new_name = code_to_graphlet(valid, new)
            if new_name != None:
                delta[graphlets.index(new_name)] += 1
            
            blacklist3.add(n3)
            
            if N == 2:
                continue
            
            for s2 in (n1,n2,n3):
                for n4 in G.neighbors(s2):
                    if n4 == n1 or n4 == n2 or n4 == n3:
                        continue
                    pair = tuple(sorted((n3,n4)))
                    if pair in blacklist4:
                        continue
                    
                    old4 = update_code(old3, G, [n1,n2,n3,n4])
                    old_name = code_to_graphlet(valid, old4)
                    if old_name != None:
                        delta[graphlets.index(old_name)] -= 1
                    
                    new = old4^1
                    new_name = code_to_graphlet(valid, new)
                    if new_name != None:
                        delta[graphlets.index(new_name)] += 1
                    
                    blacklist4.add(pair)
                    
                    if N == 8:
                        continue
                    
                    for s3 in (n1,n2,n3,n4):
                        for n5 in G.neighbors(s3):
                            if n5 == n1 or n5 == n2 or n5 == n3 or n5 == n4:
                                continue
                            triplet = tuple(sorted((n3,n4,n5)))
                            if triplet in blacklist5:
                                continue
                            
                            old5 = update_code(old4, G, [n1,n2,n3,n4,n5])
                            old_name = code_to_graphlet(valid, old5)
                            if old_name != None:
                                delta[graphlets.index(old_name)] -= 1
                            
                            new = old5^1
                            new_name = code_to_graphlet(valid, new)
                            if new_name != None:
                                delta[graphlets.index(new_name)] += 1
                            
                            blacklist5.add(triplet)
    
    return np.array(delta, dtype=int)
