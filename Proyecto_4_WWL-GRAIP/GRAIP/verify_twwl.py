import os
import numpy as np
import networkx as nx
from generator.twwl import Twwl

# Parámetros
H = 5
input_dir = "/home/gabarca/Escritorio/Main/Práctica III/Código/twwl/mutag_graphs" # misma carpeta que exportó Julia
output_file = "dist_python.csv"

def load_graph(filename):
    """Carga un grafo desde el formato exportado por Julia."""
    with open(filename, 'r') as f:
        n = int(f.readline().strip())
        labels = f.readline().strip().split()
        edges = []
        for line in f:
            u, v = map(int, line.strip().split())
            edges.append((u, v))
    
    G = nx.Graph()
    # Añadir nodos con etiquetas (los identificadores son 1..n)
    for i in range(1, n+1):
        G.add_node(i, label=labels[i-1])
    G.add_edges_from(edges)
    return G

# Cargar todos los grafos
graph_files = sorted([os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.txt')])
graphs = [load_graph(f) for f in graph_files]

# Calcular matriz de distancias
twwl = Twwl(H=H, label_attr='label', measure_attr=None)
twwl.fit(graphs, update_measure=True, verbose=True)
np.savetxt("seq_python_g1.csv", twwl._seqs[0], delimiter=',', fmt='%d')
dist = twwl.tree_metric(verbose=True)

# Guardar como CSV
np.savetxt(output_file, dist, delimiter=',')
print(f"Matriz de distancias guardada en {output_file}")