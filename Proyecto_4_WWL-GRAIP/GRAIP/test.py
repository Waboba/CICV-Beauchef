import os
import networkx as nx
from random import random

from generator.generator import GraphGenerator
from evaluate import mmd_degree, mmd_graphlets, spread_diameter, spread_cc, mmd_twwl



#### An example script to show how to use the code. In this script, we will:
####    1. Generate a small probabilistic BA graph as target network (you will skip this step in practice, as you will already have a target network).
####    2. Derive its properties through sampling.
####    3. Generate graphs with GRAIP and SwapCon.
####    4. Evaluate the generated graphs.


os.mkdir("test")    # A new foldGuarda dist_julia en un archivo (por ejemplo, CSV o .npy mediante NPZ.jl).er to store the output
os.chdir("test")


# Construct target network

n = 100
m = 5
G = nx.barabasi_albert_graph(n, m)
for edge in G.edges:
    G.edges[edge]["probability"] = random()
nx.write_gml(G, "probabilistic_graph.gml")

# Guarda dist_julia en un archivo (por ejemplo, CSV o .npy mediante NPZ.jl).

# Now we will sample this network and generate new graphs.
# This can be done through the GraphGenerator class.

G = nx.read_gml("probabilistic_graph.gml", label=None)      # Read target graph
max_gl_size = 4     # We will consider graphlets of order 3 and 4

gen = GraphGenerator(G=G, max_gl_size=max_gl_size)

# Now perform sampling. We will construct 1000 samples and save the data for later use.

sample_path = "graph_data"

gen.get_properties(N=1000, savepath=sample_path, plot=True)    # Plots can be found under test/degree.PNG and test/graphlets.PNG


# Generate graphs with GRAIP and SwapCon, 10 for each model

graphs = 10
savepath = "new_graphs"     # Save new graphs in folder new_graphs. Graph files (edge list) will be named [model]_graph[n].txt

os.mkdir(savepath)

# Generate graphs with GRAIP. We will use mostly the default settings, but still provide them explicitly for clarity.
    
print("Generating GRAIP graphs...")

GRAIP_graphs = gen.generate(model='GRAIP', graphs=graphs, savepath=savepath, max_steps=10000, node_step=5, w=2/3, max_rej=None)

# Generate graphs with SwapCon.

print("Generating SwapCon graphs...")

SwapCon_graphs = gen.generate(model='SwapCon', graphs=graphs, savepath=savepath, temperature=0.01, cooling=0.99, threshold=0.05, max_reject=None)


# Now we will evaluate the graphs using the metrics used in our paper.
# We did not generate enough graphs to draw meaningful conclusions, but this is just an example.

print("Evaluating generated graphs...")

print("GRAIP:")

# Quality metrics:
mmd_degree(GRAIP_graphs, target_network=G, batches=1)
mmd_graphlets(GRAIP_graphs, max_gl_size, target_network=G, batches=1)
mmd_twwl(GRAIP_graphs, target_network=G, H=5, lambda_=1.0, batches=1)   # <--- NUEVA LÍNEA

# Randomness metrics:
spread_diameter(GRAIP_graphs, target_network=G, samples=1000)
spread_cc(GRAIP_graphs, target_network=G, samples=1000)


print("SwapCon:")

# Quality metrics:
mmd_degree(SwapCon_graphs, target_network=G, batches=1)
mmd_graphlets(SwapCon_graphs, max_gl_size, target_network=G, batches=1)
mmd_twwl(SwapCon_graphs, target_network=G, H=5, lambda_=1.0, batches=1)   # <--- NUEVA LÍNEA

# Randomness metrics:
spread_diameter(SwapCon_graphs, target_network=G, samples=1000)
spread_cc(SwapCon_graphs, target_network=G, samples=1000)

