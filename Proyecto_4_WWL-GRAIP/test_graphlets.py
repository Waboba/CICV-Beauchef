import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from graphlets_counts import get_graphlet_report 

# 1. Definición de los grafos
aris1 = [(1,2),(2,3),(3,1),(3,4),(3,5),(5,6),(7,8),(8,5),(7,5)]
G1 = nx.Graph()
G1.add_edges_from(aris1)

aris2 = [(6,1),(1,2),(2,6),(3,5),(3,4),(3,5),(5,6),(4,7),(4,8),(7,8)]
G2 = nx.Graph()
G2.add_edges_from(aris2)

# 2. Obtener reportes M_1 a M_29
report_G1 = get_graphlet_report(G1)
report_G2 = get_graphlet_report(G2)

# 3. Crear una tabla comparativa limpia
df_comp = pd.DataFrame({
    'Graphlet': report_G1.keys(),
    'G1': report_G1.values(),
    'G2': report_G2.values()
})

# Filtrar para mostrar solo filas donde haya al menos un graphlet detectado
df_results = df_comp[(df_comp['G1'] > 0) | (df_comp['G2'] > 0)]

print("\n" + "="*35)
print(" COMPARATIVA ESTRUCTURAL M1-M29 ")
print("="*35)
print(df_results.to_string(index=False))
print("="*35)

# 4. Visualización de los grafos
plt.figure(figsize=(10, 4))
plt.subplot(121)
nx.draw(G1, with_labels=True, node_color='skyblue', font_weight='bold')
plt.title("Grafo G1")

plt.subplot(122)
nx.draw(G2, with_labels=True, node_color='orange', font_weight='bold')
plt.title("Grafo G2")

plt.show()