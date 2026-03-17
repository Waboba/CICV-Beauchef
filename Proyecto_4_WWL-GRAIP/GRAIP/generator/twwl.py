import networkx as nx
import numpy as np
from collections import defaultdict
from tqdm import tqdm

class Twwl:
    def __init__(self, H, label_attr='label', measure_attr='measure'):
        self.H = H
        self.label_attr = label_attr
        self.measure_attr = measure_attr
        self.fitted = False
        self.N = None
        self.parent = None
        self.measures = None
        self._label_to_id = None
        self._seqs = None

    def fit(self, graphs, update_measure=True, verbose=False):
        self.N = len(graphs)
        N = self.N
        H = self.H

        if update_measure:
            for G in graphs:
                measure = 1.0 / G.number_of_nodes()
                for node in G.nodes:
                    G.nodes[node][self.measure_attr] = measure
        else:
            for G in graphs:
                for node in G.nodes:
                    if self.measure_attr not in G.nodes[node]:
                        raise ValueError(f"Node {node} lacks attribute '{self.measure_attr}'")

        for G in graphs:
            for node in G.nodes:
                if self.label_attr not in G.nodes[node]:
                    raise ValueError(f"Node {node} lacks attribute '{self.label_attr}'")

        self._label_to_id = {}
        self.parent = []
        self.measures = {}
        self.parent.append(-1)
        self.measures[0] = np.zeros(N)

        seqs = [np.zeros((G.number_of_nodes(), H+1), dtype=np.uint64) for G in graphs]

        pbar = tqdm(total=(H+1)*N, desc="Fitting TWWL", disable=not verbose)
        for gid, G in enumerate(graphs):
            node_to_idx = {node: i for i, node in enumerate(G.nodes)}
            for node in G.nodes:
                label = G.nodes[node][self.label_attr]
                # Convertir a uint64
                label_hash = hash(label) & 0xffffffffffffffff
                seqs[gid][node_to_idx[node], 0] = label_hash

                if label_hash not in self._label_to_id:
                    new_id = len(self.parent)
                    self._label_to_id[label_hash] = new_id
                    self.parent.append(0)
                    self.measures[new_id] = np.zeros(N)
            pbar.update(1)

        for h in range(H):
            for gid, G in enumerate(graphs):
                node_to_idx = {node: i for i, node in enumerate(G.nodes)}
                for node in G.nodes:
                    u_idx = node_to_idx[node]
                    current_hash = seqs[gid][u_idx, h]
                    current_id = self._label_to_id[current_hash]

                    if h == H-1:
                        self.measures[current_id][gid] += G.nodes[node][self.measure_attr]
                    else:
                        neighbor_hashes = []
                        for nb in G.neighbors(node):
                            nb_idx = node_to_idx[nb]
                            neighbor_hashes.append(seqs[gid][nb_idx, h])
                        neighbor_hashes.sort()
                        combined = (current_hash, tuple(neighbor_hashes))
                        new_hash = hash(combined) & 0xffffffffffffffff
                        seqs[gid][u_idx, h+1] = new_hash

                        if new_hash not in self._label_to_id:
                            new_id = len(self.parent)
                            self._label_to_id[new_hash] = new_id
                            self.parent.append(current_id)
                            self.measures[new_id] = np.zeros(N)
                pbar.update(1)
        pbar.close()

        self.fitted = True
        self._seqs = seqs

    def tree_metric(self, verbose=False):
        if not self.fitted:
            raise RuntimeError("El modelo no ha sido ajustado. Ejecute fit() primero.")

        N = self.N
        dist = np.zeros((N, N))
        M = len(self.parent)
        pbar = tqdm(total=M-1, desc="Computing tree metric", disable=not verbose)

        # Copia de medidas para modificar (excluimos la raíz virtual)
        measures = {nid: self.measures[nid].copy() for nid in range(1, M)}

        for nid in range(M-1, 0, -1):
            measure = measures[nid]
            # Índices donde la medida es no nula y sus valores
            nz_idx = np.nonzero(measure)[0]
            vals = measure[nz_idx]
            if len(nz_idx) == 0:
                continue

            # Contribución a filas
            dist[nz_idx, :] += vals[:, np.newaxis] / 2.0
            dist[nz_idx, nz_idx] -= vals / 2.0

            # Contribución a columnas
            dist[:, nz_idx] += vals[np.newaxis, :] / 2.0
            dist[nz_idx, nz_idx] -= vals / 2.0

            # Contribución de diferencias absolutas entre pares
            diff_mat = np.abs(vals[:, np.newaxis] - vals[np.newaxis, :])
            dist[nz_idx[:, np.newaxis], nz_idx] += diff_mat / 2.0

            # Propagar medida al padre (si no es la raíz)
            pid = self.parent[nid]
            if pid != 0:
                measures[pid] += measure

            pbar.update(1)
        pbar.close()

        # Asegurar simetría (por posibles errores numéricos)
        dist = (dist + dist.T) / 2.0
        return dist