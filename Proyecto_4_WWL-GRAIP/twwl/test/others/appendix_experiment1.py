import os
import time
from collections import OrderedDict
from multiprocessing import TimeoutError

import certifi
import numpy as np
import pandas as pd
from grakel.datasets import fetch_dataset
from grakel.kernels import (
    GraphletSampling,
    ShortestPath,
    WeisfeilerLehman,
    WeisfeilerLehmanOptimalAssignment,
    MultiscaleLaplacian
)
from pyarrow import feather
from sklearn.preprocessing import OneHotEncoder
from tqdm import tqdm

def main(outfile, n_trial=3, seed=42):
    datasets = [
        "MUTAG", "PTC_MR", "ENZYMES", "PROTEINS_full", "DD", "NCI1", "COLLAB",
        "REDDIT-BINARY", "REDDIT-MULTI-5K", "REDDIT-MULTI-12K", "DBLP_v1", "github_stargazers"
    ]
    methods = ["SP", "WL", "GL", "WL-OA", "MLG"]

    if os.path.exists(outfile):
        df = pd.read_csv(outfile, dtype={"method": str, "dataset": str, "trial": int, "time": float})
    else:
        df = pd.DataFrame({
            "method": pd.Series(dtype=str),
            "dataset": pd.Series(dtype=str),
            "trial": pd.Series(dtype=int),
            "time": pd.Series(dtype=float)
        })

    bar = tqdm(total=len(methods) * len(datasets) * n_trial)
    for method in methods:
        if method == "SP":
            kernel = ShortestPath(normalize=True)
        elif method == "WL":
            kernel = WeisfeilerLehman(n_iter=5, normalize=True)
        elif method == "WL-OA":
            kernel = WeisfeilerLehman(n_iter=5, normalize=True)
        elif method == "GL":
            kernel = GraphletSampling(normalize=True, k=3, random_state=seed)
        elif method == "MLG":
            kernel = MultiscaleLaplacian(normalize=True, L=2, heta=0.01, n_samples=50, random_state=seed)
        else:
            raise KeyError("unknown method")

        for dataset in datasets:
            inp = fetch_dataset(dataset, produce_labels_nodes=True, verbose=False)
            G = inp.data
            G = np.array(G)

            if method == "MLG":
                uniq_labels = np.unique(sum([list(g[1].values()) for g in G], []))
                encoder = OneHotEncoder(sparse_output=False)
                encoder.fit([[label] for label in uniq_labels])
                for g in G:
                    for v in g[1]:
                        g[1][v] = [float(val) for val in encoder.transform([[g[1][v]]])[0]]

            for trial in range(1, n_trial + 1):
                if (
                    (df["method"] == method)
                    & (df["dataset"] == dataset)
                    & (df["trial"] == trial)
                ).any():
                    continue

                st = time.perf_counter()
                K = kernel.fit_transform(G)
                ed = time.perf_counter()

                df.loc[len(df)] = {
                    "method": method,
                    "dataset": dataset,
                    "trial": trial,
                    "time": ed - st,
                }
                df.to_csv(outfile, header=True, index=False)

                bar.set_postfix(OrderedDict(method=method, dataset=dataset, trial=trial))
                bar.update(1)

if __name__ == "__main__":
    main(os.environ["FILENAME"], n_trial=1)
