import os, time
import numpy as np
import pandas as pd
from collections import OrderedDict
from pyarrow import feather
from grakel.datasets import fetch_dataset
from experiments.svm import grakel_to_nx
from utils.distances import wl_lower_bound
from tqdm import tqdm

def compute_wllb(nx_G, k, verbose=False):
    n = len(nx_G)
    ret = np.zeros((n, n))

    bar = tqdm(total=0.5 * n * (n-1), desc=f"[compute_wllb, {k}]", disable=not verbose)
    for i in range(n):
        for j in range(n):
            if i >= j:
                continue
            ret[i, j] = wl_lower_bound(nx_G[i], nx_G[j], k)
            ret[j, i] = ret[i, j]
            bar.update(1)
    return ret

def main(outfile, n_trial=3, seed=42):
    datasets = [
        "MUTAG", "PTC_MR",
        # "ENZYMES", "PROTEINS_full", "NCI1",
        # "DD", "COLLAB",
        # "REDDIT-BINARY", "REDDIT-MULTI-5K", "REDDIT-MULTI-12K", "DBLP_v1", "github_stargazers"
    ]

    if os.path.exists(outfile):
        df = pd.read_csv(outfile, dtype={"method": str, "dataset": str, "trial": int, "time": float})
    else:
        df = pd.DataFrame({
            "method": pd.Series(dtype=str),
            "dataset": pd.Series(dtype=str),
            "trial": pd.Series(dtype=int),
            "time": pd.Series(dtype=float)
        })

    bar = tqdm(total=len(datasets) * n_trial)
    for dataset in datasets:
        inp = fetch_dataset(dataset, as_graphs=True, produce_labels_nodes=True)
        G, y = inp.data, inp.target
        nx_G = grakel_to_nx(G)

        for trial in range(1, n_trial + 1):
            if (
                (df["method"] == "R-WL")
                & (df["dataset"] == dataset)
                & (df["trial"] == trial)
            ).any():
                continue

            st = time.perf_counter()
            compute_wllb(nx_G, 3, verbose=True)
            ed = time.perf_counter()

            df.loc[len(df)] = {
                "method": "R-WL",
                "dataset": dataset,
                "trial": trial,
                "time": ed - st,
            }
            df.to_csv(outfile, header=True, index=False)

            bar.set_postfix(OrderedDict(method="R-WL", dataset=dataset, trial=trial))
            bar.update(1)

if __name__ == "__main__":
    main(os.environ["FILENAME"], n_trial=1)
