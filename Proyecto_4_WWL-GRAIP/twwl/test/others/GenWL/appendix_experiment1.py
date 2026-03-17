import os
import time
from collections import OrderedDict

import numpy as np
import pandas as pd
from tqdm import tqdm

import GlobalVariables
import GraphDataToGraphList as gc
from GraphSet import GraphSet

def main(outfile, n_trial=3, seed=42, runs=3, h_max=4, cluster_type="barycenter", k="sqrt", cluster_iters=3):
    datasets = [
        "MUTAG", "PTC_MR", "ENZYMES", "PROTEINS_full", "NCI1",
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
        dataset_path = os.path.join("./Data", dataset)
        graph_data = gc.graph_data_to_graph_list(dataset_path)

        for trial in range(1, n_trial + 1):
            if (
                (df["method"] == "R-WL")
                & (df["dataset"] == dataset)
                & (df["trial"] == trial)
            ).any():
                continue

            st = time.perf_counter()
            graphset = GraphSet(graph_data)
            for _ in range(runs):
                for i in range(h_max):
                    graphset.compute_next_gen_labels_type2(i<h_max-1, cluster_type, k, cluster_iters)
                graphset.flush()
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
