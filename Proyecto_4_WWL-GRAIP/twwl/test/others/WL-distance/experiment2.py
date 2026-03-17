import os
import numpy as np
import pandas as pd
from collections import OrderedDict
from pyarrow import feather
from grakel.datasets import fetch_dataset
from experiments.svm import grakel_to_nx
from utils.distances import wl_lower_bound
from tqdm import tqdm

from sklearn.metrics import accuracy_score
from sklearn.model_selection import ParameterGrid, StratifiedKFold
from sklearn.svm import SVC

param_grid = {
    "WLLB": {
        "k": [1,2,3,4],
        "lambda": np.logspace(-4, 1, num=6),
        "C": np.logspace(-3, 3, num=9),
    }
}

def init():
    if "FILENAME" not in os.environ:
        raise KeyError("`FILENAME` is required")

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

cache = dict()

def eval(dataset, train_idx, test_idx, ncv=10, seed=42):
    inp = fetch_dataset(dataset, as_graphs=True, produce_labels_nodes=True)
    G, y = inp.data, inp.target
    y_train = y[train_idx]
    y_test = y[test_idx]
    nx_G = grakel_to_nx(G)
    for k in param_grid["WLLB"]["k"]:
        if ("WLLB", dataset, k) not in cache:
            cache[("WLLB", dataset, k)] = compute_wllb(nx_G, k, verbose=True)

    bar = tqdm(total=len(ParameterGrid(param_grid["WLLB"])) * ncv, desc=f"WLLB, {dataset}")
    best_param = dict()
    best_score = 0.0
    for outer_iter, param in enumerate(ParameterGrid(param_grid["WLLB"])):
        cv_accs = []
        cv = StratifiedKFold(n_splits=ncv, shuffle=True, random_state=seed)
        for inner_iter, (train, valid) in enumerate(cv.split(train_idx, y_train)):
            D = cache[("WLLB", dataset, param["k"])]
            D_tr = D[np.ix_(train_idx[train], train_idx[train])]
            D_val = D[np.ix_(train_idx[valid], train_idx[train])]
            y_tr = y_train[train]
            y_val = y_train[valid]

            K_tr = np.exp(-param["lambda"] * D_tr)
            K_val = np.exp(-param["lambda"] * D_val)
            model = SVC(kernel="precomputed", C=param["C"])
            model.fit(K_tr, y_tr)
            y_pred = model.predict(K_val)
            cv_accs.append(accuracy_score(y_val, y_pred))

            bar.set_postfix(OrderedDict(outer_iter=outer_iter, inner_iter=inner_iter))
            bar.update(1)

        if np.mean(cv_accs) > best_score:
            best_score = np.mean(cv_accs)
            best_param = param

    D = cache[("WLLB", dataset, best_param["k"])]
    K = np.exp(-best_param["lambda"] * D)
    K_train = K[np.ix_(train_idx, train_idx)]
    K_test = K[np.ix_(test_idx, train_idx)]
    model = SVC(kernel="precomputed", C=best_param["C"])
    model.fit(K_train, y_train)
    y_pred = model.predict(K_test)

    return accuracy_score(y_test, y_pred)

def main(indexfile="../../out/idx.arrow", ncv=10, seed=42):
    if not os.path.exists(indexfile):
        raise FileNotFoundError(f"Index file {indexfile} not found")
    index_df = feather.read_feather(indexfile)

    outfile = os.environ["FILENAME"]
    if os.path.exists(outfile):
        df = pd.read_csv(
            outfile,
            dtype={"method": str, "dataset": str, "trial": int, "accuracy": float},
        )
    else:
        df = pd.DataFrame(
            {
                "method": pd.Series(dtype=str),
                "dataset": pd.Series(dtype=str),
                "trial": pd.Series(dtype=int),
                "accuracy": pd.Series(dtype=float),
            }
        )

    for _, row in index_df.iterrows():
        dataset = row.dataset
        trial = row.trial
        train_idx = row.trainIdx - 1
        test_idx = row.testIdx - 1
        if (
            (df["method"] == "WLLB")
            & (df["dataset"] == dataset)
            & (df["trial"] == trial)
        ).any():
            continue

        accuracy = eval(dataset, train_idx, test_idx, ncv=ncv, seed=seed)
        df.loc[len(df)] = {
            "method": "WLLB",
            "dataset": dataset,
            "trial": trial,
            "accuracy": accuracy,
        }
        df.to_csv(outfile, header=True, index=False)

if __name__ == "__main__":
    init()
    main()
