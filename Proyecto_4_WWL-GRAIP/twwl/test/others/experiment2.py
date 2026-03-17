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
from sklearn.metrics import accuracy_score
from sklearn.model_selection import ParameterGrid, StratifiedKFold
from sklearn.preprocessing import OneHotEncoder
from sklearn.svm import SVC
from tqdm import tqdm

os.environ["SSL_CERT_FILE"] = certifi.where()
os.environ["ssl_cert_file"] = certifi.where()


param_grid = {
    "TWWL": {
        "H": range(1, 7 + 1),
        "lambda": np.logspace(-4, 1, num=6),
        "C": np.logspace(-3, 3, num=9),
    },
    "WWL": {
        "H": range(1, 7 + 1),
        "gamma": [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 10.0],
        "lambda": np.logspace(-4, 1, num=6),
        "C": np.logspace(-3, 3, num=9),
    },
    "WL": {"H": range(1, 7 + 1), "C": np.logspace(-3, 3, num=9)},
    "WL-OA": {"H": range(1, 7 + 1), "C": np.logspace(-3, 3, num=9)},
    "SP": {"C": np.logspace(-3, 3, num=9)},
    "GL": {"k": [3, 4, 5], "C": np.logspace(-3, 3, num=9)},
    "MLG": {
        "C": np.logspace(-3, 3, num=9),
        "L": [1, 2, 3],
        "heta": [0.01, 0.1, 1.0],
    }
}


def init():
    if "FILENAME" not in os.environ:
        raise KeyError("`FILENAME` is required")
    if "METHOD" not in os.environ:
        raise KeyError("`METHOD` is required")


cache = dict()


def eval(method, dataset, train_idx, test_idx, ncv=10, seed=42):
    inp = fetch_dataset(dataset, produce_labels_nodes=True, verbose=False)
    G, y = inp.data, inp.target
    G = np.array(G)
    y_train = y[train_idx]
    y_test = y[test_idx]

    if method == "MLG": # convert node labels to one-hot
        uniq_labels = np.unique(sum([list(g[1].values()) for g in G], []))
        encoder = OneHotEncoder(sparse_output=False)
        encoder.fit([[label] for label in uniq_labels])
        for g in G:
            for v in g[1]:
                g[1][v] = [float(val) for val in encoder.transform([[g[1][v]]])[0]]

    bar = tqdm(
        total=len(ParameterGrid(param_grid[method])) * ncv,
        desc=f"{method}, {dataset}",
    )
    best_param = dict()
    best_score = 0.0
    for outer_iter, param in enumerate(ParameterGrid(param_grid[method])):
        cv_accs = []
        cv = StratifiedKFold(n_splits=ncv, shuffle=True, random_state=seed)
        for inner_iter, (train, valid) in enumerate(cv.split(train_idx, y_train)):
            if method == "WL":
                kernel = WeisfeilerLehman(n_iter=param["H"], normalize=True)
            elif method == "WL-OA":
                kernel = WeisfeilerLehmanOptimalAssignment(
                    n_iter=param["H"], normalize=True
                )
            elif method == "SP":
                kernel = ShortestPath(normalize=True)
            elif method == "GL":
                kernel = GraphletSampling(
                    normalize=True, k=param["k"], random_state=seed
                )
            elif method == "MLG":
                kernel = MultiscaleLaplacian(
                    normalize=True, L=param["L"],
                    heta=param["heta"], n_samples=100,
                    random_state=seed
                )
            else:
                raise KeyError("unknown method")

            if method == "SP" and (method, dataset) in cache:
                K = cache[(method, dataset)]
            elif method in ["WL", "WL-OA"] and (method, dataset, param["H"]) in cache:
                K = cache[(method, dataset, param["H"])]
            elif method == "GL" and (method, dataset, param["k"]) in cache:
                K = cache[("GL", dataset, param["k"])]
            elif method == "MLG" and (method, dataset, param["L"], param["heta"]) in cache:
                K = cache[("MLG", dataset, param["L"], param["heta"])]
            else:
                K = kernel.fit_transform(G)

                if method == "SP":
                    cache[(method, dataset)] = K
                elif method == "GL":
                    cache[("GL", dataset, param["k"])] = K
                elif method == "MLG":
                    cache[("MLG", dataset, param["L"], param["heta"])] = K
                else:
                    cache[(method, dataset, param["H"])] = K

            K_tr = K[np.ix_(train_idx[train], train_idx[train])]  # type: ignore
            K_val = K[np.ix_(train_idx[valid], train_idx[train])]  # type: ignore
            y_tr = y_train[train]
            y_val = y_train[valid]

            model = SVC(kernel="precomputed", C=param["C"])
            model.fit(K_tr, y_tr)
            y_pred = model.predict(K_val)
            cv_accs.append(accuracy_score(y_val, y_pred))

            bar.set_postfix(OrderedDict(outer_iter=outer_iter, inner_iter=inner_iter))
            bar.update(1)

        if np.mean(cv_accs) > best_score:
            best_score = np.mean(cv_accs)
            best_param = param

    if method == "WL":
        kernel = WeisfeilerLehman(n_iter=best_param["H"], normalize=True)
    elif method == "WL-OA":
        kernel = WeisfeilerLehmanOptimalAssignment(
            n_iter=best_param["H"], normalize=True
        )
    elif method == "SP":
        kernel = ShortestPath(normalize=True)
    elif method == "GL":
        kernel = GraphletSampling(normalize=True, k=best_param["k"], random_state=seed)
    elif method == "MLG":
        kernel = kernel = MultiscaleLaplacian(
            normalize=True, L=best_param["L"],
            heta=best_param["heta"], n_samples=100,
            random_state=seed
        )
    else:
        raise KeyError("unknown method")

    K = (
        cache[(method, dataset, best_param["H"])]
        if method in ["WL", "WL-OA"]
        else cache[(method, dataset, best_param["k"])]
        if method == "GL"
        else cache[(method, dataset, best_param["L"], best_param["heta"])]
        if method == "MLG"
        else cache[(method, dataset)]
    )
    K_train = K[np.ix_(train_idx, train_idx)]
    K_test = K[np.ix_(test_idx, train_idx)]
    model = SVC(kernel="precomputed", C=best_param["C"])
    model.fit(K_train, y_train)
    y_pred = model.predict(K_test)

    return accuracy_score(y_test, y_pred)


def main(indexfile="../out/idx.arrow", ncv=10, seed=42):
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

    method = os.environ["METHOD"]
    for _, row in index_df.iterrows():
        dataset = row.dataset
        trial = row.trial
        train_idx = row.trainIdx - 1
        test_idx = row.testIdx - 1
        if (
            (df["method"] == method)
            & (df["dataset"] == dataset)
            & (df["trial"] == trial)
        ).any():
            continue
        accuracy = eval(method, dataset, train_idx, test_idx, ncv=ncv, seed=seed)
        df.loc[len(df)] = {
            "method": method,
            "dataset": dataset,
            "trial": trial,
            "accuracy": accuracy,
        }
        df.to_csv(outfile, header=True, index=False)


if __name__ == "__main__":
    init()
    main()
