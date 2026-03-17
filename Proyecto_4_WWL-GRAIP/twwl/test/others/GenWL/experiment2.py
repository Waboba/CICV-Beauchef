import os
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score
from sklearn.model_selection import ParameterGrid, StratifiedKFold
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from pyarrow import feather
from tqdm import tqdm
import GlobalVariables
import GraphDataToGraphList as gc
from GraphSet import GraphSet

cache = dict()

def eval(
    dataset, train_idx, test_idx,
    runs=3, h_max=4, cluster_type="barycenter", k="sqrt", cluster_iters=3
):
    if ("R-WL", dataset) in cache:
        graphset = cache[("R-WL", dataset)]
    else:
        dataset_path = os.path.join("./Data", dataset)
        graph_data = gc.graph_data_to_graph_list(dataset_path)
        graphset = GraphSet(graph_data)
        for r in range(runs):
            for i in range(h_max):
                graphset.compute_next_gen_labels_type2(i<h_max-1, cluster_type, k, cluster_iters)
            graphset.flush()
        cache[("R-WL", dataset)] = graphset

    X = []
    y = []
    vec_len = graphset.wl.idxc
    for g in graphset.graph_set:
        vec = np.zeros(vec_len)
        g_lbl = g.get_graph_label()
        for lbl in g.svm_lbls:
            vec[lbl] +=1
        X.append(vec)
        y.append(g_lbl)

    X_train = [X[i] for i in train_idx]
    X_test = [X[i] for i in test_idx]
    y_train = [y[i] for i in train_idx]
    y_test = [y[i] for i in test_idx]

    model = SVC(kernel="linear", max_iter=1000000)
    clf = GridSearchCV(model, {"C": [2**i for i in [-12,-8,-5,-3,-1,1,3,5,8,12]]}, scoring="accuracy", cv=3, n_jobs=GlobalVariables.threads)
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
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
            (df["method"] == "R-WL")
            & (df["dataset"] == dataset)
            & (df["trial"] == trial)
        ).any():
            continue

        accuracy = eval(dataset, train_idx, test_idx)
        df.loc[len(df)] = {
            "method": "R-WL",
            "dataset": dataset,
            "trial": trial,
            "accuracy": accuracy,
        }
        df.to_csv(outfile, header=True, index=False)



def init():
    if "FILENAME" not in os.environ:
        raise KeyError("`FILENAME` is required")

if __name__ == "__main__":
    init()
    main()
