#!/usr/bin/env python3
import os
import sys

n = 16
for i, a in enumerate(sys.argv):
    if a == "--ncores" and i + 1 < len(sys.argv):
        n = int(sys.argv[i + 1])
        break
for k in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[k] = str(n)

import argparse
import logging
from pathlib import Path
import numpy as np
from sklearn.manifold import SpectralEmbedding, trustworthiness
import umap

# Override logging.INFO to display as [PYTHON-INFO]
logging.addLevelName(logging.INFO, "[PYTHON-INFO]")

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
log = logging.getLogger()

CONFIG = {
    "n_components": 2,
    "random_state": 42,
    "n_neighbors": 15,
    "min_dist": 0.1,
    "metric": "euclidean",
    "method": "umap",
    "affinity": "nearest_neighbors",
    "gamma": None,
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--npz", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--ncores", type=int, default=16)
    parser.add_argument("--n_components", type=int, default=CONFIG["n_components"])
    parser.add_argument("--n_neighbors", type=int, default=CONFIG["n_neighbors"])
    parser.add_argument("--min_dist", type=float, default=CONFIG["min_dist"])
    parser.add_argument("--metric", default=CONFIG["metric"])
    parser.add_argument("--random_state", type=int, default=CONFIG["random_state"])
    parser.add_argument(
        "--method", choices=["umap", "spectral"], default=CONFIG["method"]
    )
    parser.add_argument(
        "--affinity", choices=["nearest_neighbors", "rbf"], default=CONFIG["affinity"]
    )
    parser.add_argument("--gamma", type=float, default=CONFIG["gamma"])
    args = parser.parse_args()

    log.info("loading %s", args.npz)
    data = np.load(args.npz)
    keys = data.files
    emb_key = [k for k in keys if data[k].ndim == 2][0]
    id_key = [k for k in keys if data[k].ndim == 1][0]
    embeddings = np.asarray(data[emb_key], dtype=np.float32)
    ids = data[id_key]

    log.info("shape %s dtype %s", embeddings.shape, embeddings.dtype)
    log.info(
        "method=%s n_components=%d n_neighbors=%d",
        args.method,
        args.n_components,
        args.n_neighbors,
    )

    if args.method == "umap":
        reducer = umap.UMAP(
            n_components=args.n_components,
            n_neighbors=args.n_neighbors,
            min_dist=args.min_dist,
            metric=args.metric,
            random_state=args.random_state,
            n_jobs=args.ncores,
        )
        reduced = reducer.fit_transform(embeddings)
    else:
        gamma = (
            args.gamma if args.gamma is not None else 1.0 / embeddings.shape[1]
        )
        reducer = SpectralEmbedding(
            n_components=args.n_components,
            affinity=args.affinity,
            n_neighbors=args.n_neighbors,
            gamma=gamma,
            eigen_solver="arpack",
            random_state=args.random_state,
            n_jobs=args.ncores,
        )
        reduced = reducer.fit_transform(embeddings)

    tw = trustworthiness(embeddings, reduced, n_neighbors=args.n_neighbors)
    log.info("trustworthiness: %.4f", tw)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        ids=ids,
        embedding=reduced,
        trustworthiness=np.float32(tw),
        method=args.method,
    )
    log.info("saved %s", args.out)


if __name__ == "__main__":
    main()