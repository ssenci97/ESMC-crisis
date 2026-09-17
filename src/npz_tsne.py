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
from sklearn.manifold import TSNE

# Override logging.INFO to display as [PYTHON-INFO]
logging.addLevelName(logging.INFO, "[PYTHON-INFO]")

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
log = logging.getLogger()

CONFIG = {"n_components": 2, "random_state": 42}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--npz", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--ncores", type=int, default=16)
    parser.add_argument("--n_components", type=int, default=CONFIG["n_components"])
    args = parser.parse_args()

    log.info("loading %s", args.npz)
    data = np.load(args.npz)
    keys = data.files
    emb_key = [k for k in keys if data[k].ndim == 2][0]
    id_key = [k for k in keys if data[k].ndim == 1][0]
    embeddings = np.asarray(data[emb_key], dtype=np.float32)
    ids = data[id_key]

    log.info("shape %s dtype %s", embeddings.shape, embeddings.dtype)
    log.info("tsne to %d components", args.n_components)

    tsne = TSNE(
        n_components=args.n_components,
        random_state=CONFIG["random_state"],
        n_jobs=args.ncores,
    )
    reduced = tsne.fit_transform(embeddings)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, ids=ids, tsne=reduced)
    log.info("saved %s", args.out)


if __name__ == "__main__":
    main()