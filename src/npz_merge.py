#!/usr/bin/env python3
import argparse
import numpy as np
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Merge multiple NPZ embedding files into one.")
    parser.add_argument("--npzs", required=True, help="Comma-separated list of NPZ files to merge")
    parser.add_argument("--out_npz", required=True, help="Path to save the merged NPZ output")
    args = parser.parse_args()

    file_paths = [p.strip() for p in args.npzs.split(",")]
    
    all_embeddings = []
    all_ids = []

    for fp in file_paths:
        data = np.load(fp)
        keys = data.files
        
        # Automatically detect 2D embedding and 1D ID keys matching your original script
        emb_key = [k for k in keys if data[k].ndim == 2][0]
        id_key = [k for k in keys if data[k].ndim == 1][0]
        
        all_embeddings.append(np.asarray(data[emb_key], dtype=np.float32))
        all_ids.append(data[id_key])

    merged_embeddings = np.concatenate(all_embeddings, axis=0)
    merged_ids = np.concatenate(all_ids, axis=0)

    Path(args.out_npz).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out_npz, embeddings=merged_embeddings, ids=merged_ids)
    print(f"Successfully merged {len(file_paths)} files into {args.out_npz} with shape {merged_embeddings.shape}")

if __name__ == "__main__":
    main()

