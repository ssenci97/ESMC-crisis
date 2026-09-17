#!/usr/bin/env python3
import argparse, resource, time, warnings, urllib.parse, urllib.request
from datetime import datetime
from io import StringIO
from pathlib import Path
import pandas as pd

CONFIG_PATH = Path("config/query_configs.tsv")
OUT_DIR = Path("data/UniProtKB")
DEFAULT_PROTEOME = "UP000005640"
FIELDS = ["accession", "sequence", "length"]

def load_config_proteomes() -> pd.DataFrame:
    if CONFIG_PATH.exists():
        return pd.read_csv(CONFIG_PATH, sep="\t", dtype=str).replace(r"^\s*$", pd.NA, regex=True).dropna()
    return pd.DataFrame(columns=["proteome_id", "label"])

def fetch_uniprot(query: str, fields: list[str]) -> pd.DataFrame:
    url = f"https://rest.uniprot.org/uniprotkb/stream?query={urllib.parse.quote(query, safe='')}&format=tsv&fields={','.join(fields)}"
    with urllib.request.urlopen(url) as resp:
        return pd.read_csv(StringIO(resp.read().decode("utf-8")), sep="\t")

def get_feature_map() -> dict[str, str]:
    mock = fetch_uniprot("(accession:A0A0R4I9Y1)", FIELDS)
    return dict(zip(mock.columns, FIELDS))

def save_tsv(df: pd.DataFrame, path: Path) -> None:
    if path.exists():
        warnings.warn(f"Overwriting existing file: {path}", UserWarning)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)

def save_fasta(df: pd.DataFrame, path: Path) -> None:
    if path.exists():
        warnings.warn(f"Overwriting existing file: {path}", UserWarning)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for _, row in df.iterrows():
            if pd.notna(row.get("accession")) and pd.notna(row.get("sequence")):
                acc, seq = str(row["accession"]).strip(), str(row["sequence"]).strip()
                f.write(f"{acc if acc.startswith('>') else f'>{acc}'}\n{seq}\n")

def main() -> None:
    start_time = time.time()
    parser = argparse.ArgumentParser(description="Fetch UniProt sequences.")
    parser.add_argument("--query", default=None)
    parser.add_argument("--label", default="example")
    parser.add_argument("--fp_seqs", type=Path)
    parser.add_argument("--fp_fasta", type=Path)
    args = parser.parse_args()
    feature_map = get_feature_map()
    processed_count, total_fetched = 0, 0
    if args.query:
        df_sub = fetch_uniprot(args.query, FIELDS).rename(columns=feature_map)
        save_tsv(df_sub, args.fp_seqs or (OUT_DIR / f"{args.label}_seqs.tsv"))
        save_fasta(df_sub, args.fp_fasta or (OUT_DIR / f"{args.label}_seqs.fasta"))
        processed_count, total_fetched = 1, len(df_sub)
    else:
        df_config = load_config_proteomes()
        items = df_config.to_dict("records") if not df_config.empty else [{"proteome_id": DEFAULT_PROTEOME, "label": "example"}]
        for item in items:
            q, lbl = f"upid:{item['proteome_id']}", item["label"]
            df_sub = fetch_uniprot(q, FIELDS).rename(columns=feature_map)
            save_tsv(df_sub, args.fp_seqs or (OUT_DIR / f"{lbl}_seqs.tsv"))
            save_fasta(df_sub, args.fp_fasta or (OUT_DIR / f"{lbl}_seqs.fasta"))
            processed_count += 1
            total_fetched += len(df_sub)
    mem_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"[PYTHON-INFO] Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | Processed Queries: {processed_count} | Fetched Rows: {total_fetched} | Elapsed: {time.time() - start_time:.2f}s | Max Memory: {mem_mb:.2f} MB | Status: Success")

if __name__ == "__main__":
    main()