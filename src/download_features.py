#!/usr/bin/env python3
import argparse, re, sys, resource, time, warnings, urllib.parse, urllib.request
from datetime import datetime
from io import StringIO
from pathlib import Path
import pandas as pd

CONFIG_PATH = Path("config/query_configs.tsv")
OUT_DIR = Path("data/UniProtKB")
ALL_FIELDS = ["accession", "id", "uniparc_id", "sequence_version", "version", "date_created", "date_modified", "date_sequence_modified", "gene_names", "gene_primary", "gene_synonym", "gene_oln", "gene_orf", "protein_name", "protein_existence", "reviewed", "fragment", "organism_name", "organism_id", "lineage", "lineage_ids", "virus_hosts", "sequence", "length", "mass", "encoded_in", "ft_act_site", "ft_binding", "ft_carbohyd", "ft_chain", "ft_coiled", "ft_compbias", "ft_conflict", "ft_crosslnk", "ft_disulfid", "ft_dna_bind", "ft_domain", "ft_helix", "ft_init_met", "ft_intramem", "ft_lipid", "ft_mod_res", "ft_motif", "ft_mutagen", "ft_non_cons", "ft_non_std", "ft_non_ter", "ft_peptide", "ft_propep", "ft_region", "ft_repeat", "ft_signal", "ft_site", "ft_strand", "ft_topo_dom", "ft_transit", "ft_transmem", "ft_turn", "ft_unsure", "ft_var_seq", "ft_variant", "ft_zn_fing", "ec", "rhea", "cc_catalytic_activity", "cc_cofactor", "cc_activity_regulation", "kinetics", "ph_dependence", "redox_potential", "temp_dependence", "absorption", "cc_function", "cc_pathway", "cc_miscellaneous", "cc_caution", "cc_sequence_caution", "cc_polymorphism", "cc_rna_editing", "cc_alternative_products", "cc_sc_epred", "cc_interaction", "cc_subunit", "cc_subcellular_location", "cc_developmental_stage", "cc_induction", "cc_tissue_specificity", "cc_allergen", "cc_biotechnology", "cc_disruption_phenotype", "cc_disease", "cc_pharmaceutical", "cc_toxic_dose", "cc_ptm", "cc_domain", "protein_families", "go", "go_id", "go_p", "go_c", "go_f", "xref_proteomes", "structure_3d", "lit_pubmed_id", "keywordid", "keyword", "annotation_score", "comment_count", "feature_count", "tools"]
KEYS_NOSEQ = [col for col in ALL_FIELDS if col != "sequence"]

def fetch_uniprot(query: str, fields: list[str]) -> pd.DataFrame:
    fields_query = fields[0] if len(fields) == 1 else ",".join(fields)
    encoded = urllib.parse.quote(query, safe="")
    url = f"https://rest.uniprot.org/uniprotkb/stream?query={encoded}&format=tsv&fields={fields_query}"
    with urllib.request.urlopen(url) as response:
        tsv = response.read().decode("utf-8")
    return pd.read_csv(StringIO(tsv), sep="\t")

def _build_feature_map() -> dict[str, str]:
    mock = fetch_uniprot("(accession:A0A0R4I9Y1)", fields=ALL_FIELDS)
    return dict(zip(mock.columns, ALL_FIELDS))

def rename_columns(df: pd.DataFrame, feature_map: dict[str, str]) -> pd.DataFrame:
    return df.rename(columns=feature_map)

def save_tsv(df: pd.DataFrame, path: Path) -> None:
    if path.exists():
        warnings.warn(f"Overwriting existing file: {path}", UserWarning)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)

def main() -> None:
    start_time = time.time()
    parser = argparse.ArgumentParser(description="Fetch UniProt features.")
    parser.add_argument("--query", default=None)
    parser.add_argument("--label", default="example")
    parser.add_argument("--fp_feats", type=Path)
    args = parser.parse_args()
    feature_map = _build_feature_map()
    if args.query:
        df_feats = rename_columns(fetch_uniprot(args.query, fields=KEYS_NOSEQ), feature_map)
        save_tsv(df_feats, args.fp_feats or (OUT_DIR / f"{args.label}_feats.tsv"))
    else:
        if not CONFIG_PATH.exists():
            raise FileNotFoundError(f"{CONFIG_PATH} does not exist")
        df_config = pd.read_csv(CONFIG_PATH, sep="\t", dtype=str)
        for _, row in df_config.iterrows():
            pid, label = row["proteome_id"], row["label"]
            df_feats = rename_columns(fetch_uniprot(f"upid:{pid}", fields=KEYS_NOSEQ), feature_map)
            save_tsv(df_feats, OUT_DIR / f"{label}_feats.tsv")
    mem_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"[PYTHON-INFO] Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | Elapsed: {time.time() - start_time:.2f}s | Max Memory: {mem_mb:.2f} MB | Status: Success")

if __name__ == "__main__":
    main()