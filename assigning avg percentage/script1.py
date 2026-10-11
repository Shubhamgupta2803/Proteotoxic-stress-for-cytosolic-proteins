#!/usr/bin/env python3
import h5py
import pandas as pd
import numpy as np
import os
import json
import argparse
from tqdm import tqdm
from functools import reduce
from multiprocessing import Pool, cpu_count, freeze_support
import re

def natural_sort_key(s):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', s)]



# ============================================================
# 2. GLOBAL K% CDS CALCULATOR
# ============================================================
def compute_kpct_position(df, count_col, k_fraction):

    total = df[count_col].sum()

    k = int(np.round(k_fraction * total))   # dynamic K%
    k = max(k, 1)

    df_sorted = df.sort_values("cds_position", ascending=False)
    cum = df_sorted[count_col].cumsum()

    posK = df_sorted.loc[cum >= k, "cds_position"].iloc[0]
    cds_len = df_sorted["cds_length"].iloc[0]

    pct = (posK / cds_len) * 100
    return int(posK), float(pct)


# ============================================================
# 3. WORKER FUNCTION FOR MULTIPROCESSING
# ============================================================
def process_gene(args):
    gene, filtered, k_fraction = args

    match = None

    for df in filtered.values():

        gdf = df[df["gene"] == gene]

        if not gdf.empty:
            match = gdf.iloc[0]
            break

    row = {
        "gene": gene,
        "gene_id": match["gene_id"],
        "protein_id": match["protein_id"],
        "tag": match["tag"],
        "cds_length": int(match["cds_length"])
    }

    # %CDS for each replicate
    for short, df in filtered.items():

        col = [c for c in df.columns if c.startswith("count_")][0]

        gdf = df[df["gene"] == gene]

        if gdf.empty:
            row[f"pct_{short}"] = 0

        else:
            _, pct = compute_kpct_position(gdf, col, k_fraction)
            row[f"pct_{short}"] = pct

    return row


# ============================================================
# 4. PIPELINE FOR A SINGLE CONDITION
# ============================================================
def run_from_config(config_path, k_fraction):

    with open(config_path, "r") as f:
        config = json.load(f)

    condition_name = config["condition_name"]
    replicate_files = config["replicate_files"]

    log_print(f"\n=== Processing: {condition_name.upper()} ===")

    all_files = replicate_files["interactome"] + replicate_files["translatome"]

    file_data = {}
    gene_sets = []

    # ----------------------------- LOAD HDF5 -----------------------------
    for fname in all_files:
        short = os.path.basename(fname).replace(".h5", "")
        log_print(f"Loading: {short}")

        rows = []

        with h5py.File(fname, "r") as f:

            for gene in f.keys():

                ds = f[gene]

                cds_len = int(ds.attrs["cds_length"])

                gene_id = ds.attrs.get("gene_id", "")
                protein_id = ds.attrs.get("protein_id", "")

                # decode bytes if stored as bytes
                if isinstance(gene_id, bytes):
                    gene_id = gene_id.decode()

                if isinstance(protein_id, bytes):
                    protein_id = protein_id.decode()
                tag = ds.attrs.get("tag", "")

                if isinstance(tag, bytes):
                    tag = tag.decode()

                positions = ds[0]
                counts = ds[1]

                for pos, cnt in zip(positions, counts):

                    rows.append([
                        gene,
                        gene_id,
                        protein_id,
                        tag,
                        cds_len,
                        int(pos),
                        int(cnt)
                    ])

        df = pd.DataFrame(
            rows,
            columns=[
                "gene",
                "gene_id",
                "protein_id",
                "tag",
                "cds_length",
                "cds_position",
                f"count_{short}"
            ]
        )

        file_data[short] = df
        gene_sets.append(set(df["gene"]))

    # -------------------------- COMMON GENES -----------------------------
    common_genes = set.intersection(*gene_sets)

    filtered = {
        short: df[df["gene"].isin(common_genes)].copy()
        for short, df in file_data.items()
    }

    # ---------------------- MERGE TOTAL COUNTS ---------------------------
    merged_parts = []

    for short, df in filtered.items():

        col = [c for c in df.columns if c.startswith("count_")][0]

        merged_parts.append(
            df.groupby("gene")[col]
            .sum()
            .reset_index()
            .rename(columns={col: short})
        )

    merged = reduce(lambda x, y: pd.merge(x, y, on="gene"), merged_parts)

    # ---------------------- MEDIAN FILTER -------------------------------
    interactome_cols = [
        os.path.basename(f).replace(".h5", "")
        for f in replicate_files["interactome"]
    ]

    medians = {
        col: merged[col].median()
        for col in interactome_cols
    }

    mask = np.logical_and.reduce(
        [merged[col] > medians[col] for col in interactome_cols]
    )

    genes_keep = merged.loc[mask, "gene"].tolist()

    log_print(f"Genes passing interactome median cutoff: {len(genes_keep)}")

    # ---------------------- PARALLEL PROCESSING --------------------------
    N_CORES = min(4, cpu_count())

    log_print(
        f"\nRunning {int(k_fraction*100)}% CDS computation using {N_CORES} CPU cores…\n"
    )

    args_list = [
        (gene, filtered, k_fraction)
        for gene in genes_keep
    ]

    rows_out = []

    with Pool(N_CORES) as p:

        for result in tqdm(
            p.imap_unordered(process_gene, args_list, chunksize=20),
            total=len(args_list),
            desc=f"{condition_name} {int(k_fraction*100)}% CDS",
            ncols=100,
            dynamic_ncols=True,
            mininterval=0.4,
            smoothing=0.1,
            leave=True
        ):
            rows_out.append(result)

    df_result = pd.DataFrame(rows_out)
    # ---------------------- SUMMARY STATS ---------------------------------

    interactome_pct_cols = [
        c for c in df_result.columns
        if c.startswith("pct_") and "interactome" in c
    ]

    translatome_pct_cols = [
        c for c in df_result.columns
        if c.startswith("pct_") and "Translatome" in c
    ]

    # ============================================================
    # CLEAN REPLICATE COLUMN NAMES
    # ============================================================

    for i, col in enumerate(interactome_pct_cols, start=1):

        df_result.rename(
            columns={
                col: f"interactome_rep{i}_{condition_name}"
            },
            inplace=True
        )

    for i, col in enumerate(translatome_pct_cols, start=1):

        df_result.rename(
            columns={
                col: f"translatome_rep{i}_{condition_name}"
            },
            inplace=True
        )

    # refresh lists after renaming
    interactome_pct_cols = [
        c for c in df_result.columns
        if c.startswith("interactome_rep")
    ]

    translatome_pct_cols = [
        c for c in df_result.columns
        if c.startswith("translatome_rep")
    ]

    # ================= INTERACTOME =================

    df_result[f"avg_interactome_{condition_name}"] = df_result[
        interactome_pct_cols
    ].mean(axis=1)

    df_result[f"sd_interactome_{condition_name}"] = df_result[
        interactome_pct_cols
    ].std(axis=1)

    # ================= TRANSLATOME =================

    df_result[f"avg_translatome_{condition_name}"] = df_result[
        translatome_pct_cols
    ].mean(axis=1)

    df_result[f"sd_translatome_{condition_name}"] = df_result[
        translatome_pct_cols
    ].std(axis=1)

    return df_result
# ========================= MAIN =============================
# ============================================================
if __name__ == "__main__":
    from io import StringIO
    freeze_support()

    parser = argparse.ArgumentParser()
    parser.add_argument("--with_chx", required=True)
    parser.add_argument("--without_chx", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--k_fraction", type=float, default=0.9,
                        help="Fraction of CDS cutoff (0.9=90%, 0.8=80%, 0.7=70%)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ================== SUMMARY LOG BUFFER ==================
    print_buffer = StringIO()

    def log_print(*args, **kwargs):
        text = " ".join(str(a) for a in args)
        print(text)                          # print to terminal normally
        print_buffer.write(text + "\n")      # store in buffer

    df_with = run_from_config(args.with_chx, args.k_fraction)

    scatter_with_path = os.path.join(args.output_dir, "df_with_scatter_input.csv")
    df_with.to_csv(scatter_with_path, index=False)
    log_print(f"✔ Saved scatter-ready file → {scatter_with_path}")

    # -------------------------------------------------------
    # RUN FOR WITHOUT CHX
    # -------------------------------------------------------
    df_without = run_from_config(args.without_chx, args.k_fraction)

    scatter_without_path = os.path.join(args.output_dir, "df_without_scatter_input.csv")
    df_without.to_csv(scatter_without_path, index=False)
    log_print(f"✔ Saved scatter-ready file → {scatter_without_path}")

    # -------------------------------------------------------
    # FINAL MERGED OUTPUT (for Script 3 or analysis)
    # -------------------------------------------------------
    # collect all rep columns dynamically
    with_rep_cols = [
        c for c in df_with.columns
        if c.startswith("interactome_rep") or c.startswith("translatome_rep")
    ]

    without_rep_cols = [
        c for c in df_without.columns
        if c.startswith("interactome_rep") or c.startswith("translatome_rep")
    ]

    final = pd.merge(
        df_with[
            [
                "gene",
                "gene_id",
                "protein_id",
                "tag",
                "cds_length",
                "avg_interactome_with_chx",
                "sd_interactome_with_chx",
                "avg_translatome_with_chx",
                "sd_translatome_with_chx",
            ] + with_rep_cols
        ],
        df_without[
            [
                "gene",
                "avg_interactome_without_chx",
                "sd_interactome_without_chx",
                "avg_translatome_without_chx",
                "sd_translatome_without_chx",
            ] + without_rep_cols
        ],
        on="gene",
        how="inner"
    )
    final_path = os.path.join(args.output_dir, "final_merged_output.csv")
    final.to_csv(final_path, index=False)

    log_print(f"\n✔ Saved final merged CSV → {final_path}  (shape = {final.shape})\n")

    # ================== WRITE FINAL SUMMARY ==================
    summary_path = os.path.join(args.output_dir, "run_summary.txt")
    with open(summary_path, "w") as f:
        f.write(print_buffer.getvalue())

    log_print(f"✔ Summary saved → {summary_path}")
    # =========================================================
