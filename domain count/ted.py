#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import requests
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO

# ============================================================
# 0. PARSE COMPARISON EXPRESSIONS
# ============================================================
def parse_comparison(expr):
    expr = expr.strip()

    if ":" in expr:
        comp, mode = expr.split(":")
        mode = mode.lower()
    else:
        comp = expr
        mode = "both"

    comp = comp.replace(" ", "")

    if comp.startswith(">="):
        op = ">="
        val = float(comp[2:])
    elif comp.startswith("<"):
        op = "<"
        val = float(comp[1:])
    else:
        raise ValueError(f"Invalid comparison {comp}")

    if mode not in ("both", "interactome", "translatome"):
        raise ValueError(f"Invalid mode '{mode}'")

    return op, val, mode


def apply_filter(series, op, val):
    if op == ">=":
        return series >= val
    elif op == "<":
        return series < val
    else:
        raise ValueError("Only < and >= allowed.")


# ============================================================
# 1. TED DOMAIN FUNCTIONS
# ============================================================
def get_ted_domains(uniprot_acc, limit=100):
    try:
        skip = 0
        total_count = 0

        while True:
            url = f"https://ted.cathdb.info/api/v1/uniprot/summary/{uniprot_acc}"
            params = {"skip": skip, "limit": limit}
            r = requests.get(url, params=params, timeout=20)

            if r.status_code != 200:
                return 0

            data = r.json()
            count = data.get("count", 0)

            if skip == 0:
                total_count = count

            if skip + limit >= count:
                break

            skip += limit

        return total_count

    except Exception:
        return 0


def fetch_ted_parallel(uniprot_ids, workers=12):
    results = {}

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(get_ted_domains, uid): uid for uid in uniprot_ids}

        for f in tqdm(as_completed(futures), total=len(futures), desc="Fetching TED domains"):
            uid = futures[f]
            results[uid] = f.result()

    return results


# ============================================================
# 2. PROCESS GROUP (UPDATED: NO UNIPROT ANNOTATION)
# ============================================================
def process_group(df, label, outdir):
    df = df.copy()

    # Keep only valid UniProt IDs
    df = df[df["UniProt"].notna()]
    df = df[df["UniProt"].astype(str).str.strip() != ""]

    print(f"\n===============================")
    print(f"GROUP: {label}")
    print(f"Rows: {df.shape[0]}")
    print("===============================")

    # Unique UniProt IDs
    uniprot_ids = df["UniProt"].unique().tolist()
    print("Unique UniProt IDs:", len(uniprot_ids))

    # Fetch TED domains
    ted_counts = fetch_ted_parallel(uniprot_ids)

    # Map back
    df["ted_domain_count"] = df["UniProt"].map(ted_counts)

    # Save
    out_csv = os.path.join(outdir, f"{label}_annotated.csv")
    df.to_csv(out_csv, index=False)

    print(f"✔ Saved → {out_csv} (shape={df.shape})\n")


# ============================================================
# 3. MAIN
# ============================================================
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Comparison filtering + TED domain extraction"
    )
    parser.add_argument("--indir", required=True)
    parser.add_argument("--comparisons", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.indir)

    interact = "avg_interactome_with_chx"
    translat = "avg_translatome_with_chx"

    comparisons = [parse_comparison(c) for c in args.comparisons]

    # ----------------------------------------------------------
    # PROCESS ALL COMPARISONS
    # ----------------------------------------------------------
    for op, thr, mode in comparisons:
        thr_int = int(thr)

        # === INTERACTOME ===
        if mode in ("both", "interactome"):
            df_i = df[apply_filter(df[interact], op, thr)]

            process_group(
                df_i[df_i["localization"] == "mitochondrial"],
                f"I_mito_{op}{thr_int}",
                args.outdir
            )

            process_group(
                df_i[df_i["localization"] == "non-mitochondrial"],
                f"I_non_{op}{thr_int}",
                args.outdir
            )

        # === TRANSLATOME ===
        if mode in ("both", "translatome"):
            df_t = df[apply_filter(df[translat], op, thr)]

            process_group(
                df_t[df_t["localization"] == "mitochondrial"],
                f"T_mito_{op}{thr_int}",
                args.outdir
            )

            process_group(
                df_t[df_t["localization"] == "non-mitochondrial"],
                f"T_non_{op}{thr_int}",
                args.outdir
            )

    print("\n✔ ALL COMPLETED SUCCESSFULLY.\n")