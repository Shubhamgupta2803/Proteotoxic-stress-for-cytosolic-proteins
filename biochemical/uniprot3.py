#!/usr/bin/env python3
import os
import argparse
import subprocess
import pandas as pd
from collections import Counter
from io import StringIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

__author__ = "The source script solely belongs to SHUBHAM GOYAL"
if __author__ != "The source script solely belongs to SHUBHAM GOYAL":
    raise RuntimeError("Source integrity check failed. Author signature missing.")

# ===========================================================
# LOGGING
# ===========================================================
print_buffer = StringIO()

def log_print(*args):
    text = " ".join(str(a) for a in args)
    print(text)
    print_buffer.write(text + "\n")


# ===========================================================
# HYDROPATHY
# ===========================================================
hydropathy_scale = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5,
    'C': 2.5, 'Q': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9,
    'L': 3.8, 'M': 1.9, 'P': -1.6, 'S': -0.8,
    'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

def hydropathy_avg(sequence):
    seq = sequence.upper()
    values = [hydropathy_scale.get(aa) for aa in seq if aa in hydropathy_scale]
    return sum(values) / len(values) if values else float("nan")


# ===========================================================
# ALIPHATIC INDEX
# ===========================================================
def aliphatic_index(sequence):
    seq = sequence.upper()
    L = len(seq)
    if L == 0:
        return float("nan")

    counts = Counter(seq)

    return (
        counts["A"] / L * 100 +
        2.9 * counts["V"] / L * 100 +
        3.9 * (counts["I"] + counts["L"]) / L * 100
    )


# ===========================================================
# PROTPARAM FEATURES (SAFE)
# ===========================================================
def compute_protein_features(sequence):
    seq = str(sequence).replace("*", "").upper()

    # Keep only valid amino acids
    valid = set("ACDEFGHIKLMNPQRSTVWY")
    seq = "".join([c for c in seq if c in valid])

    if len(seq) == 0:
        return pd.Series({
            "instability_index": float("nan"),
            "aliphatic_index": float("nan"),
            "gravy": float("nan")
        })

    try:
        P = ProteinAnalysis(seq)
        return pd.Series({
            "instability_index": P.instability_index(),
            "aliphatic_index": aliphatic_index(seq),
            "gravy": P.gravy()
        })
    except Exception:
        return pd.Series({
            "instability_index": float("nan"),
            "aliphatic_index": float("nan"),
            "gravy": float("nan")
        })


# ===========================================================
# IUPRED
# ===========================================================
def run_iupred(sequence, iupred_path, mode="long"):
    seq = sequence.strip().upper()
    if len(seq) == 0:
        return []

    with open("temp_seq.fasta", "w") as f:
        f.write(">seq\n" + seq)

    result = subprocess.run(
        ["python3", iupred_path, "temp_seq.fasta", mode],
        capture_output=True, text=True
    )

    scores = []
    for line in result.stdout.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        try:
            scores.append(float(line.split()[2]))
        except:
            continue

    return scores


def idr_summary(scores, threshold=0.5):
    if not scores:
        return {"idr_percent": float("nan"), "num_idrs": float("nan")}

    L = len(scores)
    disordered = [i for i, s in enumerate(scores) if s > threshold]

    regions = []
    if disordered:
        start = disordered[0]
        for i in range(1, len(disordered)):
            if disordered[i] != disordered[i - 1] + 1:
                regions.append((start, disordered[i - 1]))
                start = disordered[i]
        regions.append((start, disordered[-1]))

    return {
        "idr_percent": len(disordered) / L * 100,
        "num_idrs": len(regions)
    }


# ===========================================================
# MAIN
# ===========================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--iupred_path", required=True)
    parser.add_argument("--logfile", default="feature_log.txt")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    files = [f for f in os.listdir(args.indir) if f.endswith("_annotated.csv")]
    if not files:
        raise ValueError("No *_annotated.csv files found.")

    log_print(f"\nFound {len(files)} annotated files.\n")

    total_skipped = 0

    for fname in files:
        in_path = os.path.join(args.indir, fname)
        out_path = os.path.join(args.output, fname.replace("_annotated.csv", "_features.csv"))

        log_print("\n===============================")
        log_print(f"PROCESSING: {fname}")
        log_print("===============================\n")

        df = pd.read_csv(in_path)

        df = df.loc[:, ~df.columns.str.contains("^Unnamed")]
        df = df.dropna(how="all")
        df = df.reset_index(drop=True)

        df["sequence"] = df["sequence"].fillna("").astype(str).str.upper()

        total = df.shape[0]

        # Count separately
        bad_u = df["sequence"].str.contains("U", na=False).sum()
        bad_x = df["sequence"].str.contains("X", na=False).sum()
        bad_ux = df["sequence"].str.contains("[UX]", na=False).sum()

        total_skipped += bad_ux

        log_print(f"Total rows: {total}")
        log_print(f"Rows with U: {bad_u}")
        log_print(f"Rows with X: {bad_x}")
        log_print(f"Rows skipped (U or X): {bad_ux}\n")

        df_valid = df[~df["sequence"].str.contains("[UX]", na=False)].copy()
        df_valid = df_valid.reset_index(drop=False)

        # FEATURES
        df_valid["hydropathy_avg"] = df_valid["sequence"].apply(hydropathy_avg)

        pp = df_valid["sequence"].apply(compute_protein_features)
        df_valid[["instability_index", "aliphatic_index", "gravy"]] = pp

        idr = df_valid["sequence"].apply(
            lambda s: idr_summary(run_iupred(s, args.iupred_path))
        )
        df_valid[["idr_percent", "num_idrs"]] = pd.DataFrame(idr.tolist())

        # MERGE BACK
        new_cols = [
            "hydropathy_avg",
            "instability_index",
            "aliphatic_index",
            "gravy",
            "idr_percent",
            "num_idrs"
        ]

  
        df_valid = df_valid.drop(columns=["sequence"], errors="ignore")

        df_valid.to_csv(out_path, index=False)

        log_print(f"✔ Saved: {out_path}")
        log_print(f"✔ Processed: {df_valid.shape[0]}")
        log_print(f"✔ Skipped (U+X): {bad_ux}")

    # SAVE LOG
    log_path = os.path.join(args.output, args.logfile)
    with open(log_path, "w") as f:
        f.write(print_buffer.getvalue())

    log_print("\n===============================")
    log_print("✔ ALL FILES COMPLETED")
    log_print(f"✔ TOTAL SKIPPED (U+X): {total_skipped}")
    log_print("===============================\n")