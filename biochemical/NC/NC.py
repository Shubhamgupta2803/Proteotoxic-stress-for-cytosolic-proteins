#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob
import os
import math
import argparse


# ==========================================================
# ARGUMENTS
# ==========================================================

parser = argparse.ArgumentParser()

parser.add_argument(
    "--indir",
    required=True,
    help="Folder containing *annotated.csv files"
)

parser.add_argument(
    "--outdir",
    required=True,
    help="Output directory"
)

args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# ==========================================================
# FIND FILES
# ==========================================================

files = glob.glob(os.path.join(args.indir, "*annotated.csv"))

print(f"\nFound {len(files)} annotated files\n")

# ==========================================================
# PROCESS FILES
# ==========================================================

for file in files:

    print(f"Processing: {os.path.basename(file)}")

    df = pd.read_csv(file)

    # ------------------------------------------------------
    # ensure numeric
    # ------------------------------------------------------

    df["avg_interactome_with_chx"] = pd.to_numeric(
        df["avg_interactome_with_chx"],
        errors="coerce"
    )

    df["avg_translatome_with_chx"] = pd.to_numeric(
        df["avg_translatome_with_chx"],
        errors="coerce"
    )

    df["Length"] = pd.to_numeric(
        df["Length"],
        errors="coerce"
    )

    # ------------------------------------------------------
    # amino acid counts
    # ------------------------------------------------------

    interactome_counts = np.round(
        df["Length"] *
        df["avg_interactome_with_chx"] / 100
    ).astype("Int64")

    translatome_counts = np.round(
        df["Length"] *
        df["avg_translatome_with_chx"] / 100
    ).astype("Int64")

    # ------------------------------------------------------
    # extract sequences
    # ------------------------------------------------------

    seq_interactome = []
    seq_translatome = []

    for seq, i_count, t_count in zip(
        df["sequence"],
        interactome_counts,
        translatome_counts
    ):

        if pd.isna(seq):

            seq_interactome.append(np.nan)
            seq_translatome.append(np.nan)
            continue

        seq = str(seq)

        # interactome subsequence
        if pd.isna(i_count):
            seq_interactome.append(np.nan)

        else:
            seq_interactome.append(seq[:int(i_count)])

        # translatome subsequence
        if pd.isna(t_count):
            seq_translatome.append(np.nan)

        else:
            seq_translatome.append(seq[:int(t_count)])

    # ------------------------------------------------------
    # add columns
    # ------------------------------------------------------

    df["sequence_interactome"] = seq_interactome

    df["sequence_translatome"] = seq_translatome

    # optional count columns
    df["interactome_aa_count"] = interactome_counts

    df["translatome_aa_count"] = translatome_counts

    # ------------------------------------------------------
    # save
    # ------------------------------------------------------

    outname = os.path.basename(file).replace(
       "_annotated.csv",
       "_NC_sequence.csv"
    )

    outfile = os.path.join(args.outdir, outname)

    df.to_csv(outfile, index=False)

    print(f"Saved → {outfile}")

print("\nDONE\n")