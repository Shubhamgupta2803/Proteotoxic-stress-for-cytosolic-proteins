# Biochemical Feature Analysis (`uniprot3.py` + `biochemical_plot.py`)

Two scripts that together compute biochemical and biophysical properties for each protein sequence — hydropathy, instability, aliphatic index, GRAVY, and intrinsic disorder — and then plot comparisons between interactome/translatome groups as boxplots and ECDFs.

---

## Where They Fit in the Full Pipeline

```
ted.py  →  I_mito_>=90_annotated.csv
            T_mito_>=90_annotated.csv  ...
                        │
                        ▼
                  uniprot3.py
                        │
                        ├──  I_mito_>=90_features.csv
                        ├──  T_mito_>=90_features.csv  ...
                        │
                        ▼
              biochemical_plot.py
                        │
                        ├──  instability_index_*.svg / _ecdf.svg
                        ├──  aliphatic_index_*.svg / _ecdf.svg
                        └──  ALL_PLOTS.pdf
```

---

## `uniprot3.py` — Compute Biochemical Features

### What It Does

Reads every `*_annotated.csv` file in a folder (the output of `ted.py`), computes a set of sequence-derived biochemical properties for each protein, drops sequences containing non-standard amino acid codes (`U`, `X`), and saves a `*_features.csv` per input file.

---

### Input

| File pattern | Description |
|--------------|-------------|
| `*_annotated.csv` | Any file from `ted.py`'s output folder — must contain a `sequence` column |

### Required Column

| Column | Description |
|--------|-------------|
| `sequence` | Full protein amino acid sequence |

### Additional Required Argument

| Argument | Description |
|----------|-------------|
| `--iupred_path` | Path to a local IUPred script (`iupred2a.py` or similar) used for intrinsic disorder prediction. Must be installed/obtained separately — not bundled with this script. |

---

### Output

| File | Description |
|------|-------------|
| `{group}_features.csv` | Each input `*_annotated.csv`, with the `sequence` column replaced by 6 new biochemical feature columns |
| `feature_log.txt` | Run log: row counts, skipped sequences, and per-file summaries |

### New Columns Added

| Column | Description | Source |
|--------|-------------|--------|
| `hydropathy_avg` | Mean Kyte-Doolittle hydropathy across the sequence | Custom lookup table |
| `instability_index` | Predicted protein stability (Guruprasad et al. instability index) | Biopython `ProteinAnalysis` |
| `aliphatic_index` | Relative volume occupied by aliphatic side chains (A, V, I, L) | Custom calculation |
| `gravy` | Grand Average of Hydropathy | Biopython `ProteinAnalysis` |
| `idr_percent` | Percentage of residues predicted to be intrinsically disordered | IUPred |
| `num_idrs` | Number of distinct disordered regions | IUPred |

---

### Step-by-Step Breakdown

#### Step 1 — Discover and Load Annotated Files

```python
files = [f for f in os.listdir(args.indir) if f.endswith("_annotated.csv")]
```

Finds every file ending in `_annotated.csv` in the input directory (these are the outputs of `ted.py`). Raises an error if none are found.

---

#### Step 2 — Clean and Filter Invalid Sequences

```python
bad_ux = df["sequence"].str.contains("[UX]", na=False).sum()
df_valid = df[~df["sequence"].str.contains("[UX]", na=False)].copy()
```

Counts and removes any sequence containing **U** (selenocysteine) or **X** (unknown residue), since these break the hydropathy, ProtParam, and IUPred calculations downstream. Reports how many rows contained `U`, how many contained `X`, and how many were dropped overall.

---

#### Step 3 — Hydropathy Average

```python
hydropathy_avg(sequence)
```

Averages the Kyte-Doolittle hydropathy value for every residue in the sequence using a hardcoded lookup table. Higher values indicate a more hydrophobic protein on average.

---

#### Step 4 — Instability Index, Aliphatic Index, GRAVY (`compute_protein_features`)

```python
P = ProteinAnalysis(seq)
P.instability_index()
P.gravy()
aliphatic_index(seq)
```

Strips any stop codons (`*`) and non-standard residues, then uses Biopython's `ProteinAnalysis` to compute:
- **Instability index** — values above 40 suggest the protein is likely unstable in vitro
- **GRAVY** — overall hydrophobicity score (same underlying concept as hydropathy, computed via Biopython's own scale)

The **aliphatic index** is computed manually using the standard formula:

```
AI = X(Ala) + 2.9 × X(Val) + 3.9 × (X(Ile) + X(Leu))
```

where each `X` is the mole percent of that residue. Higher values indicate greater thermal stability.

If a sequence is empty or fails parsing, all three values are returned as `NaN` rather than raising an error, so one bad sequence doesn't halt the batch.

---

#### Step 5 — Intrinsic Disorder Prediction (IUPred)

```python
result = subprocess.run(["python3", iupred_path, "temp_seq.fasta", mode], ...)
```

For each sequence:
1. Writes it to a temporary FASTA file.
2. Calls the external IUPred script as a subprocess (`mode="long"` by default — IUPred's long-disorder prediction mode).
3. Parses the per-residue disorder scores from IUPred's stdout output.

Then `idr_summary()` processes the per-residue scores:
- Marks any residue with a score **> 0.5** as disordered.
- Computes `idr_percent` = percentage of disordered residues.
- Groups consecutive disordered residues into discrete regions and counts them as `num_idrs`.

**Note:** this step calls IUPred once per sequence sequentially (not parallelized), so it is the slowest part of the script for large gene sets.

---

#### Step 6 — Merge, Drop Sequence Column, Save

```python
df_valid = df_valid.drop(columns=["sequence"], errors="ignore")
df_valid.to_csv(out_path, index=False)
```

Attaches all 6 new feature columns to the DataFrame, removes the now-redundant `sequence` column (keeping file sizes manageable), and saves as `{original_name}_features.csv`.

---

### Usage

```bash
python uniprot3.py \
  --indir       ted_output/ \
  --output      features_output/ \
  --iupred_path /path/to/iupred2a.py \
  --logfile     feature_log.txt
```

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--indir` | Yes | — | Folder containing `*_annotated.csv` files from `ted.py` |
| `--output` | Yes | — | Folder to save `*_features.csv` and the log |
| `--iupred_path` | Yes | — | Path to local IUPred script |
| `--logfile` | No | `feature_log.txt` | Log filename within the output folder |

---

## `biochemical_plot.py` — Plot Boxplots and ECDFs

### What It Does

Identical in structure to `domain.py` / `length.py` (see those READMEs for the shared plotting machinery). Reads the `*_features.csv` files from `uniprot3.py`, pairs interactome vs translatome groups at matching thresholds, and plots **instability index** and **aliphatic index** as boxplots and ECDFs.

---

### Input

`*_features.csv` files from `uniprot3.py`, following the naming convention `{I|T}_{mito|non}_{op}{threshold}_features.csv`.

### Metrics Plotted

```python
metrics = ["instability_index", "aliphatic_index"]
```

Both metrics are plotted automatically for every matched group pair. (`gravy`, `hydropathy_avg`, `idr_percent`, and `num_idrs` are computed by `uniprot3.py` but not plotted by this script by default — add them to the `metrics` list to include them.)

---

### Usage

```bash
python biochemical_plot.py \
  --features   features_output/ \
  --output     biochemical_plots/ \
  --k_fraction 0.9 \
  --compare    "I_mito_>=90,T_mito_>=90;I_non_<90,T_non_<90"
```

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--features` | Yes | — | Folder containing `*_features.csv` from `uniprot3.py` |
| `--output` | Yes | — | Folder to save plots and PDF |
| `--k_fraction` | No | `0.9` | K% fraction — used only in plot title formatting |
| `--compare` | No | — | Optional custom cross-group comparisons |

---

### Output Files

| File | Description |
|------|-------------|
| `ALL_PLOTS.pdf` | All boxplots and ECDFs in one PDF |
| `instability_index_I_*_VS_T_*.svg` | Boxplot of instability index per group pair |
| `instability_index_I_*_VS_T_*_ecdf.svg` | ECDF of instability index per group pair |
| `aliphatic_index_I_*_VS_T_*.svg` | Boxplot of aliphatic index per group pair |
| `aliphatic_index_I_*_VS_T_*_ecdf.svg` | ECDF of aliphatic index per group pair |

---

### Plot Behavior (Same Machinery as `domain.py` / `length.py`)

- **Boxplots** include a Mann-Whitney U significance bracket, dashed mean line, effect size `r`, and sample size `n` per group.
- **Colors** follow the same blue (interactome) / orange (translatome) scheme, shaded light → medium → dark by threshold rank.
- **Log Y-axis** is used automatically only if the max value exceeds 1500 (unlikely for instability/aliphatic index, which are typically in the 0–100 range, so plots are linear by default).
- **`--compare`** allows arbitrary group pairings beyond the automatic same-threshold pairing.

For full details on the shared color logic, statistical test, and pairing mechanism, see the `domain.py` / `length.py` READMEs — the plotting code in this script is byte-for-byte identical except for the `metrics` list.

---

## Biological Interpretation

| Metric | High Value Suggests | Low Value Suggests |
|--------|---------------------|---------------------|
| `hydropathy_avg` | Hydrophobic / membrane-associated protein | Hydrophilic / soluble protein |
| `instability_index` (>40) | Likely unstable in vitro | Likely stable |
| `aliphatic_index` | Greater thermostability | Lower thermostability |
| `gravy` | More hydrophobic overall | More hydrophilic overall |
| `idr_percent` | More intrinsically disordered regions | More structured/folded |
| `num_idrs` | Multiple separate disordered segments | Few or no disordered segments |

These properties are commonly examined in the context of mitochondrial protein import, since precursor proteins destined for mitochondrial import often carry distinct biophysical signatures (e.g., amphipathic, positively charged presequences) compared to cytosolic proteins.

---

## Dependencies

### `uniprot3.py`
```
pandas
biopython     # Bio.SeqUtils.ProtParam
```
Plus an external, separately-installed IUPred script (Python 2 or 3 depending on version) referenced via `--iupred_path`.

### `biochemical_plot.py`
```
pandas
seaborn
matplotlib
pingouin
```

Install with:

```bash
pip install pandas biopython seaborn matplotlib pingouin
```

---

## Notes on Script Internals

The `uniprot3.py` script contains a hardcoded string-equality check at import time (`__author__` variable compared against itself) that always evaluates to true and has no effect on functionality — it exists purely as a static assertion and does not gate any actual behavior. It does not need to be modified or removed to run the script.

IUPred is called once per sequence via subprocess rather than in batch, which is the main performance bottleneck in this script. For large gene sets (hundreds to thousands of proteins), expect this step to dominate total runtime.
