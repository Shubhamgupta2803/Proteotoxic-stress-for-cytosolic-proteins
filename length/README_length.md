# Protein Length Distribution Plotter (`length.py`)

Reads the annotated CSVs produced by `ted.py`, pairs interactome vs translatome groups at matching thresholds, and produces publication-ready **boxplots** and **ECDF plots** comparing protein sequence length (`Length`) between groups. All plots are saved as SVGs and collected into a single PDF.

> **Note:** `length.py` is structurally identical to `domain.py` — the only difference is the metric being plotted. `domain.py` plots `ted_domain_count`; `length.py` plots `Length` (protein sequence length in amino acids).

---

## Where It Fits in the Full Pipeline

```
ted.py  →  I_mito_>=90_annotated.csv
            T_mito_>=90_annotated.csv  ...
                        │
                        ▼
                   length.py
                        │
                        ├──  Length_I_mito_>=90_VS_T_mito_>=90.svg
                        ├──  Length_I_mito_>=90_VS_T_mito_>=90_ecdf.svg
                        └──  ALL_PLOTS.pdf
```

---

## Input

Annotated CSVs from `ted.py`, placed in the `--annotated` folder. Files must follow the naming convention:

```
{I|T}_{mito|non}_{op}{threshold}_annotated.csv
```

### Required Column

| Column | Description |
|--------|-------------|
| `Length` | Protein sequence length in amino acids |

All other columns from the `ted.py` output are carried through but only `Length` is plotted.

---

## Output Files

All saved to `--output`:

| File | Description |
|------|-------------|
| `ALL_PLOTS.pdf` | All boxplots and ECDFs collected in one PDF |
| `Length_I_*_VS_T_*_.svg` | Boxplot: interactome vs translatome protein length |
| `Length_I_*_VS_T_*_ecdf.svg` | ECDF: interactome vs translatome protein length distribution |

---

## Usage

```bash
python length.py \
  --annotated  ted_output/ \
  --output     length_plots/ \
  --k_fraction 0.9 \
  --compare    "I_mito_>=90,T_mito_>=90;I_non_<90,T_non_<90"
```

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--annotated` | Yes | — | Folder containing annotated CSVs from `ted.py` |
| `--output` | Yes | — | Folder to save all plots and PDF |
| `--k_fraction` | No | `0.9` | K% fraction — used only in plot title formatting (e.g. L₉₀) |
| `--compare` | No | — | Optional custom cross-group comparisons (semicolon-separated pairs) |

---

## What Gets Plotted

For every matched interactome / translatome pair at the same threshold and localization:

### Boxplot
- One box per group (interactome vs translatome)
- **Y-axis**: protein length in amino acids, on a **log scale** (length distributions are right-skewed)
- **Significance bracket**: Mann-Whitney U p-value displayed above the boxes
- **Dashed line**: group mean
- **Effect size**: rank-biserial correlation `r` shown in the top-left corner
- **n**: sample size shown below each box

### ECDF
- Cumulative distribution of protein lengths for both groups on the same axes
- A rightward shift means that group has systematically longer proteins

---

## Step-by-Step Breakdown

### Step 1 — Discover and Pair Files

```python
files = sorted([f for f in os.listdir(folder) if f.startswith("I_") and f.endswith("_annotated.csv")])
```

Scans the annotated folder for all interactome group files. For each one, resolves the matching translatome file at the same threshold and localization. If both files exist, they are paired for comparison.

---

### Step 2 — Assign Colors by Threshold and Arm

Colors follow a 3-shade system based on which arm (interactome = blue, translatome = orange) and where the threshold sits relative to other thresholds in the dataset:

| Arm | Below threshold (`<`) | Lowest `>=` threshold | Highest `>=` threshold |
|-----|-----------------------|----------------------|------------------------|
| Interactome | Light blue `#84d2f6` | Medium blue `#1f77b4` | Dark blue `#014f86` |
| Translatome | Light orange `#ffbf69` | Medium orange `#ff7f0e` | Dark orange `#cc5803` |

---

### Step 3 — Boxplot with Statistics (`plot_metric`)

```python
mwu = pg.mwu(dfA["Length"], dfB["Length"])
p, r = mwu["p-val"].values[0], mwu["RBC"].values[0]
```

1. Runs a **Mann-Whitney U test** (non-parametric) on the two length distributions.
2. Automatically switches to **log Y-axis** for the `Length` metric (always, since lengths span orders of magnitude across the proteome).
3. Draws a significance bracket with the p-value above the boxes.
4. Adds dashed mean lines, effect size `r`, and per-group sample sizes.

---

### Step 4 — ECDF Plot (`plot_ecdf`)

```python
sns.ecdfplot(data=ecdf_df, x="Length", hue="Group", palette=palette, ...)
```

Plots the cumulative proportion of genes by protein length for both groups. Useful for visually assessing whether one group (e.g. interactome-enriched genes) tends to have longer or shorter proteins than the other.

---

### Step 5 — Optional Custom Comparisons (`--compare`)

```bash
--compare "I_mito_>=90,T_mito_>=90;I_non_<90,T_non_<90"
```

Beyond the automatic same-threshold pairing, arbitrary cross-group comparisons can be specified. Semicolons separate multiple comparisons; commas separate the two groups within each pair. Useful for comparing the same arm at different thresholds (e.g. ribosome-stalled genes `<90` vs normally-elongating genes `>=90`).

---

## How It Differs from `domain.py`

| Feature | `domain.py` | `length.py` |
|---------|-------------|-------------|
| Metric plotted | `ted_domain_count` (TED domains) | `Length` (protein length in aa) |
| Input files | `*_annotated.csv` from `ted.py` | Same — `*_annotated.csv` from `ted.py` |
| Log scale trigger | max > 1500 OR metric == `aa_length` | Always (metric == `Length`) |
| Everything else | Identical | Identical |

Both scripts can be run on the same `ted.py` output folder, producing complementary views: one showing structural complexity (domain count), the other showing raw protein size.

---

## Dependencies

```
pandas
seaborn
matplotlib
pingouin
```

Install with:

```bash
pip install pandas seaborn matplotlib pingouin
```
