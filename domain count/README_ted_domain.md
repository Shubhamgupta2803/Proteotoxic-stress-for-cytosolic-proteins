# TED Domain Analysis (`ted.py` + `domain.py`)

Two scripts that together fetch structural domain counts from the TED (The Encyclopedia of Domains) database for each gene's protein, split by localization and enrichment threshold, and then plot comparisons between groups as boxplots and ECDFs.

---

## Where They Fit in the Full Pipeline

```
final_merged_output.csv   +   all_genes_location_with_chx.csv
                │
                ▼
            ted.py          →   I_mito_>=90_annotated.csv
                                I_non_<90_annotated.csv
                                T_mito_>=90_annotated.csv  ...
                │
                ▼
           domain.py        →   boxplots + ECDFs + ALL_PLOTS.pdf
```

---

## `ted.py` — Filter, Group, and Fetch TED Domains

### What It Does

Takes the merged gene CSV (output of `script1.py` merged with localization labels), applies user-defined K% CDS thresholds to split genes into groups, and for each group fetches the number of structural domains from the TED API per UniProt accession. Saves one annotated CSV per group.

---

### Input

| File | Description |
|------|-------------|
| `final_merged_output.csv` (or equivalent) | Must contain `avg_interactome_with_chx`, `avg_translatome_with_chx`, `localization`, and `UniProt` columns |

### Required Columns

| Column | Description |
|--------|-------------|
| `avg_interactome_with_chx` | Mean K% CDS position in the interactome |
| `avg_translatome_with_chx` | Mean K% CDS position in the translatome |
| `localization` | `mitochondrial` or `non-mitochondrial` |
| `UniProt` | UniProt accession ID (used to query TED) |

---

### Usage

```bash
python ted.py \
  --indir   final_merged_output.csv \
  --comparisons ">=90:both" "<90:both" ">=80:interactome" \
  --outdir  ted_output/
```

| Argument | Required | Description |
|----------|----------|-------------|
| `--indir` | Yes | Path to the merged input CSV |
| `--comparisons` | Yes | One or more threshold expressions (see format below) |
| `--outdir` | Yes | Folder to write annotated CSVs |

### Comparison Expression Format

Each expression follows the pattern: `<operator><value>:<mode>`

| Part | Options | Example |
|------|---------|---------|
| operator | `>=` or `<` | `>=90` means K% CDS ≥ 90% |
| value | Any number | `90`, `80`, `70` |
| mode | `both`, `interactome`, `translatome` | `:both` applies the filter to both arms |

Examples:
```
">=90:both"         → genes with K% ≥ 90 in BOTH interactome and translatome
"<90:interactome"   → genes with K% < 90 in interactome only
">=80"              → same as ">=80:both" (default mode is both)
```

---

### Output Files

One CSV per group, named using the pattern `{arm}_{localization}_{op}{threshold}_annotated.csv`:

| Example File | What It Contains |
|--------------|-----------------|
| `I_mito_>=90_annotated.csv` | Mitochondrial genes where interactome K% ≥ 90, with TED domain counts |
| `I_non_>=90_annotated.csv` | Non-mitochondrial genes, interactome K% ≥ 90 |
| `T_mito_<90_annotated.csv` | Mitochondrial genes where translatome K% < 90 |
| `T_non_<90_annotated.csv` | Non-mitochondrial genes, translatome K% < 90 |

Each output CSV is the filtered input rows with one new column added:

| Column | Description |
|--------|-------------|
| `ted_domain_count` | Number of structural domains reported by TED for that protein's UniProt ID |

---

### Step-by-Step Breakdown

#### Step 1 — Parse Comparison Expressions

```python
op, val, mode = parse_comparison(">=90:both")
```

Parses each threshold string into an operator (`>=` or `<`), a numeric value, and a mode (`both`, `interactome`, or `translatome`). Raises a clear error if the format is invalid.

---

#### Step 2 — Apply Filter and Split by Localization

```python
df_i = df[df["avg_interactome_with_chx"] >= thr]
process_group(df_i[df_i["localization"] == "mitochondrial"], "I_mito_>=90", outdir)
process_group(df_i[df_i["localization"] == "non-mitochondrial"], "I_non_>=90", outdir)
```

For each comparison, filters the DataFrame by the K% threshold on the relevant column (interactome, translatome, or both), then splits the result into mitochondrial and non-mitochondrial subsets. Produces up to 4 groups per comparison (I_mito, I_non, T_mito, T_non).

---

#### Step 3 — Fetch TED Domain Counts (Parallel)

```python
ted_counts = fetch_ted_parallel(uniprot_ids, workers=12)
```

For each unique UniProt ID in the group, queries the TED API:

```
GET https://ted.cathdb.info/api/v1/uniprot/summary/{uniprot_acc}
```

Returns the total number of structural domains annotated for that protein. Requests run in parallel across 12 threads with a progress bar. If the API returns a non-200 status or times out, the count defaults to 0. Handles paginated responses automatically.

---

#### Step 4 — Map Counts Back and Save

```python
df["ted_domain_count"] = df["UniProt"].map(ted_counts)
df.to_csv(out_csv, index=False)
```

Maps the fetched domain counts back onto the DataFrame by UniProt accession and saves the annotated CSV.

---

---

## `domain.py` — Plot Boxplots and ECDFs

### What It Does

Reads the annotated CSVs produced by `ted.py`, pairs interactome vs translatome groups at the same threshold, and produces publication-ready boxplots and ECDF plots comparing `ted_domain_count` (or any other metric) between groups. All plots are saved as SVGs and collected into a single PDF.

---

### Input

All annotated CSVs from `ted.py`, organized in the `--annotated` folder. Files must follow the naming convention: `{I|T}_{mito|non}_{op}{threshold}_annotated.csv`

---

### Usage

```bash
python domain.py \
  --annotated  ted_output/ \
  --output     domain_plots/ \
  --k_fraction 0.9 \
  --compare    "I_mito_>=90,T_mito_>=90;I_non_<90,T_non_<90"
```

| Argument | Required | Description |
|----------|----------|-------------|
| `--annotated` | Yes | Folder containing annotated CSVs from `ted.py` |
| `--output` | Yes | Folder to save all plots and the PDF |
| `--k_fraction` | No | K% fraction used — controls axis label formatting (default `0.9`) |
| `--compare` | No | Optional custom pair comparisons (semicolon-separated pairs) |

---

### Output Files

| File | Description |
|------|-------------|
| `ALL_PLOTS.pdf` | All plots collected in a single PDF |
| `ted_domain_count_I_*_VS_T_*_ecdf.svg` | ECDF comparison per group pair |
| `ted_domain_count_I_*_VS_T_*.svg` | Boxplot comparison per group pair |

---

### Step-by-Step Breakdown

#### Step 1 — Discover and Pair Files

```python
files = sorted([f for f in os.listdir(folder) if f.startswith("I_") and f.endswith("_annotated.csv")])
```

Scans the annotated folder for all interactome files (`I_*`). For each interactome file found, resolves the matching translatome file (`T_*`) at the same threshold and localization. If both files exist, they form a comparison pair.

---

#### Step 2 — Assign Colors by Threshold

```python
colA = get_color(labelA, threshold_map)
```

Colors are assigned based on the experimental arm (interactome = blue, translatome = orange) and threshold level (lower threshold = lighter shade, higher threshold = darker shade):

| Arm | Below threshold | Lowest `>=` | Highest `>=` |
|-----|----------------|-------------|--------------|
| Interactome | Light blue `#84d2f6` | Medium blue `#1f77b4` | Dark blue `#014f86` |
| Translatome | Light orange `#ffbf69` | Medium orange `#ff7f0e` | Dark orange `#cc5803` |

---

#### Step 3 — Boxplot with Statistics (`plot_metric`)

```python
mwu = pg.mwu(dfA[metric], dfB[metric])
```

For each group pair and metric:

1. Runs a **Mann-Whitney U test** (non-parametric, appropriate for skewed domain count distributions) using `pingouin`.
2. Plots a **boxplot** with outliers shown as open circles.
3. Annotates the plot with:
   - **p-value** in a significance bracket above the boxes
   - **Mean** as a dashed horizontal line inside each box
   - **Effect size r** (rank-biserial correlation) in the top-left corner
   - **n** (sample size) below each box
4. Uses log scale on the Y-axis when values span a large range (max > 1500) or for sequence length.

---

#### Step 4 — ECDF Plot (`plot_ecdf`)

```python
sns.ecdfplot(data=ecdf_df, x=metric, hue="Group", palette=palette, ...)
```

Plots the cumulative distribution of the metric for both groups on the same axes. A rightward shift of one curve indicates that group tends to have higher domain counts. Saved as SVG alongside the boxplot.

---

#### Step 5 — Custom Comparisons (`--compare`)

```bash
--compare "I_mito_>=90,T_mito_>=90;I_non_<90,T_non_<90"
```

Allows arbitrary cross-group comparisons beyond the automatic interactome-vs-translatome pairing. Semicolons separate multiple comparisons; commas separate the two groups within each comparison. Useful for comparing the same arm at different thresholds (e.g. `>=90` vs `<90`).

---

## Color Reference

| Group | Shade | Hex |
|-------|-------|-----|
| Interactome (below threshold) | Light blue | `#84d2f6` |
| Interactome (low `>=` threshold) | Medium blue | `#1f77b4` |
| Interactome (high `>=` threshold) | Dark blue | `#014f86` |
| Translatome (below threshold) | Light orange | `#ffbf69` |
| Translatome (low `>=` threshold) | Medium orange | `#ff7f0e` |
| Translatome (high `>=` threshold) | Dark orange | `#cc5803` |

---

## Dependencies

```
pandas
requests
tqdm
statsmodels
seaborn
matplotlib
pingouin
```

Install with:

```bash
pip install pandas requests tqdm statsmodels seaborn matplotlib pingouin
```

---

## Key Design Decisions

**Parallel TED fetching.** UniProt accessions within a group are fetched concurrently (12 threads) since each API call is network-bound. This reduces total fetch time roughly 12-fold compared to sequential requests.

**Mann-Whitney U, not t-test.** Domain count distributions are non-normal and right-skewed. The Mann-Whitney U test makes no normality assumption. Effect size is reported as rank-biserial correlation `r`, which is directly interpretable (0 = no effect, ±1 = complete separation).

**Log Y-axis when values are large.** The boxplot function automatically switches to a log scale if the maximum value exceeds 1,500, preventing small-count proteins from being compressed into an invisible band at the bottom.

**SVG + PDF.** Individual plots are saved as SVGs for vector editing. All plots are simultaneously written into a single `ALL_PLOTS.pdf` for easy review.
