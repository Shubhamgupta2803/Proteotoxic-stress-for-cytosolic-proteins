# TOM Interactome Enrichment Analysis (`Data_analysis_1.ipynb`)

Computes per-gene enrichment of ribosome footprints in the TOM complex interactome pull-down relative to the total translatome, across two treatment conditions (with CHX and without CHX). Produces RPKM-normalized enrichment scores, per-codon gene trace plots, and cross-condition correlation figures.

---

## Where It Fits in the Full Pipeline

```
all_genes.ipynb         →  all_genes_metadata_without_chx.csv
localization.ipynb      →  all_genes_location_with/without_chx.csv
script1.py              →  K% CDS positions (ribosome position metric)
        │
        ▼
Data_analysis_1.ipynb   →  RPKM enrichment scores + plots
```

---

## Input Files

| File | Condition | Description |
|------|-----------|-------------|
| `enrichment/Translatome_chx_rep1.h5` | With CHX | Translatome ribosome counts, replicate 1 |
| `enrichment/Translatome_chx_rep2.h5` | With CHX | Translatome ribosome counts, replicate 2 |
| `enrichment/TOM_interactome_chx_rep1.h5` | With CHX | TOM pull-down ribosome counts, replicate 1 |
| `enrichment/TOM_interactome_chx_rep2.h5` | With CHX | TOM pull-down ribosome counts, replicate 2 |
| `enrichment/Translatome_rep1.h5` | Without CHX | Translatome ribosome counts, replicate 1 |
| `enrichment/TOM_interactome_rep1.h5` | Without CHX | TOM pull-down ribosome counts, replicate 1 |
| `all_genes_location_with_chx.csv` | — | Gene localization labels (output of `localization.ipynb`) |
| `all_genes_location_without_chx.csv` | — | Gene localization labels, without CHX condition |

### HDF5 File Structure (per gene)

| Element | Description |
|---------|-------------|
| `ds[0]` | Array of nucleotide positions along the CDS |
| `ds[1]` | Ribosome footprint counts at each position |
| `attrs["cds_length"]` | CDS length in nucleotides |

---

## Output Files

| File | Description |
|------|-------------|
| `with_chx/gene_enrichment_list_threshold_rep*_median.csv` | Per-replicate filtered enrichment table (after read-count threshold) |
| `with_chx/merged_avg_RPKM_median.csv` | Averaged RPKM + Log2 enrichment across replicates, with CHX |
| `without_chx/merged_avg_RPKM_median.csv` | Same for without CHX condition |
| `with_chx/scatter_avg_RPKM_median.png` | Log-log scatter: translatome vs interactome RPKM |
| `log2_correlation.csv` | Per-gene Log2 enrichment from both conditions merged |
| `log2_correlation.png` | Scatter comparing Log2 enrichment: with vs without CHX |

---

## Step-by-Step Breakdown

### Step 1 — Load HDF5 Files and Find Common Genes

```python
common_genes = set.intersection(*[set(h.keys()) for h in handles])
```

Loads all 4 HDF5 files for the with-CHX condition and finds the intersection of gene keys. Only genes present in **every** replicate of both the translatome and interactome are carried forward. This guarantees every downstream comparison has a complete set of counts.

---

### Step 2 — Compute Total Read Counts Per File

```python
p1_Total_Reads = calculate_total_counts("p1", p1, common_genes)
p2_Total_Reads = calculate_total_counts("p2", p2, common_genes)
```

Sums all ribosome counts across the common gene set for each file. These totals are used as the library size denominators for RPM and RPKM normalization in subsequent steps.

For the without-CHX condition, totals are loaded from a cached file if one exists, or calculated from the gene list in `all_genes_location_without_chx.csv` and saved for reuse.

---

### Step 3 — Per-Codon Normalization and Smoothing (`gene_trace`)

```python
codon_idx = math.ceil(nt_position / 3)
Dict_c = {k: 0.05 + (Dict_c.get(k, 0) / Reads * 1e6) for k in range(1, codon_length + 1)}
```

For a single gene and condition:

1. **Convert nucleotide positions to codon positions** by dividing by 3 (ceiling). Counts at the same codon are summed.
2. **Normalize to RPM** (reads per million) using the library total, with a pseudocount of `0.05` added to every codon to avoid zero-division in the enrichment ratio.
3. **Apply a 15-codon sliding window** (7-codon flank on each side) to smooth the signal before plotting.

Returns a list of `(codon_position, smoothed_RPM)` pairs covering the full CDS.

---

### Step 4 — Compute Per-Codon Enrichment (`gene_enrichment`)

```python
y3 = [(y2[i] / y1[i]) if y1[i] != 0 and y2[i] != 0 else np.nan for i in range(len(y1))]
```

Divides the interactome RPM trace by the translatome RPM trace at every codon position. The result is a position-by-position enrichment ratio — values above 1 indicate local ribosome enrichment at TOM relative to the translatome background.

---

### Step 5 — Plot Individual Gene Traces (`plot_gene_enrichment`)

```python
plot_gene_enrichment('LARS2')
```

Produces a two-panel figure for any gene:

- **Top panel**: enrichment ratio (interactome / translatome) along the CDS
- **Bottom panel**: raw RPM traces for translatome (black) and interactome (dark red)

X-axis is codon / amino acid position. Used for quality control and inspecting genes of interest.

---

### Step 6 — Whole-Gene RPKM and Enrichment Score (`process_gene_data`)

```python
p1_RPKM = (count / total_reads / length) * 1e9
log2_ratio = np.log2(p2_RPKM / p1_RPKM)
```

For each gene, collapses the positional data into a single enrichment score:

1. Sums all ribosome counts for the gene in each condition.
2. Computes **RPKM** (Reads Per Kilobase per Million) normalizing by both CDS length and library size.
3. Computes **ratio** = interactome RPKM ÷ translatome RPKM.
4. Computes **Log2(ratio)** — the enrichment score used in all downstream plots.
5. Attaches the gene's localization label from the localization CSV.

---

### Step 7 — Filter by Read Count Threshold

```python
filtered_df = df.loc[(df['p1_count'] > 0) & (df['p2_count'] > threshold_p2)]
```

Removes genes with insufficient read coverage. Default threshold: interactome count > 16, translatome count > 0. Genes below this threshold have unreliable RPKM estimates and are excluded from further analysis. The filtered table is saved per replicate.

---

### Step 8 — Merge Replicates and Average RPKM

```python
merged['avg_p1_RPKM'] = merged[p1_cols].mean(axis=1)
merged['avg_p2_RPKM'] = merged[p2_cols].mean(axis=1)
merged['Log2_ratio']  = np.log2(merged['avg_p2_RPKM'] / merged['avg_p1_RPKM'])
```

Loads all per-replicate filtered CSVs. Keeps only genes common across all replicates. Averages the RPKM values across replicates for both the translatome and interactome, then recomputes the Log2 enrichment ratio on the averaged values. Saves the final merged table.

---

### Step 9 — Scatter Plot: Translatome vs Interactome RPKM

```python
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot([mn, mx], [mn, mx], ...)   # diagonal = no enrichment
```

Log-log scatter plot with one dot per gene:

- **Grey dots** — non-mitochondrial genes
- **Blue dots** — mitochondrial genes

The diagonal line represents equal RPKM in both conditions (no enrichment). Genes above the diagonal are enriched in the interactome. Mitochondrial genes are expected to cluster above the diagonal.

Saved as `scatter_avg_RPKM_median.png` at 600 dpi.

---

### Step 10 — ECDF of Log2 Enrichment

```python
sns.ecdfplot(data=data, x='Log2_ratio', hue='localization', ...)
```

Cumulative distribution plot of Log2 enrichment scores split by localization group (mitochondrial, non-mitochondrial, readthrough). A rightward shift of the mitochondrial curve confirms that mitochondrial genes are systematically more enriched in the TOM pull-down than non-mitochondrial genes.

---

### Step 11 — Cross-Condition Correlation: With CHX vs Without CHX

```python
r, p = pearsonr(merged['without_chx'], merged['with_chx'])
```

Merges the Log2 enrichment tables from both CHX conditions on gene name. Computes the Pearson correlation coefficient between the two conditions. Plots a scatter with:

- **X-axis**: Log2 enrichment without CHX
- **Y-axis**: Log2 enrichment with CHX
- **Color**: mitochondrial (blue), non-mitochondrial (grey), readthrough (black)

High correlation confirms that CHX treatment does not fundamentally alter which genes are enriched at TOM. Saved as `log2_correlation.png` at 300 dpi.

---

## Key Metrics Explained

| Metric | Formula | Meaning |
|--------|---------|---------|
| RPM | `(count / total_reads) × 1e6` | Read depth normalized by library size |
| RPKM | `(count / total_reads / cds_length) × 1e9` | Further normalized by gene length |
| Ratio | `interactome RPKM / translatome RPKM` | Fold enrichment at TOM vs background |
| Log2 ratio | `log2(ratio)` | Enrichment score — positive = enriched, 0 = no change, negative = depleted |

---

## Dependencies

```
h5py
numpy
pandas
matplotlib
seaborn
scipy
scikit-learn
```

Install with:

```bash
pip install h5py numpy pandas matplotlib seaborn scipy scikit-learn
```

---

## Usage

Run all cells sequentially in JupyterLab or Jupyter Notebook. Before running, set the path variables at the top of the notebook:

```python
directory_path = "enrichment/"   # folder containing all .h5 files
Genome_path    = "genome/"       # folder containing localization CSVs
group1         = "p1"            # label for translatome
group2         = "p2"            # label for interactome
```
