# Mitochondrial Protein Import Analysis — Master Pipeline README

This pipeline identifies genes whose mRNAs are selectively translated at the mitochondrial outer membrane (via the TOM complex), characterizes where ribosomes sit along their CDS, classifies their subcellular localization, and compares their structural and biochemical properties against non-mitochondrial genes — across two treatment conditions (with and without cycloheximide, CHX).

It is organized as **10 sequential stages**, each consuming the output of the stage before it. This document describes the full pipeline end to end: every input, every output, and every step in between.

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│  STAGE 1   all_genes.ipynb                                              │
│  Raw HDF5 (6 files) → unique gene metadata CSV                          │
└──────────────────────────────┬────────────────────────────────────────┘
                                 ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STAGE 2   localization.ipynb                                           │
│  Gene metadata + MitoCarta + UniProt → localization labels + sequences  │
└──────────────────────────────┬────────────────────────────────────────┘
                                 ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STAGE 3   script1.py                                                   │
│  Raw HDF5 (12 files, 2 conditions) → K% CDS ribosome position per gene  │
└──────────────────────────────┬────────────────────────────────────────┘
                                 ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STAGE 4   plot.py                                                      │
│  K% CDS positions → frequency distributions + ECDFs                     │
└──────────────────────────────┬────────────────────────────────────────┘
                                 ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STAGE 5   Data_analysis_1.ipynb                                        │
│  Raw HDF5 + localization → RPKM enrichment (interactome vs translatome) │
└──────────────────────────────┬────────────────────────────────────────┘
                                 ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STAGE 6   ted.py                                                       │
│  Merged data → threshold-filtered groups + TED structural domain counts │
└──────────────────────────────┬────────────────────────────────────────┘
                                 ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STAGE 7   domain.py                                                    │
│  TED-annotated groups → domain count boxplots + ECDFs                   │
└──────────────────────────────┬────────────────────────────────────────┘
                                 ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STAGE 8   length.py                                                    │
│  TED-annotated groups → protein length boxplots + ECDFs                 │
└──────────────────────────────┬────────────────────────────────────────┘
                                 ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STAGE 9   uniprot3.py                                                  │
│  TED-annotated groups → biochemical features (hydropathy, IDR, etc.)    │
└──────────────────────────────┬────────────────────────────────────────┘
                                 ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STAGE 10  biochemical_plot.py                                          │
│  Feature-annotated groups → instability/aliphatic index boxplots+ECDFs  │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Master Input/Output Map

| Stage | Script/Notebook | Primary Input | Primary Output |
|-------|-----------------|----------------|-----------------|
| 1 | `all_genes.ipynb` | 6 HDF5 files (interactome + translatome) | `all_genes_metadata_without_chx.csv` |
| 2 | `localization.ipynb` | Stage 1 CSV + MitoCarta + UniProt proteome | `all_genes_location_without_chx.csv` |
| 3 | `script1.py` | 12 HDF5 files (2 conditions × interactome/translatome × 3 reps) | `final_merged_output.csv` |
| 4 | `plot.py` | Stage 3 CSV | Frequency + ECDF plots |
| 5 | `Data_analysis_1.ipynb` | HDF5 files + Stage 2 CSV | `merged_avg_RPKM_median.csv` (per condition) |
| 6 | `ted.py` | Stage 3/5 merged CSV with localization | `{group}_annotated.csv` (TED domain counts) |
| 7 | `domain.py` | Stage 6 annotated CSVs | TED domain count boxplots/ECDFs + PDF |
| 8 | `length.py` | Stage 6 annotated CSVs | Protein length boxplots/ECDFs + PDF |
| 9 | `uniprot3.py` | Stage 6 annotated CSVs | `{group}_features.csv` (biochemical properties) |
| 10 | `biochemical_plot.py` | Stage 9 feature CSVs | Instability/aliphatic index boxplots/ECDFs + PDF |

---

## Stage 1 — Gene Metadata Extraction (`all_genes.ipynb`)

**What it does:** Reads 6 ribosome profiling HDF5 files (3 TOM interactome replicates + 3 translatome replicates) and extracts a deduplicated list of every unique gene observed, with basic metadata.

**Input:**
- `TOM_interactome_rep1/2/3.h5`
- `Translatome_rep1/2/3.h5`

Each file's top-level HDF5 keys are gene names; each gene dataset carries `protein_id` and `cds_length` attributes.

**Process:**
1. Iterate over all 6 files in order (interactome first, then translatome).
2. For each gene, if not already recorded, store `gene`, `protein_id`, and `cds_length`.
3. First-file-wins deduplication — a gene seen in interactome rep1 is never overwritten by a later file.
4. Sort alphabetically by gene name and save.

**Output:**
- `all_genes_metadata_without_chx.csv` — one row per unique gene (16,170 genes in the reference run), columns: `gene`, `protein_id`, `cds_length`.

---

## Stage 2 — Localization Classification (`localization.ipynb`)

**What it does:** Classifies every gene as mitochondrial, non-mitochondrial, readthrough, or intermediate, and assigns each one a UniProt accession and full protein sequence — cascading through 4 reference sources (MitoCarta name+length match, UniProt proteome lookup, NCBI Entrez, and fuzzy matching).

**Input:**
- `all_genes_metadata_without_chx.csv` (Stage 1 output)
- `mitocarta.xls` — MitoCarta reference database
- `uniprot_proteome.tsv` — full UniProt human proteome

**Process (15 steps, condensed):**
1. Derive protein length from CDS length (`⌈cds_length/3⌉ − 1`).
2. Classify by MitoCarta name + length match → `mitochondrial`, `name_only`, `length_only`, or `unclassified`.
3. Label the UniProt proteome itself as mitochondrial/non-mitochondrial via MitoCarta UniProt ID overlap.
4. Rescue genes via MitoCarta synonyms.
5. Resolve remaining unclassified genes via UniProt name+length lookup.
6. Fill missing UniProt IDs by name+length matching.
7. Assign protein sequences from the UniProt proteome.
8. Fetch missing sequences via the UniProt REST API (length-validated).
9. Mark hyphenated gene names as `readthrough`.
10–12. For genes still unresolved, fetch sequences and GeneIDs from NCBI Entrez, then map GeneIDs back to UniProt accessions via the UniProt ID Mapping API.
13. Match fetched sequences by length and assign final localization.
14. Fuzzy gene-name match against UniProt (±1 amino acid tolerance) as a last resort.
15. Restore original gene names alongside canonical names and save.

**Output:**
- `all_genes_location_without_chx.csv` (and a `with_chx` equivalent run separately) — adds `localization`, `proteinLength`, `UniProt`, `sequence`, and `gene_old` columns to the Stage 1 table.

**Localization labels:** `mitochondrial`, `non-mitochondrial`, `readthrough`, `name_only`, `length_only`, `unclassified`.

---

## Stage 3 — K% CDS Position Computation (`script1.py`)

**What it does:** For every gene, computes the position along the CDS (as a % of total length) where K% (default 90%) of all ribosome footprint counts are accounted for, working from the 3′ end inward. Computed separately for 12 HDF5 files spanning 2 conditions (with/without CHX) × 2 arms (interactome/translatome) × 3 replicates.

**Input:**
- 2 JSON config files (one per condition), each pointing to 6 HDF5 files (3 interactome + 3 translatome replicates)

**Process:**
1. Load all 6 HDF5 files per condition; keep only genes common to all 6.
2. Apply an interactome median-count filter (drop low-expression genes).
3. For each gene and replicate, sort CDS positions 3′→5′, cumulate counts, and find the position where the cumulative sum reaches K% of the total. Express as % of CDS length.
4. Run this in parallel across up to 4 CPU cores.
5. Average the K% value across replicates; compute standard deviation.
6. Repeat for the second condition and merge both on gene name.

**Output:**
- `df_with_scatter_input.csv`, `df_without_scatter_input.csv` — per-condition K% tables
- `final_merged_output.csv` — both conditions merged, with per-replicate, average, and SD columns for K% CDS position
- `run_summary.txt` — execution log

**Command:**
```bash
python script1.py --with_chx config_with.json --without_chx config_without.json \
  --output_dir results/ --k_fraction 0.9
```

---

## Stage 4 — K% CDS Distribution Plots (`plot.py`)

**What it does:** Visualizes how K% CDS positions are distributed across all genes, for all four experimental groups (translatome/interactome × with/without CHX).

**Input:**
- `final_merged_output.csv` (Stage 3 output)

**Process:**
1. Bin all genes into 20 fixed 5% windows (0–5%, 5–10%, …, 95–100%) per group.
2. Plot binned counts as a smooth line plot and as a step/histogram plot.
3. Plot ECDFs comparing interactome vs translatome, separately for each CHX condition, with median lines.

**Output:**
- `frequency_table.csv`
- `frequency_plot.svg`, `frequency_step.svg`
- `ecdf_with_chx.svg`, `ecdf_without_chx.svg`

**Command:**
```bash
python plot.py --indir results/final_merged_output.csv --outdir results/plots --k_fraction 0.9
```

---

## Stage 5 — Enrichment Analysis (`Data_analysis_1.ipynb`)

**What it does:** Computes RPKM-normalized enrichment scores quantifying how much each gene's translation is concentrated at the TOM interactome relative to the bulk translatome — both at the per-codon level (for trace plots) and at the whole-gene level (for scatter/correlation analysis).

**Input:**
- HDF5 files: `Translatome_chx_rep1/2.h5`, `TOM_interactome_chx_rep1/2.h5` (with CHX) and equivalents without CHX
- `all_genes_location_with_chx.csv` / `all_genes_location_without_chx.csv` (Stage 2 output)

**Process:**
1. Find genes common across all replicates; compute library-size totals.
2. Per-codon: convert nucleotide → codon position, normalize to RPM with pseudocount, smooth with a 15-codon sliding window, compute interactome/translatome ratio.
3. Plot individual gene traces (e.g. LARS2) showing enrichment and raw RPM along the CDS.
4. Whole-gene: compute RPKM per condition, ratio, and Log2(ratio) as the enrichment score; attach localization label.
5. Filter genes by minimum read count (interactome > 16, translatome > 0).
6. Average RPKM across replicates and recompute Log2 ratio on the averages.
7. Plot log-log scatter (translatome vs interactome RPKM, mitochondrial genes highlighted).
8. Plot ECDF of Log2 enrichment by localization group.
9. Correlate Log2 enrichment between with-CHX and without-CHX conditions (Pearson r).

**Output:**
- `with_chx/merged_avg_RPKM_median.csv`, `without_chx/merged_avg_RPKM_median.csv`
- `scatter_avg_RPKM_median.png`
- `log2_correlation.csv`, `log2_correlation.png`

---

## Stage 6 — TED Domain Annotation (`ted.py`)

**What it does:** Splits genes into groups by K% CDS threshold (e.g. ≥90% vs <90%) and localization (mitochondrial vs non-mitochondrial), then fetches the number of structural domains for each protein from the TED (The Encyclopedia of Domains) API.

**Input:**
- Merged CSV containing `avg_interactome_with_chx`, `avg_translatome_with_chx`, `localization`, `UniProt` columns (combination of Stage 3 and Stage 2 outputs)

**Process:**
1. Parse threshold expressions (e.g. `">=90:both"`, `"<90:interactome"`).
2. Filter genes by the K% threshold on the relevant arm(s); split by localization.
3. Query the TED API in parallel (12 threads) for each unique UniProt ID's domain count.
4. Map counts back and save one CSV per group.

**Output:**
- `I_mito_>=90_annotated.csv`, `I_non_>=90_annotated.csv`, `T_mito_<90_annotated.csv`, etc. — one file per (arm × localization × threshold) combination, each with a new `ted_domain_count` column.

**Command:**
```bash
python ted.py --indir final_merged_output.csv --comparisons ">=90:both" "<90:both" --outdir ted_output/
```

---

## Stage 7 — TED Domain Plots (`domain.py`)

**What it does:** Pairs interactome vs translatome groups at the same threshold and plots `ted_domain_count` as boxplots (with Mann-Whitney U statistics) and ECDFs.

**Input:**
- `*_annotated.csv` files from Stage 6

**Process:**
1. Auto-pair interactome/translatome files by matching threshold and localization.
2. Assign colors by arm (blue=interactome, orange=translatome) and threshold rank (light/medium/dark).
3. Run Mann-Whitney U test; plot boxplot with p-value bracket, mean line, effect size, and n.
4. Plot ECDF of domain counts for both groups.
5. Optionally run custom comparisons via `--compare`.

**Output:**
- `ted_domain_count_I_*_VS_T_*.svg` (boxplot), `_ecdf.svg` (ECDF)
- `ALL_PLOTS.pdf`

**Command:**
```bash
python domain.py --annotated ted_output/ --output domain_plots/ --k_fraction 0.9
```

---

## Stage 8 — Protein Length Plots (`length.py`)

**What it does:** Identical plotting machinery to Stage 7, but plots `Length` (protein sequence length in amino acids) instead of domain count, always on a log Y-axis.

**Input:**
- `*_annotated.csv` files from Stage 6 (same files as Stage 7 — both scripts can run on the same folder)

**Output:**
- `Length_I_*_VS_T_*.svg` (boxplot), `_ecdf.svg` (ECDF)
- `ALL_PLOTS.pdf`

**Command:**
```bash
python length.py --annotated ted_output/ --output length_plots/ --k_fraction 0.9
```

---

## Stage 9 — Biochemical Feature Computation (`uniprot3.py`)

**What it does:** Computes 6 biochemical/biophysical features per protein sequence — hydropathy, instability index, aliphatic index, GRAVY, and intrinsic disorder (% and region count) — for every group from Stage 6.

**Input:**
- `*_annotated.csv` files from Stage 6 (must contain a `sequence` column)
- A locally installed IUPred script (path provided via `--iupred_path`)

**Process:**
1. Load each annotated file; drop sequences containing `U` or `X` (non-standard residues).
2. Compute hydropathy average (Kyte-Doolittle scale).
3. Compute instability index, aliphatic index, and GRAVY via Biopython's `ProteinAnalysis` + a custom aliphatic index formula.
4. Run IUPred per sequence (subprocess call) to predict intrinsic disorder; summarize as % disordered and number of disordered regions.
5. Drop the `sequence` column and save.

**Output:**
- `{group}_features.csv` — one per input group, with `sequence` replaced by `hydropathy_avg`, `instability_index`, `aliphatic_index`, `gravy`, `idr_percent`, `num_idrs`
- `feature_log.txt`

**Command:**
```bash
python uniprot3.py --indir ted_output/ --output features_output/ --iupred_path /path/to/iupred2a.py
```

---

## Stage 10 — Biochemical Feature Plots (`biochemical_plot.py`)

**What it does:** Same plotting machinery as Stage 7/8, applied to `instability_index` and `aliphatic_index` from Stage 9.

**Input:**
- `*_features.csv` files from Stage 9

**Output:**
- `instability_index_I_*_VS_T_*.svg` / `_ecdf.svg`
- `aliphatic_index_I_*_VS_T_*.svg` / `_ecdf.svg`
- `ALL_PLOTS.pdf`

**Command:**
```bash
python biochemical_plot.py --features features_output/ --output biochemical_plots/ --k_fraction 0.9
```

---

## Full Dependency List

```
h5py
pandas
numpy
matplotlib
seaborn
scipy
scikit-learn
statsmodels
pingouin
biopython
tqdm
requests
openpyxl
```

Install everything with:

```bash
pip install h5py pandas numpy matplotlib seaborn scipy scikit-learn statsmodels pingouin biopython tqdm requests openpyxl
```

Additionally required, installed separately:
- **IUPred** (Stage 9) — local script, path passed via `--iupred_path`
- **NCBI API key** (Stage 2) — free, obtained at https://www.ncbi.nlm.nih.gov/account/

---

## End-to-End Execution Order

```bash
# Stage 1
jupyter nbconvert --to notebook --execute all_genes.ipynb

# Stage 2
jupyter nbconvert --to notebook --execute localization.ipynb

# Stage 3
python script1.py --with_chx config_with.json --without_chx config_without.json \
  --output_dir results/ --k_fraction 0.9

# Stage 4
python plot.py --indir results/final_merged_output.csv --outdir results/plots

# Stage 5
jupyter nbconvert --to notebook --execute Data_analysis_1.ipynb

# Stage 6
python ted.py --indir results/final_merged_output.csv --comparisons ">=90:both" "<90:both" \
  --outdir ted_output/

# Stage 7
python domain.py --annotated ted_output/ --output domain_plots/

# Stage 8
python length.py --annotated ted_output/ --output length_plots/

# Stage 9
python uniprot3.py --indir ted_output/ --output features_output/ --iupred_path /path/to/iupred2a.py

# Stage 10
python biochemical_plot.py --features features_output/ --output biochemical_plots/
```

---

## Key Concepts Reference

| Term | Meaning |
|------|---------|
| **K% CDS Position** | The codon position (as % of CDS length) where K% of ribosome footprint counts are accounted for, measured from the 3′ end. High = normal elongation, low = 5′ stalling. |
| **RPKM** | Reads Per Kilobase per Million — normalizes ribosome counts by gene length and library size. |
| **Log2 enrichment** | `log2(interactome RPKM / translatome RPKM)` — positive means selectively captured at TOM. |
| **Localization labels** | `mitochondrial`, `non-mitochondrial`, `readthrough`, `name_only`, `length_only`, `unclassified` |
| **TED domain count** | Number of structural domains annotated for a protein by the TED database. |
| **Instability index** | Predicted in vitro protein stability (Guruprasad et al.); >40 suggests instability. |
| **Aliphatic index** | Relative volume occupied by aliphatic side chains; correlates with thermostability. |
| **GRAVY** | Grand Average of Hydropathy — overall hydrophobicity score. |
| **IDR (Intrinsically Disordered Region)** | Stretch of protein predicted to lack stable 3D structure. |

---

## File Naming Convention (Stages 6–10)

All grouped output files follow this pattern:

```
{arm}_{localization}_{operator}{threshold}_{suffix}.csv

arm:          I (interactome) or T (translatome)
localization: mito or non
operator:     >= or <
threshold:    numeric K% value (e.g. 90)
suffix:       annotated (Stage 6) or features (Stage 9)
```

Example: `I_mito_>=90_annotated.csv` = interactome-enriched, mitochondrial genes with K% CDS ≥ 90%, annotated with TED domain counts.
