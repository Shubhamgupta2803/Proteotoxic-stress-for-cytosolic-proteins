# Gene Localization Classification Pipeline

A multi-stage pipeline to classify genes as **mitochondrial**, **non-mitochondrial**, **readthrough**, or intermediate categories (`name_only`, `length_only`, `unclassified`), and to assign UniProt accession IDs and protein sequences to each gene.

---

## Input Files

| File | Description |
|------|-------------|
| `all_genes_metadata_without_chx.csv` | Primary gene metadata — must contain `gene`, `cds_length`, and `protein_id` columns |
| `mitocarta.xls` | MitoCarta reference database — must contain `Symbol`, `Synonyms`, `ProteinLength`, and `UniProt` columns |
| `uniprot_proteome.tsv` | UniProt human proteome TSV — must contain `Entry`, `Gene Names`, `Gene Names (synonym)`, `Length`, `Sequence`, and `Localization` columns |

---

## Output File

| File | Description |
|------|-------------|
| `all_genes_location_without_chx.csv` | The primary input file augmented with `localization`, `proteinLength`, `UniProt`, and `sequence` columns |

---

## Pipeline Steps

### Step 1 — Load and Prepare Data

```python
df = pd.read_csv("all_genes_metadata_without_chx.csv")
df["Length"] = np.ceil(df["cds_length"] / 3).astype(int) - 1

df_mito = pd.read_excel("mitocarta.xls")
```

**What it does:**

- Loads the main gene metadata CSV into `df`.
- Derives protein length (`Length`) from the CDS nucleotide length using the formula `⌈cds_length / 3⌉ − 1` (divides by the codon size and subtracts the stop codon).
- Loads the MitoCarta reference database into `df_mito`.
- Standardizes column names and types: gene symbols and UniProt IDs are uppercased; lengths are coerced to numeric, replacing any non-numeric values with `NaN`.

---

### Step 2 — MitoCarta Name/Length Classification

```python
df[["localization", "proteinLength", "UniProt"]] = df.apply(classify, axis=1)
```

**What it does:**

- Builds a lookup dictionary (`name_to_info`) from MitoCarta mapping every gene symbol and synonym to a list of `(ProteinLength, UniProt)` tuples.
- Builds a set of all MitoCarta protein lengths (`mito_lengths`).
- Classifies each gene row with the `classify` function using the following priority order:

| Condition | Label |
|-----------|-------|
| Gene name matches MitoCarta **and** protein length matches | `mitochondrial` |
| Protein length exists in MitoCarta lengths but no name match | `length_only` |
| Gene name matches MitoCarta but lengths differ | `name_only` |
| Neither name nor length matches | `unclassified` |

- When `mitochondrial`, the matched `proteinLength` and `UniProt` are recorded.

---

### Step 3 — UniProt Proteome Localization

```python
uniprot["Localization"] = uniprot["Entry"].apply(
    lambda x: "mitochondrial" if x in mito_uniprot_ids else "non-mitochondrial"
)
```

**What it does:**

- Loads the UniProt proteome TSV (`uniprot`).
- Creates a set of all MitoCarta UniProt IDs (`mito_uniprot_ids`).
- Labels every entry in the UniProt proteome as `"mitochondrial"` or `"non-mitochondrial"` based on membership in the MitoCarta ID set.
- Reports how many MitoCarta UniProt IDs are absent from the proteome file (`missing_ids`).

---

### Step 4 — Synonym Rescue for Missing UniProt IDs

```python
df = df.apply(update_loc, axis=1)
```

**What it does:**

- Inspects the MitoCarta entries whose UniProt IDs were absent from the proteome (`missing_info`).
- Builds a synonym set (`synonym_set`) and a synonym-to-canonical-symbol mapping (`syn_to_symbol`) from those entries.
- For any gene currently labeled `name_only` or `length_only`, checks whether its gene symbol appears in the synonym set. If it does, upgrades the label to `"mitochondrial"` and replaces the gene symbol with the canonical MitoCarta symbol.

---

### Step 5 — UniProt-Based Secondary Classification

```python
df_subset = df_subset.apply(classify, axis=1)
df.update(df_subset.drop(columns=["Uniprot"], errors="ignore"))
```

**What it does:**

- Takes all rows still labeled `name_only`, `length_only`, or `unclassified` as `df_subset`.
- Builds two lookup dictionaries (`mito_lookup`, `non_mito_lookup`) from the UniProt proteome, mapping every gene name and synonym to `(canonical_gene, protein_length, UniProt_entry)` tuples, split by the proteome-assigned localization.
- For each unresolved row, searches the mitochondrial lookup first. If a gene name and protein length match exactly, the row is upgraded to `"mitochondrial"`. If no mitochondrial match, repeats the check against the non-mitochondrial lookup and labels as `"non-mitochondrial"`.
- Merges results back into the main `df`, with UniProt accessions written to the `UniProt` column.

---

### Step 6 — Fill Missing UniProt IDs by Name + Length

```python
df = fill_uniprot_by_name_length(df, uniprot)
```

**What it does:**

- Identifies rows where the `UniProt` column is empty, `"0"`, or `NaN`.
- Builds a `(gene, length) → (UniProt_entry, canonical_gene)` lookup from the full UniProt proteome (including synonyms).
- For each row missing a UniProt ID, looks up the `(gene, Length)` pair in the lookup and fills in the matched accession and canonical gene name.

---

### Step 7 — Assign Sequences from UniProt Proteome

```python
df = assign_sequence_from_uniprot(df, uniprot)
```

**What it does:**

- Builds a `UniProt_Entry → protein_sequence` dictionary from the UniProt proteome file.
- Maps each row's `UniProt` accession to its sequence, storing the result in a new `sequence` column.

---

### Step 8 — Fetch Missing Sequences from UniProt REST API

```python
seq, length = fetch_uniprot_seq(uid)
```

**What it does:**

- Identifies rows that still lack a sequence after Step 7 (empty, `None`, or `NaN`).
- For each such row that has a valid UniProt accession, queries the UniProt REST API (`https://rest.uniprot.org/uniprotkb/{uid}.json`).
- Only assigns the fetched sequence if its length matches the expected protein `Length` in the row, preventing incorrect assignments from isoform mismatches.

---

### Step 9 — Mark Readthrough Genes

```python
df.loc[df["gene"].str.contains("-", na=False), "localization"] = "readthrough"
```

**What it does:**

- Any gene whose name contains a hyphen (`-`) is labeled `"readthrough"`. These represent read-through transcripts (e.g., `AATF-ANKHD1`) rather than standard protein-coding genes.

---

### Step 10 — NCBI Entrez Sequence Fetch for Remaining Unresolved Genes

```python
records = fetch_batch(batch)
```

**What it does:**

- Collects the remaining `unclassified`, `name_only`, and `length_only` rows into `data`.
- Extracts their unique `protein_id` values and queries the NCBI Entrez `protein` database in batches of 200 using the Biopython `Entrez.efetch` API (FASTA format).
- Requires a valid NCBI API key (set in `Entrez.api_key`) to use 10 requests/second instead of the default 3.
- Implements exponential backoff: if a batch fails, it retries up to 3 times with 5 s / 10 s / 15 s delays.

---

### Step 11 — Fetch Gene IDs (GeneID) from NCBI GenBank Records

```python
records = fetch_gene_ids(batch)
```

**What it does:**

- Re-queries the NCBI `protein` database for the same protein IDs in GenBank (`.gb`) format.
- Parses each record's feature table to extract the `GeneID` cross-reference.
- Merges the resulting `protein_id → gene_id` mapping back into `df_seq`.

---

### Step 12 — Map GeneIDs to UniProt via UniProt ID Mapping API

```python
geneid_to_uniprot = run_idmapping(gene_ids)
```

**What it does:**

- Submits the collected GeneIDs to the UniProt ID Mapping REST API (`https://rest.uniprot.org/idmapping`), converting from `GeneID` to `UniProtKB` accessions.
- Polls the job status every 3 seconds (up to 3 minutes) until completion.
- Paginates through all result pages using the `Link: rel="next"` HTTP header.
- When multiple UniProt entries map to the same GeneID, prefers reviewed (Swiss-Prot) entries over unreviewed (TrEMBL) entries.
- Merges the resolved UniProt accessions back into `df_seq`.

---

### Step 13 — Match Sequences by Length and Assign Localization

```python
matched_df = temp[temp["Length"] == temp["seq_length"]]
```

**What it does:**

- Merges `data` (unresolved genes) with `df_seq` (NCBI-fetched sequences) on `protein_id`.
- Computes `seq_length` as the character length of each fetched sequence.
- Keeps only rows where the expected `Length` equals `seq_length` (`matched_df`); discards mismatches (`unmatched_df`).
- Looks up each row's UniProt-based localization from the `uniprot` reference frame and updates `localization` to `"mitochondrial"` or `"non-mitochondrial"` accordingly.
- For rows not resolved via UniProt, falls back to checking the gene symbol against the full MitoCarta gene + synonym set.
- Merges all updates back into the main `df`.

---

### Step 14 — Fuzzy Gene-Name Match Against UniProt for Remaining Unresolved Rows

```python
m = df_subset.reset_index().merge(tmp[...], left_on="gene", right_on="match", ...)
```

**What it does:**

- Targets rows still labeled `unclassified` or `length_only`.
- Explodes the UniProt `Gene Names` and `Gene Names (synonym)` columns to one token per row, creating a flat gene-name lookup table (`tmp`).
- Merges `df_subset` against this table on the gene name.
- Computes the absolute length difference between the row's `Length` and the UniProt entry's `Length`.
- Keeps the single closest match per row (minimum `len_diff`).
- Assigns the UniProt entry and sequence only when `len_diff ≤ 1` amino acid.
- Renames columns (`Localization` → `localization`, `Entry` → `UniProt`, `Sequence` → `sequence`) and merges back into the main `df`.

---

### Step 15 — Restore Original Gene Names and Save

```python
df = df.merge(unedited[["protein_id", "gene"]].rename(...), on="protein_id", how="left")
df.to_csv("all_genes_location_without_chx.csv")
```

**What it does:**

- Re-reads the original input CSV to recover the unmodified `gene` column (some steps above replaced gene symbols with canonical names).
- Merges the original gene names back as `gene_old` alongside the pipeline-updated `gene`.
- Writes the final DataFrame to `all_genes_location_without_chx.csv`.

---

## Localization Labels

| Label | Meaning |
|-------|---------|
| `mitochondrial` | Confirmed mitochondrial protein (by name + length, UniProt lookup, or MitoCarta synonym) |
| `non-mitochondrial` | Confirmed non-mitochondrial protein via UniProt |
| `readthrough` | Gene name contains `-`, indicating a read-through transcript |
| `name_only` | Gene name matches MitoCarta but length does not; may be an isoform |
| `length_only` | Protein length matches a MitoCarta entry but gene name does not |
| `unclassified` | No match by name or length in any reference source |

---

## Dependencies

```
pandas
numpy
openpyxl        # for reading .xls MitoCarta file
biopython       # Entrez / SeqIO (Steps 10–11)
requests        # UniProt REST API (Steps 8, 12)
```

Install with:

```bash
pip install pandas numpy openpyxl biopython requests
```

---

## NCBI API Key

Steps 10 and 11 use the NCBI Entrez API. Replace the placeholder credentials with your own:

```python
Entrez.email   = "your_email@example.com"
Entrez.api_key = "your_ncbi_api_key"
```

You can obtain a free API key at [https://www.ncbi.nlm.nih.gov/account/](https://www.ncbi.nlm.nih.gov/account/).

---

## Pipeline Summary

```
all_genes_metadata_without_chx.csv
        │
        └──▶  localization.ipynb   (Gene Localization Classification Pipeline)

Input CSV + MitoCarta + UniProt Proteome
        │
        ├─ Step 1:  Load & derive protein length
        ├─ Step 2:  MitoCarta name/length classification
        ├─ Step 3:  UniProt proteome localization labels
        ├─ Step 4:  Synonym rescue for missing UniProt IDs
        ├─ Step 5:  UniProt-based secondary classification
        ├─ Step 6:  Fill missing UniProt IDs by name + length
        ├─ Step 7:  Assign sequences from UniProt file
        ├─ Step 8:  Fetch missing sequences via UniProt REST API
        ├─ Step 9:  Mark readthrough genes
        ├─ Step 10: NCBI Entrez batch sequence fetch
        ├─ Step 11: NCBI GenBank GeneID extraction
        ├─ Step 12: GeneID → UniProt ID mapping
        ├─ Step 13: Match sequences by length + assign localization
        ├─ Step 14: Fuzzy gene-name match against UniProt
        └─ Step 15: Restore original gene names → save CSV
```
