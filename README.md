# Proteotoxic-stress-for-cytosolic-proteins
Proteotoxic stress leads to the accumulation of cytosolic proteins onto the outer mitochondrial surface via TOM. 

 STAGE 1   all_genes.ipynb                                             
 Raw HDF5 (6 files) → unique gene metadata CSV                          

                                 ▼

│  STAGE 2   localization.ipynb                                           │
│  Gene metadata + MitoCarta + UniProt → localization labels + sequences  │

                                 ▼

│  STAGE 3   script1.py                                                   │
│  Raw HDF5 (12 files, 2 conditions) → K% CDS ribosome position per gene  │

                                 ▼

│  STAGE 4   plot.py                                                      │
│  K% CDS positions → frequency distributions + ECDFs                     │

                                 ▼

│  STAGE 5   Data_analysis_1.ipynb                                        │
│  Raw HDF5 + localization → RPKM enrichment (interactome vs translatome) │

                                 ▼

│  STAGE 6   ted.py                                                       │
│  Merged data → threshold-filtered groups + TED structural domain counts │

                                 ▼

│  STAGE 7   domain.py                                                    │
│  TED-annotated groups → domain count boxplots + ECDFs                   │

                                 ▼

│  STAGE 8   length.py                                                    │
│  TED-annotated groups → protein length boxplots + ECDFs                 │

                                 ▼

│  STAGE 9   uniprot3.py                                                  │
│  TED-annotated groups → biochemical features (hydropathy, IDR, etc.)    │

                                 ▼

│  STAGE 10  biochemical_plot.py                                          │
│  Feature-annotated groups → instability/aliphatic index boxplots+ECDFs  │

