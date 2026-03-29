# Data

Data files are too large for git (~20 GB total). They are available on the **Iridis shared space**.

## Iridis users — copy data in one command

```bash
# From inside your cloned repo:
cd ~/aging-challenge-2026

cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/train.h5ad          data_prep/output/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/val.h5ad            data_prep/output/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/test.h5ad           data_prep/output/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/combined.h5ad       data_prep/output/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/donor_metadata.csv  data_prep/output/

mkdir -p data_prep/output/pseudobulk
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/pseudobulk/combined_pseudobulk_combined.h5ad         data_prep/output/pseudobulk/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/pseudobulk/combined_pseudobulk_donor_aggregated.h5ad data_prep/output/pseudobulk/

mkdir -p data_prep/output/geneformer
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/geneformer/geneformer_pseudobulk_train.tsv.gz data_prep/output/geneformer/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/geneformer/geneformer_pseudobulk_val.tsv.gz   data_prep/output/geneformer/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/geneformer/geneformer_pseudobulk_test.tsv.gz  data_prep/output/geneformer/
```

## What each notebook needs

| Notebook | Files required |
|----------|---------------|
| `01_anndata_and_pseudobulk` | `combined.h5ad`, `pseudobulk/combined_pseudobulk_combined.h5ad`, `pseudobulk/combined_pseudobulk_donor_aggregated.h5ad` |
| `02_baseline_model` | `pseudobulk/combined_pseudobulk_donor_aggregated.h5ad` |
| `03_evaluation_metrics` | `pseudobulk/combined_pseudobulk_donor_aggregated.h5ad` + a completed training run |
| `04_geneformer_embeddings` | `geneformer/geneformer_pseudobulk_{train,val,test}.tsv.gz` |
| **Model training** | `pseudobulk/combined_pseudobulk_donor_aggregated.h5ad` |
| **Submission** | `train.h5ad`, `val.h5ad`, `test.h5ad` (if building your own features) |

## File sizes

| File | Size |
|------|------|
| `train.h5ad` | 7.8 GB |
| `val.h5ad` | 1.0 GB |
| `test.h5ad` | 1.1 GB |
| `combined.h5ad` | 9.8 GB |
| `pseudobulk/combined_pseudobulk_donor_aggregated.h5ad` | 1.1 GB |
| `pseudobulk/combined_pseudobulk_combined.h5ad` | 1.1 GB |
| `geneformer/*.tsv.gz` | 26 MB total |

## External users

Contact [IfLSAdmin@soton.ac.uk](mailto:IfLSAdmin@soton.ac.uk) for a download link.

## Data Dictionary

### h5ad obs columns (cell level)

| Column | Description |
|--------|-------------|
| `donor_id` | Anonymised integer (1–981) |
| `celltype` | One of: CD4_T, CD8_T, NK, B_cells, monocytes |
| `age` | Donor age in years (**train/val only**; absent from test) |
| `_split` | train / val / test |
| `n_genes_by_counts` | Number of detected genes per cell |
| `pct_counts_mt` | % mitochondrial reads (QC metric) |
| `pool_number` | Sequencing pool ID |

### donor_metadata.csv columns

| Column | Description |
|--------|-------------|
| `donor_id` | Anonymised integer matching h5ad |
| `split` | train / val / test |
| `sex` | male / female |
| `sex_binary` | 0 = female, 1 = male |

Note: age is intentionally excluded from test donors in all public files.

## File Sizes

| File | Size |
|------|------|
| train.h5ad | ~4 GB |
| val.h5ad | ~0.5 GB |
| test.h5ad | ~0.6 GB |
| combined_pseudobulk_donor_aggregated.h5ad | ~1 GB |

## Building Pseudobulk Yourself

If you want to rebuild pseudobulk from scratch (e.g. with different cell types):

```bash
# First combine the splits into one file
python -c "
import scanpy as sc
splits = [sc.read_h5ad('train.h5ad'), sc.read_h5ad('val.h5ad'), sc.read_h5ad('test.h5ad')]
combined = sc.concat(splits, join='inner')
combined.write_h5ad('combined_for_pseudobulk.h5ad')
"

# Then aggregate
python data_prep/h5ad_to_pseudobulk.py combined_for_pseudobulk.h5ad -o my_pseudobulk/
```
