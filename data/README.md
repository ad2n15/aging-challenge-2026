# Data

Data files are too large for git (~20 GB total). They are available on the **Iridis shared space**.

**Layout:** place files under **`data/`** at the repository root (same layout as the shared `.../aging-challenge-2026/data` folder). Teaching notebooks read paths like `data/pseudobulk/...`, `data/geneformer/...` (donor-level TSV features), and optionally `data/geneformer_parquet/...` (large cell-level parquet — train/val include `age`; **test** has no `age`). Generated plots and training runs from notebooks go under **`results/`** (gitignored).

### Iridis shared scratch — use `data/`, not `data_prep/output/`

On the server, large files live here:

`/scratch/aazd1f17/shared_space/aging-challenge-2026/data/`

There is **no** `.../data_prep/output/combined.h5ad` on shared space. The separate folder `.../aging-challenge-2026/data_prep/` only contains a copy of `h5ad_to_pseudobulk.py` for reference, not prepared pipeline output. Always point copies and Apptainer binds at **`.../data`** (see main `README.md`).

## Iridis users — copy data in one command

```bash
# From inside your cloned repo:
cd ~/aging-challenge-2026

mkdir -p data/pseudobulk data/geneformer data/geneformer_parquet

cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/train.h5ad          data/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/val.h5ad            data/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/test.h5ad           data/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/combined_public.h5ad       data/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/donor_metadata.csv  data/

cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/pseudobulk/combined_pseudobulk_combined_public.h5ad         data/pseudobulk/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/pseudobulk/combined_pseudobulk_donor_aggregated_public.h5ad data/pseudobulk/

cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/geneformer/geneformer_pseudobulk_train.tsv.gz data/geneformer/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/geneformer/geneformer_pseudobulk_val.tsv.gz   data/geneformer/
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/geneformer/geneformer_pseudobulk_test.tsv.gz  data/geneformer/

# Optional (~6.5 GB): cell-level Geneformer parquets for notebook 04 (test file has no age column)
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/data/geneformer_parquet/*.parquet data/geneformer_parquet/
```

### Avoid copying (Apptainer bind)

To save disk, mount the shared folder **onto** `data/` instead of copying. See **“Use shared data without copying”** in the main `README.md` for the exact `--bind` lines and a quick `srun` check that `/scratch` is visible on compute nodes.

## What each notebook needs

Paths are relative to **`data/`**.

| Notebook | Files required |
|----------|---------------|
| `01_anndata_and_pseudobulk` | `combined_public.h5ad`, `pseudobulk/combined_pseudobulk_combined_public.h5ad`, `pseudobulk/combined_pseudobulk_donor_aggregated_public.h5ad` |
| `02_baseline_model` | `pseudobulk/combined_pseudobulk_donor_aggregated_public.h5ad` |
| `03_evaluation_metrics` | training run outputs under `results/` + optional `test_labels_hidden.csv` in `data/` for real test metrics |
| `04_geneformer_embeddings` | `geneformer/geneformer_pseudobulk_{train,val,test}.tsv.gz`; optional `geneformer_parquet/*.parquet` for raw cell embeddings |
| **Model training (CLI)** | Pass `--input data/pseudobulk/combined_pseudobulk_donor_aggregated_public.h5ad` (defaults in `train_age_model.py` may still point at `data_prep/output/` for internal use). |
| **Submission** | `train.h5ad`, `val.h5ad`, `test.h5ad` (if building your own features) |

## File sizes

| File | Size |
|------|------|
| `train.h5ad` | 7.8 GB |
| `val.h5ad` | 1.0 GB |
| `test.h5ad` | 1.1 GB |
| `combined_public.h5ad` | 9.8 GB |
| `pseudobulk/combined_pseudobulk_donor_aggregated_public.h5ad` | 1.1 GB |
| `pseudobulk/combined_pseudobulk_combined_public.h5ad` | 1.1 GB |
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

### Strip test `age` before publishing (organisers)

Internal pipeline objects may still contain `age` on test rows. To write **new** files with test `age` set to missing (NaN), leaving train/val unchanged:

```bash
python scripts/strip_test_age_h5ad.py --input-dir /path/to/internal --output-dir /path/to/staging
# default output when `data_prep/output/` exists: `data_prep/output_public/`
# default filenames: `combined_public.h5ad`, `pseudobulk/combined_pseudobulk_*_public.h5ad`
```

`scripts/prepare_shared_scratch.sh` copies those `*_public.h5ad` files into shared `data/` **with the same names** (no rename).

### Verify test split has no age (release check)

After `data/` is populated, organisers can confirm that rows with `_split == 'test'` have no finite `age` in the combined and pseudobulk objects:

```bash
python scripts/check_test_age_withheld.py
# optional: python scripts/check_test_age_withheld.py --data-dir /path/to/data
```

Exit code `0` means no test-age leak was detected; `1` means a file is missing, `_split` is absent, or at least one test row has a finite age value.

## File Sizes

| File | Size |
|------|------|
| train.h5ad | ~4 GB |
| val.h5ad | ~0.5 GB |
| test.h5ad | ~0.6 GB |
| combined_pseudobulk_donor_aggregated_public.h5ad | ~1 GB |

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
