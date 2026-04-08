# Data

Data files are too large for git (~20 GB total). They are available on the **Iridis shared space**.

**Layout:** mirror the shared scratch tree under **`data/`**. Top-level subfolders (created by `scripts/prepare_shared_scratch.sh`):

| Directory | Contents |
|-----------|----------|
| **`scRNA-seq_raw/`** | `train.h5ad`, `val.h5ad`, `test.h5ad` (cell × gene counts) |
| **`scRNA-seq_pseudobulk/`** | Split donor-aggregated pseudobulk `*_public.h5ad` (test ages stripped) |
| **`scRNA-seq_geneformer/`** | Cell-level Geneformer parquets (train/val include `age`; **test** has no `age`) |
| **`scRNA-seq_geneformer_pseudobulk/`** | `geneformer_pseudobulk_{train,val,test}.tsv.gz` |
| **`genotypes/`** | Competition VCF (or `.vcf.gz`), `Onek1k_competition_gt.tsv`, `pca_{train,val,test}.tsv` |
| **`metadata/`** | `donor_metadata.csv`, optional `splits_info.csv` |

There is **no** single combined cell-level `combined_public.h5ad` or long-form combined pseudobulk in the published bundle. Notebook outputs go under **`results/`** (gitignored).

### Iridis shared scratch — use `data/`, not `data_prep/output/`

On the server, large files live here:

`/scratch/aazd1f17/shared_space/aging-challenge-2026/data/`

There is **no** `.../data_prep/output/combined.h5ad` on shared space. The separate folder `.../aging-challenge-2026/data_prep/` only contains a copy of `h5ad_to_pseudobulk.py` for reference, not prepared pipeline output. Always point copies and Apptainer binds at **`.../data`** (see main `README.md`).

## Iridis users — copy data in one command

```bash
# From inside your cloned repo:
cd ~/aging-challenge-2026

SCR=/scratch/aazd1f17/shared_space/aging-challenge-2026/data

mkdir -p data/scRNA-seq_raw data/scRNA-seq_pseudobulk data/scRNA-seq_geneformer \
         data/scRNA-seq_geneformer_pseudobulk data/genotypes data/metadata

cp -a "${SCR}/scRNA-seq_raw/."                       data/scRNA-seq_raw/
cp -a "${SCR}/scRNA-seq_pseudobulk/."                data/scRNA-seq_pseudobulk/
cp -a "${SCR}/scRNA-seq_geneformer_pseudobulk/."     data/scRNA-seq_geneformer_pseudobulk/
cp -a "${SCR}/genotypes/."                           data/genotypes/
cp -a "${SCR}/metadata/."                            data/metadata/

# Optional (~6.5 GB): cell-level Geneformer parquets (notebook 04)
cp "${SCR}/scRNA-seq_geneformer/"*.parquet data/scRNA-seq_geneformer/ 2>/dev/null || true
```

### Avoid copying (Apptainer bind)

To save disk, mount the shared folder **onto** `data/` instead of copying. See **“Use shared data without copying”** in the main `README.md` for the exact `--bind` lines and a quick `srun` check that `/scratch` is visible on compute nodes.

## What each notebook needs

Paths are relative to **`data/`**.

| Notebook | Files required |
|----------|---------------|
| `00_biology_genome_and_ngs_primer` | *None* (conceptual background only) |
| `01_anndata_and_pseudobulk` | `scRNA-seq_raw/{train,val,test}.h5ad`; `scRNA-seq_pseudobulk/{train,val,test}_pseudobulk_donor_aggregated_public.h5ad` |
| `02_baseline_model` | `scRNA-seq_pseudobulk/train_pseudobulk_donor_aggregated_public.h5ad` (optionally merge train+val locally) |
| `03_evaluation_metrics` | `results/` + optional `test_labels_hidden.csv` (not on shared scratch) |
| `04_geneformer_embeddings` | `scRNA-seq_geneformer_pseudobulk/geneformer_pseudobulk_*.tsv.gz`; optional `scRNA-seq_geneformer/*.parquet` |
| **Model training (CLI)** | e.g. `--input data/scRNA-seq_pseudobulk/train_pseudobulk_donor_aggregated_public.h5ad` |
| **Submission** | `scRNA-seq_raw/{train,val,test}.h5ad` if building features from raw cells |

## File sizes

| File | Size |
|------|------|
| `scRNA-seq_raw/train.h5ad` (etc.) | ~8 / ~1 / ~1 GB |
| `scRNA-seq_pseudobulk/*_public.h5ad` | ~1.1 GB total across three files |
| `scRNA-seq_geneformer_pseudobulk/*.tsv.gz` | ~26 MB total |
| `genotypes/` (VCF + TSV + PCs) | modest (VCF size depends on sites) |

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

### metadata/donor_metadata.csv columns

| Column | Description |
|--------|-------------|
| `donor_id` | Anonymised integer matching h5ad |
| `split` | train / val / test |
| `age` | Donor age in years (**train/val only**; **empty for test** in published `data/`) |
| `sex` | male / female |
| `sex_binary` | 0 = female, 1 = male |
| `self_reported_ethnicity` | Where present in the release |

The internal pipeline CSV can include ages for test donors; `prepare_shared_scratch.sh` runs `scripts/strip_test_age_donor_metadata.py` so the published file does not.

### Strip test `age` before publishing (organisers)

Internal pipeline objects may still contain `age` on test rows. To write **new** files with test `age` set to missing (NaN), leaving train/val unchanged:

```bash
python scripts/split_pseudobulk_donor_aggregated.py   # if needed: split donor-aggregated pseudobulk
python scripts/strip_test_age_h5ad.py --input-dir /path/to/internal --output-dir /path/to/staging
# default output when `data_prep/output/` exists: `data_prep/output_public/`
# typical publish outputs: `pseudobulk/{train,val,test}_pseudobulk_donor_aggregated_public.h5ad` (no combined cell-level h5ad required)
```

`scripts/prepare_shared_scratch.sh` writes the six directories above: raw h5ad → **`scRNA-seq_raw/`**; stripped split pseudobulk → **`scRNA-seq_pseudobulk/`**; parquets → **`scRNA-seq_geneformer/`**; Geneformer TSVs → **`scRNA-seq_geneformer_pseudobulk/`**; VCF/GT/PCs → **`genotypes/`**; `donor_metadata.csv` (+ optional `splits_info.csv`) → **`metadata/`** (test `age` cleared in metadata).

### Verify test split has no age (release check)

After `data/` is populated, organisers can confirm that rows with `_split == 'test'` have no finite `age` in the published h5ad objects:

```bash
python scripts/check_test_age_withheld.py
# optional: python scripts/check_test_age_withheld.py --data-dir /path/to/data
```

Exit code `0` means no test-age leak was detected; `1` means a file is missing, `_split` is absent, or at least one test row has a finite age value.

## File Sizes

| File | Size |
|------|------|
| `scRNA-seq_raw/*.h5ad` | similar to table above |
| `scRNA-seq_pseudobulk/*_public.h5ad` | ~0.2–0.5 GB each |

## Building Pseudobulk Yourself

If you want to rebuild pseudobulk from scratch (e.g. with different cell types):

```bash
# First combine the splits into one file (adjust paths if your data live under scRNA-seq_raw/)
python -c "
import scanpy as sc
base = 'scRNA-seq_raw'
splits = [sc.read_h5ad(f'{base}/train.h5ad'), sc.read_h5ad(f'{base}/val.h5ad'), sc.read_h5ad(f'{base}/test.h5ad')]
combined = sc.concat(splits, join='inner')
combined.write_h5ad('combined_for_pseudobulk.h5ad')
"

# Then aggregate
python data_prep/h5ad_to_pseudobulk.py combined_for_pseudobulk.h5ad -o my_pseudobulk/
```
