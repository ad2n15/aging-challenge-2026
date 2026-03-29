# Data Download

## Access

The competition data is available via the download link provided in your registration confirmation email.

If you did not receive the link, contact: [IfLSAdmin@soton.ac.uk](mailto:IfLSAdmin@soton.ac.uk)

## Files

After downloading and extracting, you will have:

```
aging_challenge_data/
├── onek1k/
│   ├── train.h5ad                                    # 781 donors, age included
│   ├── val.h5ad                                      #  95 donors, age included
│   ├── test.h5ad                                     # 105 donors, age EXCLUDED
│   ├── pseudobulk/
│   │   ├── combined_pseudobulk_donor_aggregated.h5ad # 981 donors × 141,390 features (recommended)
│   │   └── combined_pseudobulk_combined.h5ad         # (donor, celltype) × genes
│   └── donor_metadata.csv                            # donor_id, split, sex, sex_binary
└── README_data.md
```

## Placement

Place data relative to this repo root. Recommended layout:

```
aging_challenge_2026_public/
├── data_prep/
│   └── output/
│       ├── train.h5ad
│       ├── val.h5ad
│       ├── test.h5ad
│       ├── donor_metadata.csv
│       └── pseudobulk/
│           └── combined_pseudobulk_donor_aggregated.h5ad   ← start here
```

This matches the default paths in `models/train_age_model.py`. Or pass `--input` directly.

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
