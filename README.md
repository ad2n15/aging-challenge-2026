# AI for Life Sciences Hackathon 2026 — Age Prediction Challenge

**13 April & 12 May 2026 · Avenue Campus, University of Southampton**

> Predict donor chronological age from single-cell RNA-seq data.  
> Contact: [IfLSAdmin@soton.ac.uk](mailto:IfLSAdmin@soton.ac.uk)

---

## The Challenge

You are given **peripheral blood mononuclear cell (PBMC) single-cell RNA-seq** data from the
[Onek1K cohort](https://www.science.org/doi/10.1126/science.abf3041) (981 donors, ages 19–79).
Your goal is to predict each donor's **chronological age** as accurately as possible using any
method you choose.

This is a **regression** task. Predictions are continuous (not integer).

---

## Timeline

| Date | Milestone |
|------|-----------|
| 8 April 2026 | Data released, competition opens |
| 13 April 2026 | Hackathon event — presentations & prizes |
| 12 May 2026 | Final event |

---

## Getting Started on Iridis X (University of Southampton HPC)

### Step 1 — Get the repository and container

The competition repo and the Bionemo container are pre-staged on the shared scratch space.
Copy them to your home directory:

```bash
cp -r /scratch/aazd1f17/shared_space/aging-challenge-2026 ~/aging-challenge-2026
```

This gives you:
```
~/aging-challenge-2026/
├── container/
│   └── bionemo-framework_nightly.sif   ← Apptainer container (all dependencies included)
├── notebooks/                           ← teaching notebooks
├── models/                              ← baseline scripts
├── data_prep/
└── ...
```

Alternatively, clone from GitHub and download the container separately:
```bash
git clone https://github.com/ad2n15/aging-challenge-2026.git ~/aging-challenge-2026
# Container is too large for git — copy it from shared space:
cp /scratch/aazd1f17/shared_space/aging-challenge-2026/container/bionemo-framework_nightly.sif \
   ~/aging-challenge-2026/container/
```

---

### Step 2 — Launch Jupyter via Open OnDemand

Open OnDemand lets you run Jupyter notebooks interactively in the Bionemo container
**without writing any SLURM scripts**.

1. Go to **[https://iridisondemand.soton.ac.uk/pun/sys/dashboard](https://iridisondemand.soton.ac.uk/pun/sys/dashboard)**
   and log in with your Iridis credentials.

2. Click the **"Jupyter with Apptainer Test"** icon.

3. Fill in the form:

   | Field | Value |
   |-------|-------|
   | **Working Directory** | `~/aging-challenge-2026` (or leave as `$HOME`) |
   | **User Interface** | `Jupyter Lab` |
   | **Submission Environment** | `container (advanced)` |
   | **Container File** | `/scratch/aazd1f17/shared_space/aging-challenge-2026/container/bionemo-framework_nightly.sif` — or, if you copied it locally: `~/aging-challenge-2026/container/bionemo-framework_nightly.sif` |
   | **Apptainer flags** | `--nv --bind $PWD:/workspace --pwd /workspace` |

4. Click **Launch** and wait for the session to start (~1–2 min).

5. In the Jupyter file browser, open any notebook under `notebooks/` to get started.

> **Note:** the `--bind /iridisfs/ddnb/Ahmed/data:...` flag mounts the competition data
> into the container. If the data has been copied to your own path, adjust accordingly.

---

### Step 3 — Run scripts from the terminal (batch jobs)

For longer runs (training, pseudobulk generation), submit via SLURM:

```bash
cd ~/aging-challenge-2026

# Train baseline model
sbatch run_binemo_AMD.sh models/train_age_model.py \
    --input PATH/TO/pseudobulk/combined_pseudobulk_donor_aggregated.h5ad

# Monitor job
tail -f slurm-JOBID.out
```

---

## Data

### Access

Download the competition data package from the link provided in the registration email.
**Iridis users:** data is already available — see Step 1 above.

```
aging_challenge_data/
├── onek1k/
│   ├── train.h5ad                               # 781 donors — age included
│   ├── val.h5ad                                 #  95 donors — age included
│   ├── test.h5ad                                # 105 donors — age EXCLUDED
│   ├── pseudobulk/
│   │   └── combined_pseudobulk_donor_aggregated.h5ad   # 981 donors × 141,390 features
│   └── donor_metadata.csv                       # donor_id, split, sex
└── README_data.md                               # data dictionary
```

### What is in the h5ad files

Each h5ad file contains:
- **`X`** — raw count matrix (cells × 35,477 genes, sparse)
- **`obs`** — cell-level metadata: `donor_id`, `celltype`, `age` *(train/val only)*, `_split`, QC columns

### Donor metadata

`donor_metadata.csv` contains: `donor_id`, `split`, `sex`, `sex_binary` (0=female, 1=male).
**Age is NOT included for test donors.**

### Gene names

Genes are stored as Ensembl IDs (e.g. `ENSG00000120163`).

---

## Baseline Approach

We provide a **Random Forest baseline** that:
1. Aggregates cell-level counts to **pseudobulk** (sum per donor × cell type)
2. Selects top-variance genes as features
3. Trains a Random Forest regressor

**Baseline validation performance (Onek1K):**

| MAE | RMSE | R² | Pearson | Spearman |
|-----|------|-----|---------|---------|
| ~9.2 years | ~12.1 | ~0.47 | ~0.69 | ~0.60 |

You are encouraged to improve on this using any technique.

---

## Repository Structure

```
aging_challenge_2026_public/
├── README.md                          ← you are here
├── data/
│   └── README.md                      ← data download instructions
├── data_prep/
│   └── h5ad_to_pseudobulk.py          ← aggregate cells → pseudobulk
├── models/
│   ├── train_age_model.py             ← baseline Random Forest (train + val eval)
│   ├── evaluate_val.py                ← evaluate your val predictions
│   └── README.md
├── notebooks/
│   ├── 01_anndata_and_pseudobulk.ipynb     ← understand the data format
│   ├── 02_baseline_model.ipynb             ← run and understand the baseline
│   ├── 03_evaluation_metrics.ipynb         ← understand MAE, R², Pearson, Spearman
│   └── 04_geneformer_embeddings.ipynb      ← advanced: foundation model features
└── submission/
    └── submission_template.csv            ← format your submission here
```

---

## Pseudobulk Pipeline (recommended starting point)

```bash
# Step 1 — build pseudobulk from the combined h5ad
python data_prep/h5ad_to_pseudobulk.py \
    path/to/combined.h5ad \
    -o path/to/pseudobulk_output/

# Step 2 — train the baseline model
python models/train_age_model.py \
    --input path/to/pseudobulk_output/combined_pseudobulk_donor_aggregated.h5ad

# Step 3 — evaluate on validation set
python models/evaluate_val.py \
    --predictions models/output/TIMESTAMP/val_predictions.csv
```

On the Iridis HPC cluster, prepend `sbatch run_binemo_AMD.sh` to each command.

---

## Submission Format

Submit a **CSV file** with exactly two columns:

```csv
sample_id,predicted_age
1,42.3
2,67.1
3,28.9
...
```

- `sample_id` = the anonymised integer `donor_id` from `test.h5ad`
- `predicted_age` = your predicted age (float, years)
- Must include **all 105 test donors** — missing donors score as MAE = 40

See `submission/submission_template.csv` for the list of test donor IDs.

Submit via email to [IfLSAdmin@soton.ac.uk](mailto:IfLSAdmin@soton.ac.uk) with subject:
`[AGE CHALLENGE SUBMISSION] Team Name`

Submissions are due before the hackathon event on **13 April 2026**.

---

## Evaluation Metrics

Submissions are ranked by **MAE (Mean Absolute Error)** — lower is better.

We also report RMSE, R², Pearson r, and Spearman ρ for reference.

See `notebooks/03_evaluation_metrics.ipynb` for a full explanation with examples.

---

## Ideas for Improvement

| Idea | Difficulty | Expected impact |
|------|-----------|-----------------|
| Tune RF hyperparameters (`--n-estimators`, `--max-depth`) | Low | Low–medium |
| Use more genes (`--n-genes 5000`) | Low | Low |
| Add sex as a feature (`--sex`) | Low | Low |
| Use all features (`--all-features`) | Low | Medium |
| Try gradient boosting (XGBoost / LightGBM) | Medium | Medium–high |
| Use Geneformer embeddings (see notebook 04) | Medium | Medium |
| Train on a different cell type subset | Medium | Medium |
| Use deep learning (MLP, attention) | High | High |
| Ensemble multiple models | High | High |

---

## Environment

The baseline is tested with:

```
python >= 3.10
scanpy >= 1.9
scikit-learn >= 1.3
pandas >= 1.5
numpy >= 1.24
scipy >= 1.10
joblib
```

On the Iridis cluster, use the Bionemo container — all dependencies are pre-installed.

---

## FAQ

**Q: Can I use external data?**  
A: Yes, as long as you disclose it in your presentation.

**Q: Can I use the validation set for training?**  
A: Yes, but your submission must include predictions for all 105 test donors.

**Q: Are predictions integers or floats?**  
A: Floats. Even though ground-truth ages are integers, your model will output continuous values.

**Q: How are ties broken?**  
A: By RMSE, then Pearson r.

**Q: Can I submit multiple times?**  
A: Yes — we score the last submission before the deadline.
