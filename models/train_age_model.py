#!/usr/bin/env python3
"""
Train age prediction model from donor-aggregated pseudobulk.

Data modalities (mix and match):
  pseudobulk (default)  – gene expression aggregated per (donor, cell type)
  --geneformer          – Geneformer cell-level embeddings (see notebook 04)
  --geneformer-only     – use Geneformer as the ONLY base (skip pseudobulk h5ad)
  --sex                 – add donor sex as a binary feature (0=female, 1=male)

Use --compare-pca to benchmark multiple feature combinations side by side.

Input: combined_pseudobulk_donor_aggregated.h5ad (donors × features).
Uses Random Forest (scikit-learn); train/val/test splits are read from the '_split' obs column.

Usage:
  python models/train_age_model.py --input PATH/TO/donor_aggregated.h5ad
  python models/train_age_model.py --input PATH/TO/donor_aggregated.h5ad --n-genes 5000
  python models/train_age_model.py --input PATH/TO/donor_aggregated.h5ad --all-features
  python models/train_age_model.py --input PATH/TO/donor_aggregated.h5ad --sex
  python models/train_age_model.py --input PATH/TO/donor_aggregated.h5ad --compare-pca
  python models/train_age_model.py --input PATH/TO/donor_aggregated.h5ad --n-estimators 500
"""

import argparse
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import pearsonr, spearmanr

warnings.filterwarnings("ignore", category=UserWarning)

# Paths — override any of these with CLI arguments
PROJ_ROOT = Path(__file__).resolve().parents[1]
PSEUDOBULK_AGG_PATH = (
    PROJ_ROOT / "data" / "scRNA-seq_pseudobulk" / "train_pseudobulk_donor_aggregated_public.h5ad"
)
GENOTYPE_PCA_PATH = PROJ_ROOT / "data" / "genotypes" / "pca_train.tsv"
PRS_PATH = PROJ_ROOT / "data" / "genotypes" / "prs_donor_aligned.tsv"
GENEFORMER_DIR = PROJ_ROOT / "data" / "scRNA-seq_geneformer_pseudobulk"
DONOR_METADATA_PATH = PROJ_ROOT / "data" / "metadata" / "donor_metadata.csv"
OUTPUT_BASE         = PROJ_ROOT / "models" / "output"


# ---------------------------------------------------------------------------
# Feature selection
# ---------------------------------------------------------------------------

def select_genes(X, feature_names, n_genes=2000, n_celltypes=5):
    """
    Select top n_genes by variance. Features are gene__celltype; we group by gene
    (each gene has n_celltypes columns) and select genes with highest mean variance.
    """
    n_features = X.shape[1]
    n_genes_total = n_features // n_celltypes
    if n_genes >= n_genes_total:
        return X, feature_names

    var = np.var(X, axis=0)
    gene_vars = []
    for g in range(n_genes_total):
        start = g * n_celltypes
        end = start + n_celltypes
        gene_vars.append(np.mean(var[start:end]))
    gene_vars = np.array(gene_vars)
    top_genes = np.argsort(gene_vars)[::-1][:n_genes]

    keep_cols = []
    for g in top_genes:
        keep_cols.extend(range(g * n_celltypes, (g + 1) * n_celltypes))

    return X[:, keep_cols], [feature_names[i] for i in keep_cols]


# ---------------------------------------------------------------------------
# Geneformer loader
# ---------------------------------------------------------------------------

def load_geneformer_as_pseudobulk(geneformer_dir: Path):
    """
    Load the three Geneformer pseudobulk TSVs and return data in the same
    shape as the pseudobulk h5ad loader, so all downstream code is unchanged:
        X          – (N_donors, N_features) float64
        y          – (N_donors,) age; NaN for test donors (not needed)
        splits     – (N_donors,) str  'train' | 'val' | 'test'
        donor_ids  – (N_donors,) str
        feat_cols  – list[str]  5,760 embedding column names
    """
    dfs = {}
    for split in ["train", "val", "test"]:
        f = geneformer_dir / ("geneformer_pseudobulk_%s.tsv.gz" % split)
        if not f.exists():
            raise FileNotFoundError(
                "Geneformer file not found: %s\n"
                "Run data_prep/aggregate_geneformer_embeddings.py first." % f
            )
        dfs[split] = pd.read_csv(str(f), sep="\t")

    feat_cols = [c for c in dfs["train"].columns
                 if c not in ("donor_id", "split", "age")]

    frames = []
    for split_name, df in dfs.items():
        row = df.copy()
        if "age" not in row.columns:
            row["age"] = np.nan
        row["_split"] = split_name
        frames.append(row)
    combined = pd.concat(frames, ignore_index=True)

    X = combined[feat_cols].values.astype(np.float64)
    y = combined["age"].values.astype(np.float64)
    splits = combined["_split"].values
    donor_ids = combined["donor_id"].astype(str).values

    return X, y, splits, donor_ids, feat_cols


# ---------------------------------------------------------------------------
# Metrics helper
# ---------------------------------------------------------------------------

def compute_metrics(y_true, y_pred):
    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r2 = r2_score(y_true, y_pred)
    pearson_r, pearson_p = pearsonr(y_true, y_pred)
    spearman_r, spearman_p = spearmanr(y_true, y_pred)
    return {"MAE": mae, "RMSE": rmse, "R2": r2,
            "Pearson": pearson_r, "Pearson_p": pearson_p,
            "Spearman": spearman_r, "Spearman_p": spearman_p}


def print_metrics(label, m, prefix="  "):
    print("\n%s:" % label)
    print("%sMAE:      %.3f" % (prefix, m["MAE"]))
    print("%sRMSE:     %.3f" % (prefix, m["RMSE"]))
    print("%sR²:       %.3f" % (prefix, m["R2"]))
    print("%sPearson:  %.3f (p=%.2e)" % (prefix, m["Pearson"], m["Pearson_p"]))
    print("%sSpearman: %.3f (p=%.2e)" % (prefix, m["Spearman"], m["Spearman_p"]))


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------

def train_and_evaluate(X_train, y_train, X_val, y_val, X_test,
                       test_donors, feat_names,
                       n_estimators, max_depth, seed):
    rf = RandomForestRegressor(
        n_estimators=n_estimators,
        max_depth=max_depth,
        random_state=seed,
        n_jobs=-1,
    )
    rf.fit(X_train, y_train)

    val_metrics = compute_metrics(y_val, rf.predict(X_val))
    y_test_pred = rf.predict(X_test)
    pred_df = pd.DataFrame({"donor_id": test_donors, "predicted_age": y_test_pred})

    # Top 20 features / genes by importance
    imp = rf.feature_importances_
    feat_imp = {}
    for i, name in enumerate(feat_names):
        gene = name.rsplit("__", 1)[0] if "__" in name else name
        feat_imp[gene] = feat_imp.get(gene, 0) + imp[i]
    top20 = sorted(feat_imp.items(), key=lambda x: -x[1])[:20]
    top20_df = pd.DataFrame(top20, columns=["feature", "importance"])
    top20_df.insert(0, "rank", range(1, len(top20_df) + 1))

    return rf, val_metrics, pred_df, top20_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Train age prediction model")
    parser.add_argument("--input", type=str, default=str(PSEUDOBULK_AGG_PATH),
                        help="Donor-aggregated pseudobulk h5ad")
    parser.add_argument("--all-features", action="store_true",
                        help="Use all features (no variance-based gene selection)")
    parser.add_argument("--n-genes", type=int, default=2000,
                        help="Number of top-variance genes when not using --all-features")
    parser.add_argument("--n-estimators", type=int, default=100)
    parser.add_argument("--max-depth", type=int, default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory. Default: models/output/YYYY-MM-DD_HH-MM-SS/")
    parser.add_argument("--genotype-pca", type=str, default=None, nargs="?",
                        const=str(GENOTYPE_PCA_PATH),
                        help="Append genotype PCA features. Optionally pass a custom TSV path "
                             "(default: data/genotypes/pca_train.tsv; val/test PCs also under data/genotypes/).")
    parser.add_argument("--prs", type=str, default=None, nargs="?",
                        const=str(PRS_PATH),
                        help="Append PRS as a single feature column. Optionally pass a custom TSV path "
                             "(default: data/genotypes/prs_donor_aligned.tsv if present). "
                             "Requires organisers to ship PRS; otherwise omit.")
    parser.add_argument("--compare-pca", action="store_true",
                        help="Run all enabled modality combinations and print comparison.")
    parser.add_argument("--geneformer", type=str, default=None, nargs="?",
                        const=str(GENEFORMER_DIR),
                        help="Append Geneformer pseudobulk features (5,760 dims). Optionally "
                             "pass a custom directory path "
                             "(default: data/scRNA-seq_geneformer_pseudobulk/). "
                             "Files: geneformer_pseudobulk_{train,val,test}.tsv.gz")
    parser.add_argument("--geneformer-only", action="store_true",
                        help="Use ONLY Geneformer features as the base (skip pseudobulk h5ad). "
                             "Implies --geneformer with the default directory.")
    parser.add_argument("--sex", action="store_true",
                        help="Append donor sex as a binary feature (0=female, 1=male). "
                             "Requires donor_metadata.csv (see --donor-metadata).")
    parser.add_argument("--donor-metadata", type=str, default=None,
                        help="Path to donor_metadata.csv containing sex_binary column. "
                             "Default: data/metadata/donor_metadata.csv.")
    args = parser.parse_args()

    # --geneformer-only implies --geneformer
    if args.geneformer_only and args.geneformer is None:
        args.geneformer = str(GENEFORMER_DIR)

    if args.output_dir:
        run_dir = Path(args.output_dir)
    else:
        run_dir = OUTPUT_BASE / datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_dir.mkdir(parents=True, exist_ok=True)
    print("Output directory: %s" % run_dir)

    # ------------------------------------------------------------------
    # Load base features  (pseudobulk OR Geneformer-only)
    # ------------------------------------------------------------------
    if args.geneformer_only:
        gf_dir = Path(args.geneformer)
        print("Loading Geneformer features from %s ..." % gf_dir)
        X_agg, y, splits, donor_ids_sc, feature_names = load_geneformer_as_pseudobulk(gf_dir)
        print("  Train: %d  Val: %d  Test: %d  Features: %d" % (
            (splits == "train").sum(), (splits == "val").sum(),
            (splits == "test").sum(), X_agg.shape[1]))
        # No gene selection or log-transform for embeddings
        X_sub = X_agg
        feat_sub = feature_names
    else:
        print("Loading %s..." % args.input)
        adata = sc.read_h5ad(args.input)
        X_agg = adata.X
        if hasattr(X_agg, "toarray"):
            X_agg = X_agg.toarray()
        X_agg = np.asarray(X_agg, dtype=np.float64)
        feature_names = adata.var_names.tolist()
        obs_agg = adata.obs
        y = obs_agg["age"].values
        splits = obs_agg["_split"].values
        donor_ids_sc = obs_agg["donor_id"].astype(str).values
        print("  Shape: %d donors × %d features" % (X_agg.shape[0], X_agg.shape[1]))

        # Gene selection
        if args.all_features:
            X_sub = X_agg
            feat_sub = feature_names
            print("Using all %d features (--all-features)" % X_sub.shape[1])
        else:
            print("Selecting top %d genes by variance..." % args.n_genes)
            X_sub, feat_sub = select_genes(X_agg, feature_names, n_genes=args.n_genes)
            print("  Using %d features" % X_sub.shape[1])

        # log1p transform
        X_sub = np.log1p(X_sub)

    # Split masks (same for both loading paths)
    train_mask = splits == "train"
    val_mask   = splits == "val"
    test_mask  = splits == "test"
    test_donors = donor_ids_sc[test_mask]

    print("Train: %d, Val: %d, Test: %d" % (train_mask.sum(), val_mask.sum(), test_mask.sum()))

    # ------------------------------------------------------------------
    # Load optional extra features (genotype PCA, PRS, Geneformer-as-extra)
    # ------------------------------------------------------------------
    def load_donor_aligned(path: str | None, col_filter, label: str) -> tuple[np.ndarray | None, list[str]]:
        """Load a donor-aligned TSV, reindex to the current donor order, return (array, col_names)."""
        if path is None:
            return None, []
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError("%s file not found: %s" % (label, p))
        sep = "\t" if str(p).endswith((".tsv", ".tsv.gz")) else ","
        df = pd.read_csv(p, sep=sep)
        df["donor_id"] = df["donor_id"].astype(str)
        cols = [c for c in df.columns if col_filter(c)]
        aligned = df.set_index("donor_id")[cols].reindex(donor_ids_sc)
        if aligned.isna().any().any():
            raise ValueError("%d donors missing in %s." % (aligned.isna().any(axis=1).sum(), label))
        print("Loaded %s: %d donors × %d features from %s" % (label, len(aligned), len(cols), p))
        return aligned.values.astype(np.float64), cols

    pca_path = args.genotype_pca if args.genotype_pca else (
        str(GENOTYPE_PCA_PATH) if args.compare_pca else None
    )
    X_pca_aligned, pca_cols = load_donor_aligned(
        pca_path, lambda c: c.startswith("pca"), "genotype PCA"
    )
    X_prs_aligned, prs_cols = load_donor_aligned(
        args.prs if args.prs else None, lambda c: c == "prs", "PRS"
    )

    # Geneformer as an *extra* appended to pseudobulk (only when not geneformer-only)
    X_gf_aligned, gf_extra_cols = None, []
    if args.geneformer and not args.geneformer_only:
        gf_dir = Path(args.geneformer)
        print("Loading Geneformer features to append to pseudobulk ...")
        X_gf_all, _y_gf, _splits_gf, donor_ids_gf, gf_feat_cols = load_geneformer_as_pseudobulk(gf_dir)
        # Build a per-donor lookup from the geneformer data
        gf_df = pd.DataFrame(X_gf_all, columns=gf_feat_cols)
        gf_df.insert(0, "donor_id", donor_ids_gf)
        gf_df = gf_df.set_index("donor_id")
        aligned_gf = gf_df.reindex(donor_ids_sc)
        fully_missing = aligned_gf.isna().all(axis=1)
        if fully_missing.any():
            raise ValueError("%d donors completely absent from Geneformer data: %s"
                             % (fully_missing.sum(),
                                donor_ids_sc[fully_missing.values][:5].tolist()))
        partial_missing = aligned_gf.isna().any(axis=1).sum()
        if partial_missing:
            print("  Warning: %d donors have NaN in some cell-type embedding columns "
                  "(cell type absent in that split); will be filled with 0." % partial_missing)
            aligned_gf = aligned_gf.fillna(0.0)
        X_gf_aligned = aligned_gf.values.astype(np.float64)
        gf_extra_cols = gf_feat_cols
        print("  Geneformer appended: %d donors × %d features" % (X_gf_aligned.shape[0], X_gf_aligned.shape[1]))

    # Sex as a binary feature (0 = female, 1 = male)
    sex_meta_path = None
    if args.sex:
        sex_meta_path = args.donor_metadata if args.donor_metadata else str(DONOR_METADATA_PATH)
    X_sex_aligned, sex_cols = load_donor_aligned(
        sex_meta_path,
        lambda c: c == "sex_binary",
        "sex",
    )

    # ------------------------------------------------------------------
    # Run model(s)
    # ------------------------------------------------------------------
    import joblib

    def run_one(X_mat, feat_names, label, save_prefix):
        X_tr = X_mat[train_mask]
        X_va = X_mat[val_mask]
        X_te = X_mat[test_mask]
        print("\n" + "=" * 60)
        print("Configuration: %s" % label)
        print("  Features: %d" % X_mat.shape[1])
        rf, val_m, pred_df, top20 = train_and_evaluate(
            X_tr, y[train_mask], X_va, y[val_mask], X_te,
            test_donors, feat_names,
            args.n_estimators, args.max_depth, args.seed,
        )
        print_metrics("Validation metrics", val_m)
        pred_df.to_csv(run_dir / ("%s_test_predictions.csv" % save_prefix), index=False)
        joblib.dump(rf, run_dir / ("%s_rf_model.joblib" % save_prefix))
        top20.to_csv(run_dir / ("%s_top_features.csv" % save_prefix), index=False)
        # Save the full ordered feature list so predict.py can align new datasets
        pd.DataFrame({"feature": feat_names}).to_csv(
            run_dir / ("%s_feature_names.csv" % save_prefix), index=False)
        return val_m, pred_df, rf

    results = {}

    import shutil

    # Determine base label (pseudobulk or geneformer)
    base_label = "geneformer" if args.geneformer_only else "pseudobulk"

    if args.compare_pca:
        # ---------------------------------------------------------------
        # Comparison mode: run multiple feature combinations
        # ---------------------------------------------------------------
        configs = {}   # key → (X_matrix, feat_names_list, human_label)

        # Sex suffix appended to every config when --sex is set
        def _add_sex(X_mat, f_names):
            if X_sex_aligned is not None:
                return np.hstack([X_mat, X_sex_aligned]), f_names + sex_cols
            return X_mat, f_names

        # Base config
        X_base, f_base = _add_sex(X_sub, feat_sub)
        base_sfx = "_sex" if X_sex_aligned is not None else ""
        configs[base_label + base_sfx] = (X_base, f_base, base_label + (" + sex" if X_sex_aligned is not None else ""))

        # Base + genotype PCA (+ PRS if also requested)
        if X_pca_aligned is not None:
            extras_pca = [X_pca_aligned] + ([X_prs_aligned] if X_prs_aligned is not None else [])
            extra_cols_pca = pca_cols + (prs_cols if X_prs_aligned is not None else [])
            label_pca = base_label + " + genotype PCA" + (" + PRS" if X_prs_aligned is not None else "")
            key_pca = base_label + "_pca" + ("_prs" if X_prs_aligned is not None else "")
            X_c, f_c = _add_sex(np.hstack([X_sub] + extras_pca), feat_sub + extra_cols_pca)
            configs[key_pca + base_sfx] = (X_c, f_c, label_pca + (" + sex" if X_sex_aligned is not None else ""))

        # Base + PRS only (no PCA)
        if X_prs_aligned is not None and X_pca_aligned is None:
            X_c, f_c = _add_sex(np.hstack([X_sub, X_prs_aligned]), feat_sub + prs_cols)
            configs[base_label + "_prs" + base_sfx] = (X_c, f_c, base_label + " + PRS" + (" + sex" if X_sex_aligned is not None else ""))

        # Base + Geneformer (only makes sense when base is pseudobulk)
        if X_gf_aligned is not None:
            X_c, f_c = _add_sex(np.hstack([X_sub, X_gf_aligned]), feat_sub + gf_extra_cols)
            configs[base_label + "_geneformer" + base_sfx] = (X_c, f_c, base_label + " + Geneformer" + (" + sex" if X_sex_aligned is not None else ""))
            # Base + Geneformer + PCA
            if X_pca_aligned is not None:
                extras_all = [X_gf_aligned, X_pca_aligned] + (
                    [X_prs_aligned] if X_prs_aligned is not None else []
                )
                extra_cols_all = gf_extra_cols + pca_cols + (
                    prs_cols if X_prs_aligned is not None else []
                )
                label_all = base_label + " + Geneformer + PCA" + (" + PRS" if X_prs_aligned is not None else "")
                X_c, f_c = _add_sex(np.hstack([X_sub] + extras_all), feat_sub + extra_cols_all)
                configs[base_label + "_geneformer_pca" + base_sfx] = (X_c, f_c, label_all + (" + sex" if X_sex_aligned is not None else ""))

        results = {}
        for key, (X_mat, f_names, human_label) in configs.items():
            val_m, pred_df, rf = run_one(X_mat, f_names, human_label, key)
            results[key] = val_m

        # Summary table
        print("\n" + "=" * 60)
        print("COMPARISON (validation set)")
        print("%-40s %8s %8s %8s %8s %8s" % ("Config", "MAE", "RMSE", "R²", "Pearson", "Spearman"))
        for name, m in results.items():
            print("%-40s %8.3f %8.3f %8.3f %8.3f %8.3f" % (
                name, m["MAE"], m["RMSE"], m["R2"], m["Pearson"], m["Spearman"]))

        rows = [{"config": k, **v} for k, v in results.items()]
        pd.DataFrame(rows).to_csv(run_dir / "comparison.csv", index=False)
        print("\nComparison saved to %s" % (run_dir / "comparison.csv"))

        best = min(results, key=lambda k: results[k]["MAE"])
        shutil.copy(run_dir / ("%s_test_predictions.csv" % best), run_dir / "test_predictions.csv")
        print("Best config by MAE: %s → copied to test_predictions.csv" % best)

    else:
        # ---------------------------------------------------------------
        # Single run — assemble feature matrix from whatever is requested
        # ---------------------------------------------------------------
        extra_X = []
        extra_cols = []
        label_parts = [base_label]
        prefix = base_label

        if X_gf_aligned is not None:
            extra_X.append(X_gf_aligned)
            extra_cols += gf_extra_cols
            label_parts.append("Geneformer")
            prefix += "_geneformer"
        if X_pca_aligned is not None:
            extra_X.append(X_pca_aligned)
            extra_cols += pca_cols
            label_parts.append("genotype PCA")
            prefix += "_pca"
        if X_prs_aligned is not None:
            extra_X.append(X_prs_aligned)
            extra_cols += prs_cols
            label_parts.append("PRS")
            prefix += "_prs"
        if X_sex_aligned is not None:
            extra_X.append(X_sex_aligned)
            extra_cols += sex_cols
            label_parts.append("sex")
            prefix += "_sex"

        X_run = np.hstack([X_sub] + extra_X) if extra_X else X_sub
        feat_run = feat_sub + extra_cols
        label = " + ".join(label_parts)

        val_m, pred_df, rf = run_one(X_run, feat_run, label, prefix)
        shutil.copy(run_dir / ("%s_test_predictions.csv" % prefix), run_dir / "test_predictions.csv")

    print("\nDone. All outputs in %s" % run_dir)


if __name__ == "__main__":
    main()
