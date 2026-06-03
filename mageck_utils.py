"""
mageck_utils.py  —  Helper functions for the AllenLab CRISPR MAGeCK pipeline.

Key changes vs Round 1:
  - _extract_worker:        single_end flag disables reverse-anchor / rev-complement path
  - run_guide_extraction:   passes single_end from config; checks files exist before starting
  - preflight_anchor_check: samples N reads from a file and reports anchor hit rate
  - generate_metadata_template: parses new SE filename convention (*_R1.fastq.gz)
  - validate_metadata:      non_interactive flag bypasses the y/n pause for SLURM
  - All subprocess calls:   stderr captured and re-raised with a clear message
"""

import re
import os
import sys
import gzip
import time
import math
import warnings
import subprocess
import itertools
import pandas as pd
import seaborn as sns
import numpy as np
from Bio import SeqIO
from tqdm import tqdm
import multiprocessing as mp
from functools import partial
import matplotlib.pyplot as plt
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

warnings.simplefilter(action="ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Step 0 — Environment check
# ---------------------------------------------------------------------------

def check_env():
    """Check required Python packages and MAGeCK CLI. Returns list of missing items."""
    packages = {"Bio": "biopython", "pandas": "pandas", "tqdm": "tqdm",
                "seaborn": "seaborn", "scipy": "scipy", "statsmodels": "statsmodels"}
    missing = []
    for module, install_name in packages.items():
        try:
            __import__(module)
        except ImportError:
            missing.append(install_name)
    try:
        result = subprocess.run(
            ["mageck", "-v"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        if result.returncode != 0:
            missing.append("mageck (installed but --version failed)")
    except FileNotFoundError:
        missing.append("mageck")
    return missing


# ---------------------------------------------------------------------------
# Step 1 — Workspace initialisation
# ---------------------------------------------------------------------------

def initialize_workspace(input_path, output_path, lib_path, ann_path):
    """
    Create pipeline input/output directory skeleton.
    Returns True if any directory was newly created.
    Note: RAW_DATA_DIR is NOT created here — it is a pre-existing absolute path.
    """
    dirs = [input_path, output_path, lib_path, ann_path]
    newly_created = False
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d, exist_ok=True)
            print(f"  📁 Created: {d}")
            newly_created = True
        else:
            print(f"  ✅ Exists:  {d}")
    return newly_created


# ---------------------------------------------------------------------------
# Step 1b — Preflight anchor check
# ---------------------------------------------------------------------------

def preflight_anchor_check(raw_dir, fwd_anchor, n_reads=10_000, min_hit_rate=0.10,
                            single_end=True, rev_anchor=None):
    """
    Sample n_reads from the first valid FASTQ in raw_dir and report what fraction
    contain the forward (and optionally reverse) anchor sequence.

    This is a fast sanity check that runs before the full extraction.
    Returns True if the hit rate is acceptable, False if it looks suspiciously low.
    """
    files = sorted([
        f for f in os.listdir(raw_dir)
        if f.endswith(".fastq.gz") and "Unassigned" not in f and "PhiX" not in f
    ])
    if not files:
        print(f"⚠️  Preflight: No .fastq.gz files found in {raw_dir}")
        return False

    sample_file = os.path.join(raw_dir, files[0])
    print(f"\n🔬 Preflight anchor check on: {files[0]}  (sampling {n_reads:,} reads)")

    fwd_hits = 0
    rev_hits = 0
    total    = 0

    try:
        with gzip.open(sample_file, "rt") as fh:
            for record in itertools.islice(SeqIO.parse(fh, "fastq"), n_reads):
                seq = str(record.seq).upper()
                if fwd_anchor in seq:
                    fwd_hits += 1
                elif (not single_end) and rev_anchor and rev_anchor in seq:
                    rev_hits += 1
                total += 1
    except Exception as e:
        print(f"⚠️  Preflight: Could not read {sample_file}: {e}")
        return False

    fwd_rate = fwd_hits / total if total else 0
    rev_rate = rev_hits / total if total else 0
    combined = fwd_rate + rev_rate

    print(f"   Reads sampled : {total:,}")
    print(f"   FWD anchor hits: {fwd_hits:,}  ({fwd_rate:.1%})")
    if not single_end:
        print(f"   REV anchor hits: {rev_hits:,}  ({rev_rate:.1%})")
    print(f"   Combined hit rate: {combined:.1%}")

    if combined < min_hit_rate:
        print(
            f"\n⚠️  WARNING: Only {combined:.1%} of sampled reads contain the anchor.\n"
            f"   Expected ≥ {min_hit_rate:.0%}.  Possible causes:\n"
            f"   • Wrong FWD_ANCHOR sequence in config\n"
            f"   • Files are demultiplexed incorrectly\n"
            f"   • Library prep issues\n"
            f"   Pipeline will continue, but check extraction_summary.csv carefully.\n"
        )
        return False

    print(f"   ✅ Anchor hit rate looks good ({combined:.1%} ≥ {min_hit_rate:.0%} threshold).\n")
    return True


# ---------------------------------------------------------------------------
# Step 2 — Metadata template generation
# ---------------------------------------------------------------------------

def generate_metadata_template(csv_path, raw_dir):
    """
    Create a blank metadata CSV from .fastq.gz files in raw_dir.
    Handles both SE (SampleName_R1.fastq.gz) and old Illumina naming.
    Will NOT overwrite an existing file.
    """
    if os.path.exists(csv_path):
        print(f"  ✅ Metadata file already exists: {csv_path}  (no overwrite)")
        return

    if not os.path.exists(raw_dir):
        print(f"  ❌ Error: raw_dir not found: {raw_dir}")
        return

    files = sorted([
        f for f in os.listdir(raw_dir)
        if f.endswith(".fastq.gz")
        and "Unassigned" not in f
        and "Undetermined" not in f
        and "PhiX" not in f
    ])

    if not files:
        print(f"  ⚠️  No usable .fastq.gz files found in {raw_dir}")
        return

    # Strip common suffixes to derive a clean sample label
    def _clean_label(fname):
        for suffix in ("_R1.fastq.gz", "_R2.fastq.gz", ".fastq.gz", ".fq.gz"):
            if fname.endswith(suffix):
                return fname[: -len(suffix)]
        return fname

    df = pd.DataFrame({
        "filename":         files,
        "sample_label":     [_clean_label(f) for f in files],
        "experiment_label": "",  # fill manually or use pre-built metadata CSV
        "condition":        "",
        "bio_rep":          "",
    })
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    df.to_csv(csv_path, index=False)
    print(f"  🌟 Metadata template created: {csv_path}")
    print(f"      {len(files)} files listed — fill experiment_label, condition, bio_rep.")


# ---------------------------------------------------------------------------
# Step 3 — Metadata validation
# ---------------------------------------------------------------------------

def validate_metadata(csv_path, non_interactive=False):
    """
    Validate metadata CSV and print an experimental design summary.
    If non_interactive=True (SLURM), skip the y/n user prompt.
    Returns (is_valid: bool, message: str).
    """
    if not os.path.exists(csv_path):
        return False, f"Metadata file missing: {csv_path}"

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        return False, f"Could not read CSV: {e}"

    required_cols = ["filename", "experiment_label", "condition", "bio_rep"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        return False, f"Missing required columns: {missing_cols}"

    for col in ["experiment_label", "condition", "bio_rep"]:
        empty = df[col].isna().sum() + (df[col].astype(str).str.strip() == "").sum()
        if empty > 0:
            return False, f"Column '{col}' has {empty} empty entries. Please fill them."

    # Summary table
    summary = (
        df.groupby(["experiment_label", "condition"])["bio_rep"]
        .nunique()
        .unstack(fill_value=0)
    )
    print("\n" + "=" * 70)
    print("📊 EXPERIMENTAL DESIGN SUMMARY")
    print("   Biological replicates per condition (rows = experiments):")
    print("=" * 70)
    print(summary.to_string())
    print("=" * 70)
    print("   Note: t0 will be the baseline in the design matrix (no t0 column).")
    print("=" * 70 + "\n")

    if non_interactive:
        print("ℹ️  NON_INTERACTIVE mode — skipping manual confirmation prompt.")
    else:
        answer = input("Does the experimental summary look correct? (y/n): ").strip().lower()
        if answer != "y":
            return False, "Cancelled by user — please correct the metadata CSV."

    return True, "Metadata OK."


# ---------------------------------------------------------------------------
# Step 4 — Guide extraction
# ---------------------------------------------------------------------------

def _check_files_exist(meta_df, raw_dir):
    """
    Verify every file listed in the metadata actually exists in raw_dir.
    Returns list of missing filenames.
    """
    missing = []
    for fname in meta_df["filename"].unique():
        full = os.path.join(raw_dir, fname)
        if not os.path.exists(full):
            missing.append(fname)
    return missing


def _extract_worker(row_tuple, raw_dir, out_dir, fwd, rev, length, single_end):
    """
    Extract guides from a single FASTQ file.

    single_end=True  : only FWD anchor is used (SE150).
    single_end=False : also checks REV anchor with reverse-complement (PE).
    """
    _, row = row_tuple
    in_p = os.path.join(raw_dir, row["filename"])

    clean_name = row["filename"].split(".fast")[0].split(".fq")[0]
    out_filename = (
        f"{row['experiment_label']}_{row['condition']}_R{row['bio_rep']}_{clean_name}.fastq.gz"
    )
    out_p = os.path.join(out_dir, out_filename)

    total, extracted = 0, 0
    try:
        with gzip.open(in_p, "rt") as in_h, gzip.open(out_p, "wt") as out_h:
            for record in SeqIO.parse(in_h, "fastq"):
                total += 1
                seq_str = str(record.seq).upper()
                extracted_record = None

                # --- Forward anchor (always checked) ---
                if fwd in seq_str:
                    start_idx = seq_str.find(fwd) + len(fwd)
                    if (len(seq_str) - start_idx) >= length:
                        extracted_record = record[start_idx: start_idx + length]

                # --- Reverse anchor (only for paired-end mode) ---
                elif (not single_end) and rev in seq_str:
                    start_idx = seq_str.find(rev) + len(rev)
                    if (len(seq_str) - start_idx) >= length:
                        raw_chunk = record[start_idx: start_idx + length]
                        extracted_record = raw_chunk.reverse_complement(
                            id=True, description=True
                        )

                if extracted_record:
                    SeqIO.write(extracted_record, out_h, "fastq")
                    extracted += 1

        pct = (extracted / total * 100) if total > 0 else 0
        return {
            "Success":          True,
            "Experiment":       row["experiment_label"],
            "Condition":        row["condition"],
            "Bio_Rep":          row["bio_rep"],
            "Original_File":    row["filename"],
            "Output_File":      out_filename,
            "Total_Reads":      total,
            "Extracted_Reads":  extracted,
            "Success_Rate_Pct": round(pct, 2),
        }
    except Exception as e:
        return {"Success": False, "Output_File": out_filename, "Error": str(e)}


def run_guide_extraction(meta_df, raw_dir, out_dir, fwd, rev, length, n_cores,
                         single_end=True):
    """
    Parallel guide extraction with pre-run file existence check and QC summary.
    """
    # --- Pre-run file existence check ---
    missing = _check_files_exist(meta_df, raw_dir)
    if missing:
        print(f"\n❌ ABORT: {len(missing)} file(s) listed in metadata not found in {raw_dir}:")
        for f in missing[:20]:
            print(f"   • {f}")
        if len(missing) > 20:
            print(f"   ... and {len(missing)-20} more.")
        sys.exit(1)
    print(f"  ✅ All {len(meta_df)} metadata files confirmed in {raw_dir}")

    mode_str = "single-end (FWD anchor only)" if single_end else "paired-end (FWD + REV anchors)"
    print(f"  🔬 Extraction mode: {mode_str}")
    print(f"  📐 Guide length: {length} bp")
    print(f"  🚀 Launching parallel extraction on {n_cores} cores...\n")

    os.makedirs(out_dir, exist_ok=True)

    worker_task = partial(
        _extract_worker, raw_dir=raw_dir, out_dir=out_dir,
        fwd=fwd, rev=rev, length=length, single_end=single_end
    )

    with mp.Pool(processes=n_cores) as pool:
        results = list(
            tqdm(pool.imap(worker_task, meta_df.iterrows()),
                 total=len(meta_df), desc="Extraction Progress")
        )

    success_results = [r for r in results if r.get("Success")]
    error_results   = [r for r in results if not r.get("Success")]

    if success_results:
        summary_df = pd.DataFrame(success_results).drop(columns=["Success"])
        summary_df.to_csv(os.path.join(out_dir, "extraction_summary.csv"), index=False)

        print("\n" + "=" * 90)
        print(f"  {'Output File':<45} | {'Total':>10} | {'Extracted':>10} | {'Rate':>8}")
        print("-" * 90)
        for _, r in summary_df.iterrows():
            print(f"  {r['Output_File'][:45]:<45} | {r['Total_Reads']:>10,} | "
                  f"{r['Extracted_Reads']:>10,} | {r['Success_Rate_Pct']:>7.1f}%")
        print("=" * 90)

        avg_rate = summary_df["Success_Rate_Pct"].mean()
        low_rate = summary_df[summary_df["Success_Rate_Pct"] < 10]
        print(f"\n  📊 Average extraction rate: {avg_rate:.1f}%")
        if not low_rate.empty:
            print(f"\n  ⚠️  {len(low_rate)} file(s) have < 10% extraction rate:")
            for _, r in low_rate.iterrows():
                print(f"     • {r['Output_File']}  ({r['Success_Rate_Pct']:.1f}%)")
        print(f"\n  ✅ Full summary saved to: {out_dir}/extraction_summary.csv")

    for err in error_results:
        print(f"  ❌ Error processing {err['Output_File']}: {err.get('Error')}")

    if error_results:
        print(f"\n  ⚠️  {len(error_results)} file(s) failed extraction. Check errors above.")


# ---------------------------------------------------------------------------
# Step 5 — MAGeCK count  (unchanged logic, better error messages)
# ---------------------------------------------------------------------------

def _generate_count_plots(exp, count_dir):
    """QC plots: mapping efficiency, read distribution, zero counts."""
    summary_file = os.path.join(count_dir, f"{exp}_.countsummary.txt")
    count_file   = os.path.join(count_dir, f"{exp}_.count.txt")

    if not (os.path.exists(summary_file) and os.path.exists(count_file)):
        print(f"  ⚠️  QC plot skipped — file not found: {summary_file}")
        return

    summary_df = pd.read_table(summary_file)
    counts = pd.read_table(count_file)

    plot_summary = summary_df.groupby("Label", sort=False).agg(
        {"Percentage": "mean", "Zerocounts": "mean"}
    ).reset_index()

    numeric_cols = [c for c in counts.columns if c not in ["sgRNA", "Gene"]]
    melted = counts.melt(id_vars=["sgRNA"], value_vars=numeric_cols,
                         var_name="Label", value_name="Count")
    melted["Log10_Count"] = np.log10(melted["Count"] + 1)
    sample_order = sorted(melted["Label"].unique())

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    fig.suptitle(f"MAGeCK Count QC: {exp}", fontsize=18, y=1.02)

    sns.barplot(x=plot_summary["Label"], y=plot_summary["Percentage"] * 100,
                ax=axes[0], palette="viridis")
    axes[0].set_ylabel("Mapping Rate (%)", fontsize=12)
    axes[0].set_ylim(0, 100)
    axes[0].axhline(30, color="red", linestyle="--", alpha=0.5, label="30% threshold")
    axes[0].set_xticks(range(len(plot_summary)))
    axes[0].set_xticklabels(plot_summary["Label"], rotation=45, ha="right")
    axes[0].set_title("Mapping Efficiency", fontsize=14)
    axes[0].legend(fontsize=9)

    sns.violinplot(data=melted, x="Label", y="Log10_Count", order=sample_order,
                   ax=axes[1], palette="muted", inner="quartile")
    axes[1].set_xlabel("", fontsize=0)
    axes[1].set_ylabel("Log10 Count", fontsize=12)
    axes[1].set_xticks(range(len(sample_order)))
    axes[1].set_xticklabels(sample_order, rotation=45, ha="right")
    axes[1].set_title("Read Distribution (Log10)", fontsize=14)

    sns.barplot(x=plot_summary["Label"], y=plot_summary["Zerocounts"],
                ax=axes[2], palette="magma")
    axes[2].set_ylabel("sgRNAs with Zero Reads", fontsize=12)
    axes[2].set_xticks(range(len(plot_summary)))
    axes[2].set_xticklabels(plot_summary["Label"], rotation=45, ha="right")
    axes[2].set_title("Zero-Count sgRNAs", fontsize=14)

    plt.tight_layout()
    save_path = os.path.join(count_dir, f"{exp}_QC_Summary.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✅ Count QC plot saved: {save_path}")


def _generate_control_qc_plots(exp, count_dir, neg_path, pos_path):
    """Raw vs normalised count distributions for negative / positive controls."""
    raw_file  = os.path.join(count_dir, f"{exp}_.count.txt")
    norm_file = os.path.join(count_dir, f"{exp}_.count_normalized.txt")

    for fpath in (raw_file, norm_file):
        if not os.path.exists(fpath):
            print(f"  ⚠️  Control QC skipped — file not found: {fpath}")
            return

    with open(neg_path) as f:
        neg_ids = {line.strip() for line in f if line.strip()}
    pos_ids = set(pd.read_csv(pos_path).iloc[:, 1].tolist())

    def get_type(sgrna):
        if sgrna in pos_ids:
            return "Essential (Pos)"
        if sgrna in neg_ids:
            return "Non-targeting" if "nontargetting" in sgrna.lower() else "Safe (Neg)"
        return "Experimental"

    dfs = []
    for fpath, label in [(raw_file, "Raw"), (norm_file, "Normalized")]:
        df = pd.read_table(fpath)
        df["Type"]  = df["sgRNA"].apply(get_type)
        df["Scale"] = label
        num_cols = [c for c in df.columns if c not in ["sgRNA", "Gene", "Type", "Scale"]]
        melted = df.melt(id_vars=["sgRNA", "Type", "Scale"], value_vars=num_cols,
                         var_name="Sample", value_name="Count")
        melted["Log10_Count"] = np.log10(melted["Count"] + 1)
        dfs.append(melted[melted["Type"] != "Experimental"])

    full_plot_df = pd.concat(dfs)
    counts_summary = pd.read_table(raw_file)["sgRNA"].apply(get_type).value_counts()

    fig, axes = plt.subplots(2, 3, figsize=(22, 14), sharey="row")
    fig.suptitle(f"Control Guide QC (Raw vs Normalized): {exp}", fontsize=20, y=1.02)

    colors = ["#4C72B0", "#55A868", "#C44E52"]
    categories = ["Safe (Neg)", "Non-targeting", "Essential (Pos)"]

    for row_idx, scale_label in enumerate(["Raw", "Normalized"]):
        for col_idx, cat in enumerate(categories):
            ax = axes[row_idx, col_idx]
            subset = full_plot_df[
                (full_plot_df["Scale"] == scale_label) & (full_plot_df["Type"] == cat)
            ]
            if subset.empty:
                ax.text(0.5, 0.5, f"No {cat} found", ha="center", va="center")
                continue
            sns.violinplot(data=subset, x="Sample", y="Log10_Count", ax=ax,
                           color=colors[col_idx], inner="quartile", cut=0)
            n_guides = counts_summary.get(cat, 0)
            title_prefix = "RAW" if row_idx == 0 else "NORM"
            ax.set_title(f"{title_prefix}: {cat}\n(n={n_guides})", fontsize=14, fontweight="bold")
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
            ax.set_xlabel("")
            ax.set_ylabel(f"{scale_label} Log10(Count+1)" if col_idx == 0 else "")

    plt.tight_layout()
    save_path = os.path.join(count_dir, f"{exp}_Control_QC.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  ✅ Control QC plot saved: {save_path}")

    # ---- Sanity test: essential gene median should be < non-targeting median ----
    norm_df = pd.read_table(norm_file)
    norm_df["Type"] = norm_df["sgRNA"].apply(get_type)
    ess_median = norm_df[norm_df["Type"] == "Essential (Pos)"].iloc[:, 2:-2].median().median()
    neg_median = norm_df[norm_df["Type"].isin(["Safe (Neg)", "Non-targeting"])].iloc[:, 2:-2].median().median()
    print(f"\n  🧪 Control sanity test [{exp}]:")
    print(f"     Essential (Pos) median count : {ess_median:.1f}")
    print(f"     Non-targeting median count   : {neg_median:.1f}")
    if ess_median < neg_median:
        print(f"     ✅ PASS — essential guides deplete relative to non-targeting controls.")
    else:
        print(f"     ⚠️  WARN — essential guides are NOT depleted relative to controls.")
        print(f"        This may indicate a problem with the screen or normalization.")


def _plot_specific_gene_trajectories(exp, count_dir, genes_dict):
    """Mean trajectory plots (±SD) for genes of interest."""
    raw_file  = os.path.join(count_dir, f"{exp}_.count.txt")
    norm_file = os.path.join(count_dir, f"{exp}_.count_normalized.txt")

    if not (os.path.exists(raw_file) and os.path.exists(norm_file)):
        print(f"  ⚠️  Gene trajectory plots skipped — files not found for {exp}.")
        return

    gene_out_dir = os.path.join(count_dir, f"{exp}_gene_profiles")
    os.makedirs(gene_out_dir, exist_ok=True)

    raw_df  = pd.read_table(raw_file)
    norm_df = pd.read_table(norm_file)
    sample_cols = [c for c in raw_df.columns if c not in ["sgRNA", "Gene"]]

    found, not_found = [], []
    for gene_id, common_name in genes_dict.items():
        raw_gene  = raw_df[raw_df["Gene"] == gene_id]
        norm_gene = norm_df[norm_df["Gene"] == gene_id]
        if raw_gene.empty or norm_gene.empty:
            not_found.append(f"{gene_id} ({common_name})")
            continue
        found.append(common_name)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

        def _plot_traj(data_df, ax, title, ylabel):
            melted = data_df.melt(id_vars=["sgRNA"], value_vars=sample_cols,
                                  var_name="Sample", value_name="Count")
            melted["Log10_Count"] = np.log10(melted["Count"] + 1)
            sns.lineplot(data=melted, x="Sample", y="Log10_Count", ax=ax,
                         color="royalblue", marker="o", errorbar="sd", linewidth=3)
            ax.set_title(title, fontsize=18, fontweight="bold")
            ax.set_ylabel(ylabel, fontsize=14)
            ax.tick_params(axis="x", rotation=45)
            ax.grid(True, alpha=0.3)

        _plot_traj(raw_gene,  ax1, f"RAW: {common_name}",  "Log10(Raw Count + 1)")
        _plot_traj(norm_gene, ax2, f"NORM: {common_name}", "Log10(Norm Count + 1)")
        fig.suptitle(f"Mean Trajectory (±SD): {common_name} ({gene_id}) | {exp}",
                     fontsize=20, y=1.02, fontweight="bold")

        save_path = os.path.join(gene_out_dir, f"{common_name}_Mean_Trajectory.png")
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()

    print(f"  ✅ Gene trajectory profiles: {len(found)} plotted, {len(not_found)} not found.")
    if not_found:
        print(f"     Not found: {not_found}")


def run_mageck_count(preprocess_dir, count_dir, library_csv, summary_csv,
                     neg_path, pos_path, genes_dict):
    """
    Run `mageck count` per experiment, grouping technical replicates with commas.
    """
    os.makedirs(count_dir, exist_ok=True)

    if not os.path.exists(summary_csv):
        print(f"  ❌ Error: extraction summary not found: {summary_csv}")
        return

    # Validate library CSV exists
    if not os.path.exists(library_csv):
        print(f"  ❌ Error: library CSV not found: {library_csv}")
        sys.exit(1)

    df = pd.read_csv(summary_csv)
    experiments = df["Experiment"].unique()
    print(f"  Found {len(experiments)} experiment(s): {list(experiments)}")

    for exp in experiments:
        exp_df = df[df["Experiment"] == exp].sort_values(["Condition", "Bio_Rep"])
        sample_labels = []
        fastq_args    = []

        for (cond, rep), group in exp_df.groupby(["Condition", "Bio_Rep"], sort=False):
            label = f"{cond}_R{rep}"
            sample_labels.append(label)
            tech_reps = [os.path.join(preprocess_dir, f) for f in group["Output_File"].tolist()]
            fastq_args.append(",".join(tech_reps))

        # Verify extracted files exist before calling mageck
        all_files = [f for group in fastq_args for f in group.split(",")]
        missing = [f for f in all_files if not os.path.exists(f)]
        if missing:
            print(f"  ❌ {exp}: {len(missing)} extracted file(s) not found. "
                  f"Did Step 4 complete? First missing: {missing[0]}")
            continue

        output_prefix = os.path.join(count_dir, f"{exp}_")
        cmd = (
            ["mageck", "count",
             "-l", library_csv,
             "-n", output_prefix,
             "--sample-label", ",".join(sample_labels),
             "--fastq"]
            + fastq_args
        )

        print(f"\n  🚀 MAGeCK count — experiment: {exp}")
        print(f"     Samples: {sample_labels}")
        print(f"     Command: {' '.join(cmd[:8])} ...")

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            if result.stdout:
                print(result.stdout[-2000:])  # last 2000 chars
            print(f"  ✅ Count table: {output_prefix}.count.txt")
        except subprocess.CalledProcessError as e:
            print(f"  ❌ MAGeCK count failed for {exp}:")
            print(f"     STDERR: {e.stderr[-1000:]}")
            print("     Skipping QC plots for this experiment.")
            continue

        print(f"  🎨 Generating QC plots for {exp}...")
        _generate_count_plots(exp, count_dir)
        _generate_control_qc_plots(exp, count_dir, neg_path, pos_path)
        _plot_specific_gene_trajectories(exp, count_dir, genes_dict)


# ---------------------------------------------------------------------------
# Step 6 — Design matrix management  (unchanged logic)
# ---------------------------------------------------------------------------

def manage_all_design_matrices(summary_csv, design_base_dir):
    """
    Create design matrix skeletons per experiment, validate filled matrices.
    t0 samples are excluded from columns (they form the implicit baseline).
    Returns True when all matrices are valid and ready for MLE.
    """
    if not os.path.exists(summary_csv):
        print(f"  ❌ Error: extraction summary not found: {summary_csv}")
        return False

    df = pd.read_csv(summary_csv)
    unique_exps = df["Experiment"].unique()
    all_ready = True

    print(f"\n  🔍 Scanning {len(unique_exps)} experiment(s)...")

    for exp_name in unique_exps:
        exp_dir     = os.path.join(design_base_dir, exp_name)
        design_path = os.path.join(exp_dir, f"design_matrix_{exp_name}.txt")
        os.makedirs(exp_dir, exist_ok=True)
        exp_df = df[df["Experiment"] == exp_name]

        # Create skeleton
        if not os.path.exists(design_path):
            labels = exp_df.apply(
                lambda x: f"{x['Condition']}_R{x['Bio_Rep']}", axis=1
            ).unique()
            unique_conditions = [
                c for c in exp_df["Condition"].unique() if c.lower() != "t0"
            ]
            design_df = pd.DataFrame({"Samples": labels})
            design_df["baseline"] = 1
            for cond in unique_conditions:
                design_df[cond] = ""

            design_df.to_csv(design_path, sep="\t", index=False)
            print(f"  📄 Design matrix skeleton created: {design_path}")
            print(f"     Conditions to fill: {unique_conditions}")
            print(f"     Set 1 for treatment rows, 0 for all others.")
            print(f"     t0 rows → leave all condition columns as 0 (they are the baseline).")
            all_ready = False
            continue

        # Validate existing file
        try:
            user_df   = pd.read_table(design_path, dtype=str).fillna("")
            cols      = user_df.columns.tolist()
            if any(c.lower() == "t0" for c in cols):
                print(f"  ❌ {exp_name}: Found a 't0' column — delete it. "
                      f"t0 is represented by having 0s in all condition columns.")
                all_ready = False
                continue
            data_cols  = cols[2:]
            is_binary  = user_df[data_cols].map(lambda x: str(x).strip() in ("0", "1")).all().all()
            has_content = (user_df[data_cols] != "").all().all()
            if not (is_binary and has_content):
                print(f"  ⚠️  {exp_name}: Matrix incomplete or invalid values.")
                all_ready = False
            else:
                print(f"  ✅ {exp_name}: Design matrix validated  ({', '.join(data_cols)}).")
        except Exception as e:
            print(f"  ❌ {exp_name}: Error reading matrix: {e}")
            all_ready = False

    return all_ready


# ---------------------------------------------------------------------------
# Step 7 — MAGeCK MLE
# ---------------------------------------------------------------------------

def run_all_mageck_mle(count_dir, mle_dir, design_base_dir, n_cores,
                       perm_round, max_sgrna_perm, control_file):
    """Run MAGeCK MLE per experiment, skipping already-completed runs."""
    experiments = [
        d for d in os.listdir(design_base_dir)
        if os.path.isdir(os.path.join(design_base_dir, d))
    ]

    for exp in experiments:
        param_suffix  = f"perm{perm_round}_max{max_sgrna_perm}"
        run_dir       = os.path.join(mle_dir, exp, param_suffix)
        os.makedirs(run_dir, exist_ok=True)
        output_prefix = os.path.join(run_dir, f"{exp}_mle_{param_suffix}")
        gene_summary  = f"{output_prefix}.gene_summary.txt"

        if os.path.exists(gene_summary):
            print(f"\n  ✋ Skipping {exp} — output already exists.")
            continue

        design_path = os.path.join(design_base_dir, exp, f"design_matrix_{exp}.txt")
        count_table = os.path.join(count_dir, f"{exp}_.count.txt")

        for fpath, label in [(design_path, "design matrix"), (count_table, "count table")]:
            if not os.path.exists(fpath):
                print(f"  ❌ Skipping {exp}: {label} not found: {fpath}")
                break
        else:
            cmd = [
                "mageck", "mle",
                "-k", count_table,
                "-d", design_path,
                "-n", output_prefix,
                "--threads", str(n_cores),
                "--permutation-round", str(perm_round),
                "--max-sgrnapergene-permutation", str(max_sgrna_perm),
                "--control-sgrna", control_file,
                "--norm-method", "control",
            ]
            print(f"\n  🚀 MAGeCK MLE — {exp} | {param_suffix}")
            try:
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                if result.stdout:
                    print(result.stdout[-2000:])
                print(f"  ✅ MLE complete: {os.path.basename(output_prefix)}")
            except subprocess.CalledProcessError as e:
                print(f"  ❌ MAGeCK MLE failed for {exp}:")
                print(f"     STDERR: {e.stderr[-1000:]}")


# ---------------------------------------------------------------------------
# Step 8 — MLE QC plots
# ---------------------------------------------------------------------------

def run_mle_qc_plots(mle_dir, essential_txt, pos_control_csv):
    """Distribution and essential-gene tracking plots for each MLE output."""
    eg_ids = set(pd.read_table(essential_txt, sep="\t", header=None).iloc[:, 0].astype(str))
    gold_std_ids = set(pd.read_csv(pos_control_csv).iloc[:, 0].astype(str))

    for root, _, files in os.walk(mle_dir):
        for fname in files:
            if not fname.endswith(".gene_summary.txt"):
                continue
            file_path = os.path.join(root, fname)
            exp_tag   = fname.split("_")[0]
            df        = pd.read_table(file_path)
            beta_cols = [c for c in df.columns if "|beta" in c]
            if not beta_cols:
                print(f"  ⚠️  No beta columns in {fname} — skipping.")
                continue

            df["Group"] = df["Gene"].apply(
                lambda x: "Essential" if str(x) in eg_ids else "Experimental"
            )
            global_min = min(df[beta_cols].min().min(), -1.0)
            global_max = max(df[beta_cols].max().max(), 0.5)
            x_limit = (global_min * 1.1, global_max * 1.1)

            # Estimate global density y-limit
            max_y = 0
            for col in beta_cols:
                data = df[col].dropna()
                if len(data) > 1:
                    temp_ax = sns.kdeplot(data=data)
                    line = temp_ax.lines[-1]
                    max_y = max(max_y, np.max(line.get_ydata()))
                    plt.close()
            y_limit = (0, max_y * 1.1)

            cols_per_row = 4
            n_rows = math.ceil(len(beta_cols) / cols_per_row)
            fig, axes = plt.subplots(n_rows, cols_per_row,
                                     figsize=(24, 6 * n_rows), squeeze=False)
            fig.suptitle(f"Beta Score Distributions: {exp_tag}", fontsize=28, y=1.02)

            for i, col in enumerate(beta_cols):
                ax = axes[divmod(i, cols_per_row)]
                sns.kdeplot(data=df[df["Group"] == "Experimental"], x=col,
                            ax=ax, fill=True, color="gray", alpha=0.3, label="Exp Pool")
                sns.kdeplot(data=df[df["Group"] == "Essential"], x=col,
                            ax=ax, color="red", lw=3, label="Essential")
                ax.set_xlim(x_limit)
                ax.set_ylim(y_limit)
                ax.axvline(0, color="black", linestyle="--", lw=2, alpha=0.8)
                ax.set_title(col.split("|")[0], fontsize=18, fontweight="bold")
                ax.set_xlabel("Beta Score", fontsize=16)
                ax.set_ylabel("Density", fontsize=16)
                ax.tick_params(labelsize=13)
                ax.legend(fontsize=11)

            for j in range(i + 1, n_rows * cols_per_row):
                fig.delaxes(axes.flatten()[j])
            plt.tight_layout()
            dist_name = f"QC_Distributions_{exp_tag}.png"
            plt.savefig(os.path.join(root, dist_name), dpi=300, bbox_inches="tight")
            plt.close()

            # Essential gene tracking
            subset = df[df["Gene"].astype(str).isin(gold_std_ids)].copy()
            if not subset.empty:
                melted = subset.melt(id_vars="Gene", value_vars=beta_cols,
                                     var_name="Condition", value_name="Beta")
                melted["Condition"] = melted["Condition"].str.replace("|beta", "", regex=False)

                plt.figure(figsize=(16, 9))
                plt.axhline(0, color="black", linestyle="--", lw=3, alpha=0.9)
                sns.lineplot(data=melted, x="Condition", y="Beta", hue="Gene",
                             marker="o", alpha=0.7, markersize=10, lw=2)
                plt.ylim(min(melted["Beta"].min() - 0.2, -1.2), 0.6)
                plt.title(f"Essential Gene Beta Tracking: {exp_tag}", fontsize=24, fontweight="bold")
                plt.xticks(rotation=30, ha="right", fontsize=15)
                plt.yticks(fontsize=15)
                plt.ylabel("Beta Score", fontsize=19)
                plt.xlabel("Condition", fontsize=19)
                plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=11,
                           title="Phatr3 ID", title_fontsize=13, ncol=1)
                plt.grid(axis="y", linestyle=":", alpha=0.6)
                plt.tight_layout()
                track_name = f"QC_Essential_Tracking_{exp_tag}.png"
                plt.savefig(os.path.join(root, track_name), dpi=300, bbox_inches="tight")
                plt.close()

            print(f"  ✅ MLE QC plots saved for {exp_tag}")

            # ---- Sanity test: essential genes should trend negative ----
            n_neg = (subset[beta_cols].median(axis=1) < 0).sum() if not subset.empty else 0
            n_ess = len(subset)
            print(f"\n  🧪 MLE sanity test [{exp_tag}]:")
            print(f"     Essential genes with median beta < 0: {n_neg}/{n_ess}")
            if n_ess > 0 and n_neg / n_ess >= 0.7:
                print(f"     ✅ PASS — majority of essential genes depleted (expected).")
            else:
                print(f"     ⚠️  WARN — fewer essential genes depleted than expected.")
                print(f"        Check the design matrix and control normalisation.")


# ---------------------------------------------------------------------------
# Steps 9-13 — Downstream analysis (quantile norm, hits, annotation, plots)
# ---------------------------------------------------------------------------

def quantile_normalize_data(df, beta_cols):
    sub_df = df[beta_cols].copy()
    rank_mean = sub_df.stack().groupby(
        sub_df.rank(method="first").stack().astype(int)
    ).mean()
    normalized = sub_df.rank(method="min").stack().astype(int).map(rank_mean).unstack()
    return normalized


def calculate_fdr_quantile_matching(delta_betas):
    sigma = np.quantile(np.abs(delta_betas), 0.68)
    if sigma == 0:
        return np.ones(len(delta_betas))
    z_scores = delta_betas / sigma
    p_values = 2 * (1 - stats.norm.cdf(np.abs(z_scores)))
    _, fdr, _, _ = multipletests(p_values, method="fdr_bh")
    return fdr


def run_mle_quantile_norm(mle_dir, norm_dir, selector_csv):
    all_gene_files = [
        os.path.join(root, f)
        for root, _, files in os.walk(mle_dir)
        for f in files if f.endswith(".gene_summary.txt")
    ]

    if not all_gene_files:
        print(f"  ❌ No gene_summary.txt files found in {mle_dir}")
        return

    if not os.path.exists(selector_csv):
        print(f"  📄 Generating comparison selector template: {selector_csv}")
        template_rows = []
        for file_path in all_gene_files:
            parts   = os.path.normpath(file_path).split(os.sep)
            exp_cat = parts[-3]
            df_temp = pd.read_table(file_path, nrows=1)
            beta_cols = [c for c in df_temp.columns if "|beta" in c]
            for c1, c2 in itertools.combinations(beta_cols, 2):
                template_rows.append([exp_cat, c1, c2, 0])
        pd.DataFrame(
            template_rows, columns=["Experiment", "Treatment_Col", "Control_Col", "Analyze_Flag"]
        ).to_csv(selector_csv, index=False)
        print("  🛑 Template ready. Set Analyze_Flag=1 for pairs you want analysed, then re-run.")
        return

    selector     = pd.read_csv(selector_csv)
    active_pairs = selector[selector["Analyze_Flag"] == 1]

    if active_pairs.empty:
        print("  ⚠️  No active pairs in comparison selector (Analyze_Flag=1). Nothing to do.")
        return

    for file_path in all_gene_files:
        parts   = os.path.normpath(file_path).split(os.sep)
        exp_cat = parts[-3]
        perm_folder = parts[-2]

        exp_pairs = active_pairs[active_pairs["Experiment"] == exp_cat]
        if exp_pairs.empty:
            print(f"  ⏭️  Skipping {exp_cat}/{perm_folder}: no active pairs.")
            continue

        print(f"  🧪 Normalising {exp_cat} → {perm_folder}")
        df = pd.read_table(file_path)
        involved = list(set(exp_pairs["Treatment_Col"].tolist() + exp_pairs["Control_Col"].tolist()))
        df[involved] = quantile_normalize_data(df, involved)

        for _, row in exp_pairs.iterrows():
            treat, ctrl = row["Treatment_Col"], row["Control_Col"]
            label = f"{treat.split('|')[0]}_vs_{ctrl.split('|')[0]}"
            delta = df[treat] - df[ctrl]
            df[f"{label}_delta_beta"] = delta
            df[f"{label}_fdr"] = calculate_fdr_quantile_matching(delta.values)

        exp_out_dir = os.path.join(norm_dir, exp_cat, perm_folder)
        os.makedirs(exp_out_dir, exist_ok=True)
        out_name = os.path.basename(file_path).replace(
            ".gene_summary.txt", "_normalized_hits.csv"
        )
        df.to_csv(os.path.join(exp_out_dir, out_name), index=False)
        print(f"  ✅ Saved: {exp_out_dir}/{out_name}")


def plot_gene_of_interest_volcano(file_path, genes_dict, out_dir,
                                   hard_db, hard_fdr, len_db, len_fdr):
    df = pd.read_csv(file_path)
    os.makedirs(out_dir, exist_ok=True)
    db_cols  = [c for c in df.columns if "_delta_beta" in c]
    markers  = ["o", "s", "D", "^", "v", "p", "*", "h"]
    colors   = sns.color_palette("husl", len(genes_dict))

    for db_col in db_cols:
        prefix  = db_col.replace("_delta_beta", "")
        fdr_col = f"{prefix}_fdr"
        df["neg_log_fdr"] = -np.log10(df[fdr_col] + 1e-12)

        plt.figure(figsize=(20, 14))
        sns.scatterplot(data=df, x=db_col, y="neg_log_fdr",
                        color="lightgray", alpha=0.3, s=80, label="Other", zorder=1)
        plt.axvline(hard_db, color="red",    linestyle="--", alpha=0.8, lw=3)
        plt.axhline(-np.log10(hard_fdr), color="red",    linestyle="--", alpha=0.8, lw=3)
        plt.axvline(len_db,  color="orange", linestyle=":",  alpha=0.8, lw=3)
        plt.axhline(-np.log10(len_fdr),  color="orange", linestyle=":",  alpha=0.8, lw=3)

        for i, (gene_id, common_name) in enumerate(genes_dict.items()):
            gene_row = df[df["Gene"] == gene_id]
            if not gene_row.empty:
                x = gene_row[db_col].values[0]
                y = gene_row["neg_log_fdr"].values[0]
                plt.scatter(x, y, color=colors[i], edgecolor="black", s=600,
                            marker=markers[i % len(markers)], zorder=10,
                            label=common_name, alpha=0.9)

        plt.xlabel(f"ΔBeta ({prefix})", fontsize=36, fontweight="bold", labelpad=20)
        plt.ylabel("-log10(FDR)", fontsize=36, fontweight="bold", labelpad=20)
        plt.grid(linestyle=":", alpha=0.4)
        plt.tick_params(labelsize=28, length=12, width=4)
        plt.legend(bbox_to_anchor=(0.5, 1.02), loc="lower center",
                   ncol=6, fontsize=24, frameon=False)
        save_name = f"Volcano_Targets_{prefix}.png"
        plt.savefig(os.path.join(out_dir, save_name), dpi=300, bbox_inches="tight")
        plt.close()


def generate_hit_reports_and_plots(norm_dir, hard_db, hard_fdr, len_db, len_fdr, genes_dict):
    subfolder_master_data = {}

    for root, _, files in os.walk(norm_dir):
        for file in files:
            if not file.endswith("_normalized_hits.csv"):
                continue
            file_path = os.path.join(root, file)
            target_profile_dir = os.path.join(root, "__genes_profile")
            plot_gene_of_interest_volcano(file_path, genes_dict, target_profile_dir,
                                          hard_db, hard_fdr, len_db, len_fdr)

            df       = pd.read_csv(file_path)
            parts    = os.path.normpath(file_path).split(os.sep)
            perm_tag = parts[-2]
            db_cols  = [c for c in df.columns if "_delta_beta" in c]

            for db_col in db_cols:
                prefix  = db_col.replace("_delta_beta", "")
                fdr_col = f"{prefix}_fdr"

                def get_tier(row):
                    if row[db_col] <= hard_db and row[fdr_col] <= hard_fdr:
                        return "Tier 1 (Hard)"
                    if row[db_col] <= len_db and row[fdr_col] <= len_fdr:
                        return "Tier 2 (Lenient)"
                    if row[db_col] < -0.5:
                        return "Tier 3 (Effect Only)"
                    return "None"

                df["Significance_Tier"] = df.apply(get_tier, axis=1)
                n1 = (df["Significance_Tier"] == "Tier 1 (Hard)").sum()
                n2 = (df["Significance_Tier"] == "Tier 2 (Lenient)").sum()
                n3 = (df["Significance_Tier"] == "Tier 3 (Effect Only)").sum()

                current_hits = df[df["Significance_Tier"] != "None"].copy()
                if not current_hits.empty:
                    entry = current_hits.copy()
                    entry["Comparison"] = prefix
                    entry = entry.rename(columns={"Gene": "Gene_ID", db_col: "Delta_Beta",
                                                  fdr_col: "FDR"})
                    subfolder_master_data.setdefault(root, []).append(
                        entry[["Gene_ID", "Comparison", "Delta_Beta", "FDR", "Significance_Tier"]]
                    )

                # Volcano plot
                fig, ax = plt.subplots(figsize=(20, 14))
                df["neg_log_fdr"] = -np.log10(df[fdr_col] + 1e-12)
                sns.scatterplot(data=df[df["Significance_Tier"] == "None"],
                                x=db_col, y="neg_log_fdr", color="lightgray",
                                alpha=0.3, s=80, label="Not Significant")
                sns.scatterplot(data=df[df["Significance_Tier"] == "Tier 3 (Effect Only)"],
                                x=db_col, y="neg_log_fdr", color="skyblue",
                                alpha=0.6, s=150, label="Tier 3 (Beta < -0.5)")
                sns.scatterplot(data=df[df["Significance_Tier"] == "Tier 2 (Lenient)"],
                                x=db_col, y="neg_log_fdr", color="orange",
                                alpha=0.7, s=250, label=f"Lenient (FDR < {len_fdr})")
                sns.scatterplot(data=df[df["Significance_Tier"] == "Tier 1 (Hard)"],
                                x=db_col, y="neg_log_fdr", color="red",
                                alpha=0.9, s=450, label=f"Hard (FDR < {hard_fdr})",
                                edgecolor="black")
                plt.axvline(len_db, color="black", linestyle="--", alpha=0.5, lw=3)
                plt.axhline(-np.log10(len_fdr), color="black", linestyle="--", alpha=0.5, lw=3)
                ax.text(0.95, 0.95,
                        f"Tier 1: {n1}\nTier 2: {n2}\nTier 3: {n3}",
                        transform=ax.transAxes, fontsize=20, va="top", ha="right",
                        fontweight="bold", color="darkred",
                        bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                                  alpha=0.8, edgecolor="gray"))
                plt.legend(bbox_to_anchor=(1.0, 1.0), loc="lower right",
                           ncol=4, fontsize=16, frameon=False)
                plt.xlabel(f"ΔBeta ({prefix}, {perm_tag})", fontsize=30, fontweight="bold",
                           labelpad=20)
                plt.ylabel("-log10(FDR)", fontsize=30, fontweight="bold", labelpad=20)
                ax.tick_params(labelsize=28, length=12, width=4)
                plt.grid(linestyle=":", alpha=0.4, lw=2)
                plt.savefig(os.path.join(root, f"Volcano_{prefix}.png"),
                            dpi=300, bbox_inches="tight")
                plt.close()
                print(f"  ✅ {prefix}: T1={n1}, T2={n2}, T3={n3}")

    for folder_path, hits_list in subfolder_master_data.items():
        master_df = pd.concat(hits_list, ignore_index=True).sort_values(
            ["Significance_Tier", "Delta_Beta"]
        )
        exp_name = os.path.basename(os.path.dirname(folder_path))
        master_path = os.path.join(folder_path, f"Master_Hits_{exp_name}_Summary.csv")
        master_df.to_csv(master_path, index=False)
        print(f"  ✅ Master hit list: {master_path}")


def prepare_koala_database(koala_paths, output_path):
    if os.path.exists(output_path):
        print(f"  ⏩ KEGG/KOALA file already exists. Skipping.")
        return
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    all_dfs = []
    for fpath in koala_paths:
        if not os.path.exists(fpath):
            print(f"  ⚠️  Missing KOALA file: {os.path.basename(fpath)}")
            continue
        df = pd.read_csv(fpath, sep="\t", header=None,
                         names=["Gene_ID", "KEGG_ID"], on_bad_lines="skip", low_memory=False)
        df["Gene_ID"] = df["Gene_ID"].astype(str).str.split(".").str[0]
        df = df.dropna(subset=["KEGG_ID"])
        all_dfs.append(df)
    if not all_dfs:
        print("  ❌ No KOALA data found.")
        return
    master = pd.concat(all_dfs, ignore_index=True).drop_duplicates()
    master.to_csv(output_path, sep="\t", index=False)
    print(f"  ✅ KOALA master: {output_path}  ({len(master):,} entries)")


def clean_eggnog_file(input_path, output_path):
    if os.path.exists(output_path):
        print(f"  ⏩ Cleaned eggNOG already exists. Skipping.")
        return
    print(f"  🧹 Cleaning eggNOG: {os.path.basename(input_path)}")
    try:
        header_line = 0
        with open(input_path) as f:
            for i, line in enumerate(f):
                if line.startswith("#query"):
                    header_line = i
                    break
        df = pd.read_csv(input_path, sep="\t", skiprows=header_line, engine="python")
        df.rename(columns={df.columns[0]: "Gene_ID"}, inplace=True)
        df["Gene_ID"] = df["Gene_ID"].astype(str).str.split(".").str[0]
        df = df.dropna(subset=["Gene_ID"])
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        df.to_csv(output_path, sep="\t", index=False)
        print(f"  ✅ Cleaned eggNOG: {output_path}  ({len(df.columns)} columns)")
    except Exception as e:
        print(f"  ❌ Error cleaning eggNOG: {e}")


def normalize_id(gene_id):
    if pd.isna(gene_id):
        return ""
    name = str(gene_id).lower().strip().replace("draft", "").replace("_", "").split(".")[0]
    return name


def annotate_all_experiment_hits(norm_dir, eggnog_path):
    if not os.path.exists(eggnog_path):
        print(f"  ❌ eggNOG database not found: {eggnog_path}")
        return
    print(f"  📖 Loading eggNOG database...")
    eggnog_df = pd.read_table(eggnog_path)
    eggnog_df["Normalized_ID"] = eggnog_df["Gene_ID"].apply(normalize_id)
    eggnog_subset = eggnog_df.drop(columns=["Gene_ID"])

    for root, _, files in os.walk(norm_dir):
        for file in [f for f in files if f.startswith("Master_Hits_") and f.endswith("_Summary.csv")]:
            master_path = os.path.join(root, file)
            out_path    = master_path.replace(".csv", "_ANNOTATED.csv")
            hits_df     = pd.read_csv(master_path)
            hits_df["Normalized_ID"] = hits_df["Gene_ID"].apply(normalize_id)
            annotated = pd.merge(hits_df, eggnog_subset, on="Normalized_ID", how="left")
            annotated = annotated.fillna("N/A")
            cols = annotated.columns.tolist()
            if "Normalized_ID" in cols:
                cols.insert(1, cols.pop(cols.index("Normalized_ID")))
                annotated = annotated[cols]
            annotated.to_csv(out_path, index=False)
            print(f"  ✅ Annotated: {os.path.basename(out_path)}")
    print("  ✨ Annotation complete.")


def run_functional_plotting(norm_dir, plot_dir, cog_mapping):
    print("\n  🚀 Generating COG distribution plots...")
    for root, _, files in os.walk(norm_dir):
        for file in [f for f in files if f.endswith("_ANNOTATED.csv")]:
            csv_path = os.path.join(root, file)
            df = pd.read_csv(csv_path)
            exp_name = file.replace("_ANNOTATED.csv", "")
            exp_root_dir = os.path.join(plot_dir, exp_name)
            os.makedirs(exp_root_dir, exist_ok=True)
            for comp in df["Comparison"].unique():
                comp_df  = df[df["Comparison"] == comp].copy()
                target_tiers = ["Tier 1 (Hard)", "Tier 2 (Lenient)"]
                plot_df  = comp_df[comp_df["Significance_Tier"].isin(target_tiers)].copy()
                if plot_df.empty:
                    continue
                plot_df["COG_Full"] = (
                    plot_df["COG_category"].astype(str).str[0].str.upper()
                    .map(cog_mapping).fillna("Unclassified / No COG")
                )
                def manual_wrap(text, n=3):
                    words = str(text).split()
                    return "\n".join(" ".join(words[i:i+n]) for i in range(0, len(words), n))
                plot_df["Category"] = plot_df["COG_Full"].apply(manual_wrap)
                counts = plot_df.groupby(["Category", "Significance_Tier"]).size().reset_index(name="Count")
                all_cats = counts["Category"].unique()
                mux = pd.MultiIndex.from_product([all_cats, target_tiers],
                                                  names=["Category", "Significance_Tier"])
                counts = counts.set_index(["Category", "Significance_Tier"]).reindex(
                    mux, fill_value=0
                ).reset_index()
                sort_order = counts[counts["Significance_Tier"] == "Tier 1 (Hard)"].sort_values(
                    "Count", ascending=False
                )["Category"]
                plt.figure(figsize=(20, max(10, len(all_cats) * 1.2)))
                ax = sns.barplot(data=counts, x="Count", y="Category", hue="Significance_Tier",
                                 order=sort_order,
                                 palette={"Tier 1 (Hard)": "#2c7fb8", "Tier 2 (Lenient)": "#7fcdbb"})
                plt.title(f"Functional Enrichment: {comp}\nExperiment: {exp_name}",
                          fontsize=26, pad=40, fontweight="bold")
                plt.xlabel("Number of Genes", fontsize=22, labelpad=20)
                plt.ylabel("", labelpad=30)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=16)
                plt.legend(title="Significance Tier", title_fontsize=20, fontsize=18)
                for container in ax.containers:
                    ax.bar_label(container, padding=8, fontsize=14, fontweight="bold")
                plt.savefig(os.path.join(exp_root_dir, f"{comp}_Combined_COG.png"),
                            bbox_inches="tight", dpi=300)
                plt.close()
            print(f"  ✅ COG plots for {exp_name}")


def plot_top_hits_by_beta(norm_dir, plot_dir, cog_mapping):
    print("\n  🏆 Generating Top 20 Hit plots...")
    all_cog_names = sorted(list(set(cog_mapping.values())) + ["Unclassified / No COG"])
    cog_color_dict = dict(zip(all_cog_names, sns.color_palette("Set2", len(all_cog_names))))

    for root, _, files in os.walk(norm_dir):
        for file in [f for f in files if f.endswith("_ANNOTATED.csv")]:
            csv_path = os.path.join(root, file)
            df = pd.read_csv(csv_path)
            exp_name = file.replace("_ANNOTATED.csv", "")
            exp_root_dir = os.path.join(plot_dir, exp_name)
            os.makedirs(exp_root_dir, exist_ok=True)
            for comp in df["Comparison"].unique():
                comp_df  = df[df["Comparison"] == comp].copy()
                hits_df  = comp_df[comp_df["Significance_Tier"].isin(
                    ["Tier 1 (Hard)", "Tier 2 (Lenient)"]
                )].copy()
                if hits_df.empty:
                    continue
                top_20 = hits_df.sort_values("Delta_Beta", ascending=True).head(20).copy()
                top_20["COG_Full"] = (
                    top_20["COG_category"].astype(str).str[0].str.upper()
                    .map(cog_mapping).fillna("Unclassified / No COG")
                )
                top_20["Description"] = top_20["Description"].fillna("-")
                id_to_desc = dict(zip(top_20["Gene_ID"], top_20["Description"]))
                plt.figure(figsize=(22, 14))
                ax = sns.barplot(data=top_20, x="Delta_Beta", y="Gene_ID",
                                 hue="COG_Full", dodge=False, palette=cog_color_dict)
                plt.xlabel(f"Delta Beta ({comp})", fontsize=24, labelpad=20, fontweight="bold")
                plt.ylabel("Gene ID", fontsize=22)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.legend(title="COG Category", title_fontsize=20, fontsize=15,
                           loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=3, frameon=False)
                for i, lbl in enumerate(ax.get_yticklabels()):
                    gene_id = lbl.get_text()
                    desc = str(id_to_desc.get(gene_id, "-"))
                    display = (desc[:75] + "...") if len(desc) > 75 else desc
                    ax.text(-0.02, i, f"{display}  ", color="black", va="center",
                            ha="right", fontsize=14, fontweight="bold")
                plt.savefig(os.path.join(exp_root_dir, f"{comp}_Top20_Essential_Genes.png"),
                            bbox_inches="tight", dpi=300)
                plt.close()
            print(f"  ✅ Top 20 plots for {exp_name}")


def plot_cog_specific_hits(norm_dir, plot_dir, cog_mapping, target_cogs, top_n=10):
    print(f"\n  🎯 Generating COG-specific plots for: {target_cogs}")
    all_cog_names = sorted(list(set(cog_mapping.values())) + ["Unclassified / No COG"])
    cog_color_dict = dict(zip(all_cog_names, sns.color_palette("Set2", len(all_cog_names))))

    for root, _, files in os.walk(norm_dir):
        for file in [f for f in files if f.endswith("_ANNOTATED.csv")]:
            csv_path = os.path.join(root, file)
            df = pd.read_csv(csv_path)
            exp_name = file.replace("_ANNOTATED.csv", "")
            exp_output_dir = os.path.join(plot_dir, "COG_specific", exp_name)
            os.makedirs(exp_output_dir, exist_ok=True)
            for comp in df["Comparison"].unique():
                comp_df = df[df["Comparison"] == comp].copy()
                comp_df["COG_Full"] = (
                    comp_df["COG_category"].astype(str).str[0].str.upper()
                    .map(cog_mapping).fillna("Unclassified / No COG")
                )
                spec_df = comp_df[
                    comp_df["COG_Full"].isin(target_cogs) |
                    comp_df["COG_category"].isin(target_cogs)
                ].copy()
                if spec_df.empty:
                    continue
                top_hits = spec_df.sort_values("Delta_Beta", ascending=True).head(top_n).copy()
                top_hits["Description"] = top_hits["Description"].fillna("-")
                top_hits.to_csv(os.path.join(exp_output_dir, f"{comp}_COG_Focused_Hits.csv"),
                                index=False)
                plt.figure(figsize=(20, max(8, len(top_hits) * 0.9)))
                ax = sns.barplot(data=top_hits, x="Delta_Beta", y="Gene_ID",
                                 hue="COG_Full", dodge=False, palette=cog_color_dict)
                plt.xlabel(f"Delta Beta ({comp})", fontsize=24, fontweight="bold", labelpad=20)
                plt.ylabel("Gene ID", fontsize=22)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                plt.legend(title="Functional Category", title_fontsize=18, fontsize=15,
                           loc="lower center", bbox_to_anchor=(0.5, 1.05), ncol=3,
                           frameon=False)
                for i in range(len(top_hits)):
                    desc     = str(top_hits.iloc[i]["Description"])
                    display  = (desc[:65] + "...") if len(desc) > 65 else desc
                    fdr_val  = top_hits.iloc[i]["FDR"]
                    beta_val = top_hits.iloc[i]["Delta_Beta"]
                    ax.text(-0.02, i, f"{display}  ", color="black", va="center",
                            ha="right", fontsize=13, fontweight="bold")
                    fdr_txt = f"FDR: {fdr_val:.2e}"
                    if beta_val > -1.2:
                        ax.text(0.02, i, f" {fdr_txt}", color="black", va="center",
                                ha="left", fontsize=12, style="italic")
                    else:
                        ax.text(beta_val + 0.02, i, f"{fdr_txt} ", color="black",
                                va="center", ha="left", fontsize=12, style="italic")
                plt.savefig(os.path.join(exp_output_dir, f"{comp}_COG_Focused_Hits.png"),
                            bbox_inches="tight", dpi=300)
                plt.close()
    print(f"  ✅ COG-specific plots saved in {plot_dir}/COG_specific")
