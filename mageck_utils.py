import re
import os
import sys
import gzip
import time
import math
import warnings
import subprocess
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

warnings.simplefilter(action='ignore', category=FutureWarning)


def check_env():
    """Checks for required Python packages and MAGeCK."""
    packages = {"Bio": "biopython", "pandas": "pandas", "tqdm": "tqdm"}
    missing = []
    for module, install_name in packages.items():
        try:
            __import__(module)
        except ImportError:
            missing.append(install_name)
    try:
        subprocess.run(["mageck", "-v"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (FileNotFoundError, TypeError):
        missing.append("mageck")
    return missing

def initialize_workspace(input_path, output_path, seq_parent, seq_subfolder, lib_path):
    """Creates directory structure; returns True if folders were newly made."""
    # We create the parent '2_sequencing' and the specific subfolder
    raw_data_path = os.path.join(seq_parent, seq_subfolder)
    dirs = [input_path, output_path, seq_parent, raw_data_path, lib_path]
    newly_created = False
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
            print(f"📁 Created: {d}")
            newly_created = True
    return newly_created

def generate_metadata_template(csv_path, raw_dir):
    """Creates a blank metadata CSV based on files in the RAW_DATA_DIR."""
    if os.path.exists(csv_path):
        print(f"✅ Metadata file exists: {csv_path}. (No overwrite)")
        return

    if not os.path.exists(raw_dir):
        print(f"⚠️ Error: Sequencing subfolder not found at: {raw_dir}")
        return

    # Filter out 'Undetermined' files and non-fastq files
    files = sorted([
        f for f in os.listdir(raw_dir) 
        if f.endswith(".fastq.gz") and "Undetermined" not in f
    ])
    
    if not files:
        print(f"⚠️ No valid .fastq.gz files found in {raw_dir}")
        return

    df = pd.DataFrame({
        'filename': files,
        'sample_label': [f.split('_S')[0] for f in files],
        'experiment_label': '', # manually duplicate rows for t0 if there are multiple experiments with same t0
        'condition': '',      
        'bio_rep': ''        # NEW: (e.g., 1, 2)

    })
    df.to_csv(csv_path, index=False)
    print(f"🌟 Template created successfully at: {csv_path}")

def validate_metadata(csv_path):
    """
    Checks if metadata is filled and provides an experimental design summary.
    Ensures experiment_label, condition, and bio_rep are present.
    """
    if not os.path.exists(csv_path):
        return False, f"Metadata file missing at {csv_path}."
    
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        return False, f"Could not read CSV: {e}"
    
    # 1. Check for missing values in required columns
    required_cols = ['experiment_label', 'condition', 'bio_rep']
    
    # Ensure columns exist first
    if not all(col in df.columns for col in required_cols):
        missing_cols = [col for col in required_cols if col not in df.columns]
        return False, f"Metadata is missing required columns: {missing_cols}"

    # Check for empty/NaN cells
    for col in required_cols:
        empty_count = df[col].isna().sum() + (df[col].astype(str).str.strip() == '').sum()
        if empty_count > 0:
            return False, f"Column '{col}' has {empty_count} empty entries. Please fill them."

    # 2. Generate Summary of Biological Replicates
    # We group by Experiment and Condition, then count the number of files (rows) assigned.
    # Note: If technical replicates exist, this counts total files per condition.
    summary = df.groupby(['experiment_label', 'condition'])['bio_rep'].nunique().unstack(fill_value=0)
    
    print("\n" + "="*60)
    print("📊 EXPERIMENTAL DESIGN SUMMARY")
    print("Detected Biological Replicates per Condition:")
    print("="*60)
    print(summary)
    print("="*60)
    print("Check if your shared t0 appears in all relevant experiments.")
    print("="*60 + "\n")
    
    return True, "Metadata verification successful."

def _extract_worker(row_tuple, raw_dir, out_dir, fwd, rev, length):
    """Internal helper to process one file with unique naming and QC metrics."""
    _, row = row_tuple
    in_p = os.path.join(raw_dir, row['filename'])
    
    # Strip extensions to create a unique identifier from the original filename
    clean_name = row['filename'].split('.fast')[0].split('.fq')[0]
    out_filename = f"{row['experiment_label']}_{row['condition']}_R{row['bio_rep']}_{clean_name}.fastq.gz"
    out_p = os.path.join(out_dir, out_filename)
    
    total, extracted = 0, 0
    try:
        with gzip.open(in_p, "rt") as in_h, gzip.open(out_p, "wt") as out_h:
            for record in SeqIO.parse(in_h, "fastq"):
                total += 1
                seq_str = str(record.seq).upper()
                extracted_record = None
                
                if fwd in seq_str:
                    start_idx = seq_str.find(fwd) + len(fwd)
                    if (len(seq_str) - start_idx) >= length:
                        extracted_record = record[start_idx : start_idx + length]
                elif rev in seq_str:
                    start_idx = seq_str.find(rev) + len(rev)
                    if (len(seq_str) - start_idx) >= length:
                        raw_chunk = record[start_idx : start_idx + length]
                        extracted_record = raw_chunk.reverse_complement(id=True, description=True)

                if extracted_record:
                    SeqIO.write(extracted_record, out_h, "fastq")
                    extracted += 1
        
        # Calculate percentage for this specific file
        pct = (extracted / total * 100) if total > 0 else 0
        
        return {
            'Success': True, 
            'Experiment': row['experiment_label'],
            'Condition': row['condition'],
            'Bio_Rep': row['bio_rep'],
            'Original_File': row['filename'], 
            'Output_File': out_filename, 
            'Total_Reads': total, 
            'Extracted_Reads': extracted,
            'Success_Rate_Pct': round(pct, 2)
        }
    except Exception as e:
        return {'Success': False, 'Output_File': out_filename, 'Error': str(e)}

def run_guide_extraction(meta_df, raw_dir, out_dir, fwd, rev, length, n_cores):
    """Parallel extraction with detailed QC summary."""
    os.makedirs(out_dir, exist_ok=True)
    print(f"🚀 Launching parallel extraction on {n_cores} cores...")

    worker_task = partial(_extract_worker, raw_dir=raw_dir, out_dir=out_dir, 
                          fwd=fwd, rev=rev, length=length)
    
    with mp.Pool(processes=n_cores) as pool:
        results = list(tqdm(pool.imap(worker_task, meta_df.iterrows()), 
                           total=len(meta_df), desc="Extraction Progress"))

    # Filter out successful results and create DataFrame
    success_results = [r for r in results if r.get('Success')]
    summary_df = pd.DataFrame(success_results)
    
    # Drop the internal 'Success' flag before saving
    if not summary_df.empty:
        summary_df = summary_df.drop(columns=['Success'])
        summary_df.to_csv(os.path.join(out_dir, "extraction_summary.csv"), index=False)
        
        print("\n" + "="*85)
        print(f"{'Output File':<40} | {'Total':<12} | {'Extracted':<12} | {'Success %':<10}")
        print("-" * 85)
        for _, r in summary_df.head(10).iterrows(): # Show first 10 for brevity in terminal
            print(f"{r['Output_File'][:40]:<40} | {r['Total_Reads']:<12} | {r['Extracted_Reads']:<12} | {r['Success_Rate_Pct']:<10}%")
        print("="*85)
        print(f"✅ Full summary saved to: {out_dir}/extraction_summary.csv")
    
    # Report errors
    errors = [r for r in results if not r.get('Success')]
    for err in errors:
        print(f"❌ Error in {err['Output_File']}: {err.get('Error')}")

def _generate_count_plots(exp, count_dir):
    """
    Internal helper to plot QC metrics. 
    Aggregates technical replicates to avoid FixedLocator ValueErrors.
    """
    
    # Matching your count output naming: {exp}_.count.txt
    summary_file = os.path.join(count_dir, f"{exp}_.countsummary.txt")
    count_file = os.path.join(count_dir, f"{exp}_.count.txt")
    
    print(f"📊 Attempting to plot: {summary_file}")
    
    if not (os.path.exists(summary_file) and os.path.exists(count_file)):
        print(f"⚠️ File Not Found! Check if this exists: {summary_file}")
        return

    # Load data
    summary_df = pd.read_table(summary_file)
    counts = pd.read_table(count_file)
    
    # --- FIX: Aggregate Technical Replicates ---
    # This collapses the 40 original rows into 10 biological rows
    plot_summary = summary_df.groupby('Label', sort=False).agg({
        'Percentage': 'mean', 
        'Zerocounts': 'mean'
    }).reset_index()
    
    # --- Prepare Data for Violin Plot ---
    numeric_cols = [col for col in counts.columns if col not in ['sgRNA', 'Gene']]
    melted_counts = counts.melt(id_vars=['sgRNA'], value_vars=numeric_cols, var_name='Label', value_name='Count')
    melted_counts['Log10_Count'] = np.log10(melted_counts['Count'] + 1)
    
    # Maintain consistent sample order for the X-axis
    sample_order = sorted(melted_counts['Label'].unique())

    # --- Create Figure ---
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    fig.suptitle(f"QC Summary: {exp}", fontsize=18, y=1.05)

    # 1. Mapping Efficiency (0-100%)
    sns.barplot(x=plot_summary['Label'], y=plot_summary['Percentage'] * 100, ax=axes[0], palette='viridis')
    axes[0].set_ylabel('Percentage (%)', fontsize=12)
    axes[0].set_ylim(0, 100)
    axes[0].set_xlabel('', fontsize=0)
    
    # Explicitly set 10 ticks for the 10 aggregated biological labels
    axes[0].set_xticks(range(len(plot_summary)))
    axes[0].set_xticklabels(plot_summary['Label'], rotation=45, ha='right')
    axes[0].set_title("Mapping Efficiency", fontsize=14)

    # 2. Log10 Count Distribution (Violin)
    sns.violinplot(data=melted_counts, x='Label', y='Log10_Count', order=sample_order, ax=axes[1], palette='muted', inner='quartile')
    axes[1].set_xlabel('', fontsize=0)
    axes[1].set_ylabel('Log10 Count', fontsize=12)
    
    # Explicitly set ticks to match the 10 biological samples
    axes[1].set_xticks(range(len(sample_order)))
    axes[1].set_xticklabels(sample_order, rotation=45, ha='right')
    axes[1].set_title("Read Distribution (Log10)", fontsize=14)

    # 3. Zero Counts
    sns.barplot(x=plot_summary['Label'], y=plot_summary['Zerocounts'], ax=axes[2], palette='magma')
    axes[2].set_ylabel('Count of Zero-Reads sgRNAs', fontsize=12)
    axes[2].set_xlabel('', fontsize=0)
    
    # Explicitly set 10 ticks for the 10 labels
    axes[2].set_xticks(range(len(plot_summary)))
    axes[2].set_xticklabels(plot_summary['Label'], rotation=45, ha='right')
    axes[2].set_title("Zero-Count sgRNAs", fontsize=14)

    # --- Save and Close ---
    plt.tight_layout()
    save_path = os.path.join(count_dir, f"{exp}_QC_Summary.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"✅ Successfully saved QC plot to: {save_path}")

def _generate_control_qc_plots(exp, count_dir, neg_path, pos_path):
    """
    Plots Raw vs Normalized count distributions for controls in a 2x3 grid.
    Updated for filename: {exp}_.count_normalized.txt
    """
    raw_file = os.path.join(count_dir, f"{exp}_.count.txt")
    norm_file = os.path.join(count_dir, f"{exp}_.count_normalized.txt") # Updated naming
    
    if not (os.path.exists(raw_file) and os.path.exists(norm_file)):
        print(f"⚠️ Missing files for {exp}. Check: {raw_file} or {norm_file}")
        return

    # 1. Load control IDs
    with open(neg_path, 'r') as f:
        neg_ids = [line.strip() for line in f if line.strip()]
    
    pos_df = pd.read_csv(pos_path)
    pos_ids = pos_df.iloc[:, 1].tolist()

    def get_type(sgrna):
        if sgrna in pos_ids: return 'Essential (Pos)'
        if sgrna in neg_ids:
            return 'Non-targeting' if 'nontargetting' in sgrna.lower() else 'Safe (Neg)'
        return 'Experimental'

    # 2. Process both Raw and Norm Data
    dfs = []
    for fpath, label in [(raw_file, 'Raw'), (norm_file, 'Normalized')]:
        df = pd.read_table(fpath)
        df['Type'] = df['sgRNA'].apply(get_type)
        df['Scale'] = label
        
        # Identify sample columns (exclude metadata)
        numeric_cols = [c for c in df.columns if c not in ['sgRNA', 'Gene', 'Type', 'Scale']]
        melted = df.melt(id_vars=['sgRNA', 'Type', 'Scale'], value_vars=numeric_cols, 
                         var_name='Sample', value_name='Count')
        
        # Use log10 for better visualization of count distribution
        melted['Log10_Count'] = np.log10(melted['Count'] + 1)
        dfs.append(melted[melted['Type'] != 'Experimental'])

    full_plot_df = pd.concat(dfs)
    counts_summary = pd.read_table(raw_file)['sgRNA'].apply(get_type).value_counts()

    # --- 3. Create 2x3 Figure ---
    fig, axes = plt.subplots(2, 3, figsize=(22, 14), sharey='row')
    fig.suptitle(f"Control Guide Behavior (Raw vs Normalized): {exp}", fontsize=22, y=1.02)

    colors = ['#4C72B0', '#55A868', '#C44E52']
    categories = ['Safe (Neg)', 'Non-targeting', 'Essential (Pos)']

    for row_idx, scale_label in enumerate(['Raw', 'Normalized']):
        for col_idx, cat in enumerate(categories):
            ax = axes[row_idx, col_idx]
            subset = full_plot_df[(full_plot_df['Scale'] == scale_label) & (full_plot_df['Type'] == cat)]
            
            if subset.empty:
                ax.text(0.5, 0.5, f"No {cat} found", ha='center')
                continue
                
            sns.violinplot(data=subset, x='Sample', y='Log10_Count', ax=ax, 
                           color=colors[col_idx], inner='quartile', cut=0)
            
            # Row-specific Titles
            if row_idx == 0:
                n_guides = counts_summary.get(cat, 0)
                ax.set_title(f"RAW: {cat}\n(n={n_guides})", fontsize=15, fontweight='bold')
            else:
                ax.set_title(f"NORM: {cat}", fontsize=15, fontweight='bold')

            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            ax.set_xlabel('')
            
            if col_idx == 0:
                ax.set_ylabel(f"{scale_label} Log10(Count+1)", fontsize=13)
            else:
                ax.set_ylabel("")

    plt.tight_layout()
    save_path = os.path.join(count_dir, f"{exp}_Control_QC.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"✅ Success: Comparison plot saved for {exp}")

def _plot_specific_gene_trajectories(exp, count_dir, genes_dict):
    """
    Plots Mean Trajectory (with Shaded StDev) for Raw vs Normalized counts.
    STRICT string matching for Gene IDs.
    """
    raw_file = os.path.join(count_dir, f"{exp}_.count.txt")
    norm_file = os.path.join(count_dir, f"{exp}_.count_normalized.txt")
    
    if not (os.path.exists(raw_file) and os.path.exists(norm_file)):
        print(f"⚠️ Skipping trajectories for {exp}: Files not found.")
        return

    gene_out_dir = os.path.join(count_dir, f"{exp}_gene_profiles")
    os.makedirs(gene_out_dir, exist_ok=True)

    # 1. Load Data
    raw_df = pd.read_table(raw_file)
    norm_df = pd.read_table(norm_file)
    sample_cols = [c for c in raw_df.columns if c not in ['sgRNA', 'Gene']]
    
    for gene_id, common_name in genes_dict.items():
        # STRICT MATCHING
        raw_gene = raw_df[raw_df['Gene'] == gene_id]
        norm_gene = norm_df[norm_df['Gene'] == gene_id]

        if raw_gene.empty or norm_gene.empty:
            print(f"⚠️ WARNING: Gene ID '{gene_id}' ({common_name}) not found.")
            continue

        # Create Plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
        plt.subplots_adjust(wspace=0.2)

        # Helper function to process and plot on specific axis
        def plot_agg_trajectory(data_df, ax, title, ylabel):
            # Melt to long format
            melted = data_df.melt(id_vars=['sgRNA'], value_vars=sample_cols, 
                                 var_name='Sample', value_name='Count')
            melted['Log10_Count'] = np.log10(melted['Count'] + 1)
            
            # Using seaborn lineplot with aggregation: 
            # It automatically calculates Mean (line) and StDev (shaded area)
            sns.lineplot(data=melted, x='Sample', y='Log10_Count', ax=ax, 
                         color='royalblue', marker='o', errorbar='sd', linewidth=3)
            
            ax.set_title(title, fontsize=18, fontweight='bold')
            ax.set_ylabel(ylabel, fontsize=14)
            ax.tick_params(axis='x', rotation=45)
            ax.grid(True, alpha=0.3)

        # Plot Raw (Left) and Normalized (Right)
        plot_agg_trajectory(raw_gene, ax1, f"RAW: {common_name}", "Log10(Raw Count + 1)")
        plot_agg_trajectory(norm_gene, ax2, f"NORM: {common_name}", "Log10(Norm Count + 1)")

        fig.suptitle(f"Gene Mean Trajectory (±SD): {common_name} ({gene_id}) | Exp: {exp}", 
                     fontsize=22, y=1.05, fontweight='bold')
        
        save_path = os.path.join(gene_out_dir, f"{common_name}_Mean_Trajectory.png")
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

    print(f"✅ Success: Mean profiles saved in {gene_out_dir}")
       
def run_mageck_count(preprocess_dir, count_dir, library_csv, summary_csv, neg_path, pos_path, genes_dict):
    """
    Runs MAGeCK count per Experiment and saves results to 2_mageck_count.
    Groups tech reps with commas and biological samples with spaces.
    """
    # Create the new output directory
    os.makedirs(count_dir, exist_ok=True)
    
    if not os.path.exists(summary_csv):
        print(f"❌ Error: {summary_csv} not found.")
        return

    df = pd.read_csv(summary_csv)
    experiments = df['Experiment'].unique()
    
    for exp in experiments:
        exp_df = df[df['Experiment'] == exp].sort_values(['Condition', 'Bio_Rep'])
        
        sample_labels = []
        fastq_args = []
        
        # Group by biological replicate (e.g., Nitrogen_t0_R1)
        bio_groups = exp_df.groupby(['Condition', 'Bio_Rep'], sort=False)
        
        for (cond, rep), group in bio_groups:
            label = f"{cond}_R{rep}"
            sample_labels.append(label)
            
            # Join technical replicates with commas (e.g., L001.fq,L002.fq)
            tech_reps = [os.path.join(preprocess_dir, f) for f in group['Output_File'].tolist()]
            fastq_args.append(",".join(tech_reps))
            
        # MAGeCK format: labels separated by commas, fastq groups by spaces
        label_str = ",".join(sample_labels)
        output_prefix = os.path.join(count_dir, f"{exp}_")
        
        # Build the Command
        cmd = [
            "mageck", "count",
            "-l", library_csv,
            "-n", output_prefix,
            "--sample-label", label_str,
            "--fastq"
        ] + fastq_args # Appending the list adds them with spaces
        
        # PRINT FOR VERIFICATION
        print(f"\n🚀 RUNNING MAGECK FOR EXPERIMENT: {exp}")
        print(f"COMMAND:\n{' '.join(cmd)}\n")
        
        try:
            subprocess.run(cmd, check=True)
            print(f"✅ Success! Count table saved to: {output_prefix}.count.txt")
        except subprocess.CalledProcessError as e:
            print(f"❌ MAGeCK Count failed for {exp}: {e}")

    # Plots for Count QC
    df = pd.read_csv(summary_csv)
    for exp in df['Experiment'].unique():
        print(f"🎨 Generating QC plots for {exp}...")
        _generate_count_plots(exp, count_dir)
        _generate_control_qc_plots(exp, count_dir, neg_path, pos_path) # New Control QC
        _plot_specific_gene_trajectories(exp, count_dir, genes_dict)

def manage_all_design_matrices(summary_csv, design_base_dir):
    """
    Manages design matrices for all experiments, ensuring 't0' is 
    excluded from columns and remains represented by the baseline.
    """
    if not os.path.exists(summary_csv):
        print(f"❌ Error: {summary_csv} not found.")
        return False

    df = pd.read_csv(summary_csv)
    unique_exps = df['Experiment'].unique()
    all_ready = True
    
    print(f"\n🔍 Generalized Pipeline: Scanning {len(unique_exps)} experiments...")

    for exp_name in unique_exps:
        exp_dir = os.path.join(design_base_dir, exp_name)
        design_path = os.path.join(exp_dir, f"design_matrix_{exp_name}.txt")
        os.makedirs(exp_dir, exist_ok=True)

        exp_df = df[df['Experiment'] == exp_name]
        
        # Scenario 1: Create Template (Excluding t0)
        if not os.path.exists(design_path):
            labels = exp_df.apply(lambda x: f"{x['Condition']}_R{x['Bio_Rep']}", axis=1).unique()
            
            # Filter out 't0' from the condition columns
            unique_conditions = [c for c in exp_df['Condition'].unique() if c.lower() != 't0']
            
            design_df = pd.DataFrame({'Samples': labels})
            design_df['baseline'] = 1 
            
            for cond in unique_conditions:
                design_df[cond] = "" # User fills with 0 or 1
            
            design_df.to_csv(design_path, sep='\t', index=False)
            
            print(f"📄 Skeleton created: {exp_name} -> {design_path}")
            print(f"⚠️  WARNING: 't0' has been excluded from columns. It is represented by the 'baseline' column.")
            print(f"⚠️  If you manually add columns, DO NOT include a 't0' column; it should remain as all 0s in the treatment columns.")
            all_ready = False
            continue

        # Scenario 2: Validate existing file
        try:
            user_df = pd.read_table(design_path, dtype=str).fillna("")
            cols = user_df.columns.tolist()

            # Extra Safety: Check if the user manually added a t0 column
            if any(c.lower() == 't0' for c in cols):
                print(f"❌ Error in {exp_name}: Found a 't0' condition column.")
                print("👉 Please delete the 't0' column. The 'baseline' column already handles the t0 state.")
                all_ready = False
                continue

            data_cols = cols[2:] 
            is_binary = user_df[data_cols].applymap(lambda x: str(x).strip() in ['0', '1']).all().all()
            has_content = (user_df[data_cols] != "").all().all()

            if not (is_binary and has_content):
                print(f"⚠️  {exp_name}: Matrix incomplete or invalid values.")
                all_ready = False
            else:
                print(f"✅ {exp_name}: Validated ({', '.join(data_cols)}).")

        except Exception as e:
            print(f"❌ {exp_name}: Error reading matrix: {e}")
            all_ready = False

    return all_ready

def run_all_mageck_mle(count_dir, mle_dir, design_base_dir, n_cores, perm_round, max_sgrna_perm, control_file):
    """
    Runs MAGeCK MLE in parameter-specific subfolders with unique file prefixes.
    Example: 4_mageck_mle/Nitrogen/perm10_max15/Nitrogen_mle_perm10_max15.gene_summary.txt
    """
    # Identify experiments
    experiments = [d for d in os.listdir(design_base_dir) 
                  if os.path.isdir(os.path.join(design_base_dir, d))]

    for exp in experiments:
        # 1. Create unique param string
        param_suffix = f"perm{perm_round}_max{max_sgrna_perm}"
        
        # 2. Define directory and prefix
        run_dir = os.path.join(mle_dir, exp, param_suffix)
        os.makedirs(run_dir, exist_ok=True)
        
        # Filename includes experiment, step, and parameters
        output_prefix = os.path.join(run_dir, f"{exp}_mle_{param_suffix}")
        gene_summary = f"{output_prefix}.gene_summary.txt"

        # --- SAFETY CHECK ---
        if os.path.exists(gene_summary):
            print(f"\n✋ SKIPPING: '{os.path.basename(gene_summary)}' already exists.")
            continue

        design_path = os.path.join(design_base_dir, exp, f"design_matrix_{exp}.txt")
        count_table = os.path.join(count_dir, f"{exp}_.count.txt")

        if not os.path.exists(design_path) or not os.path.exists(count_table):
            print(f"⚠️ Skipping {exp}: Input files not found.")
            continue

        cmd = [
            "mageck", "mle",
            "-k", count_table,
            "-d", design_path,
            "-n", output_prefix,
            "--threads", str(n_cores),
            "--permutation-round", str(perm_round),
            "--max-sgrnapergene-permutation", str(max_sgrna_perm),
            "--control-sgrna", control_file,
            "--norm-method", "control"
        ]

        print(f"\n🚀 RUNNING: {exp} | {param_suffix}")
        
        try:
            subprocess.run(cmd, check=True)
            print(f"✅ Success. Files saved with prefix: {os.path.basename(output_prefix)}")
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed: {e}")

def run_mle_qc_plots(mle_dir, essential_txt, pos_control_csv):
    """
    Unified QC: Generates Distributions and Tracking with extra-large text 
    and a 4-column layout for better visibility.
    """
    # 1. Load Essential Lists
    eg_df = pd.read_table(essential_txt, sep='\t', header=None)
    eg_ids = set(eg_df.iloc[:, 0].astype(str).unique().tolist())
    
    pos_df = pd.read_csv(pos_control_csv)
    gold_std_ids = set(pos_df.iloc[:, 0].astype(str).unique().tolist())

    for root, dirs, files in os.walk(mle_dir):
        for file in files:
            if not file.endswith(".gene_summary.txt"): continue
            
            file_path = os.path.join(root, file)
            exp_tag = file.split('_')[0]
            df = pd.read_table(file_path)
            
            beta_cols = [c for c in df.columns if '|beta' in c]
            if not beta_cols: continue

            # --- PLOT 1: SYNCHRONIZED DISTRIBUTIONS (4 COLUMNS) ---
            df['Group'] = df['Gene'].apply(lambda x: 'Essential' if str(x) in eg_ids else 'Experimental')
            
            # Global Axis Scaling
            global_min = min(df[beta_cols].min().min(), -1.0)
            global_max = max(df[beta_cols].max().max(), 0.5)
            x_limit = (global_min * 1.1, global_max * 1.1)

            # Find global max density for Y-axis sync
            max_y = 0
            for col in beta_cols:
                data = df[col].dropna()
                if len(data) > 1:
                    temp_ax = sns.kdeplot(data=data)
                    line = temp_ax.lines[-1]
                    max_y = max(max_y, np.max(line.get_ydata()))
                    plt.close()
            y_limit = (0, max_y * 1.1)

            # 4-Column Layout
            cols_per_row = 4
            n_rows = math.ceil(len(beta_cols) / cols_per_row)
            fig, axes = plt.subplots(n_rows, cols_per_row, figsize=(24, 6 * n_rows), squeeze=False)
            fig.suptitle(f"CRISPR Screen Distributions: {exp_tag}", fontsize=30, y=1.05)

            for i, col in enumerate(beta_cols):
                ax = axes[divmod(i, cols_per_row)]
                sns.kdeplot(data=df[df['Group']=='Experimental'], x=col, ax=ax, fill=True, color='gray', alpha=0.3, label='Exp Pool')
                sns.kdeplot(data=df[df['Group']=='Essential'], x=col, ax=ax, color='red', lw=3, label='Essential')
                
                ax.set_xlim(x_limit)
                ax.set_ylim(y_limit)
                ax.axvline(0, color='black', linestyle='--', lw=2, alpha=0.8)
                
                # Big Text Adjustments
                ax.set_title(col.split('|')[0], fontsize=20, fontweight='bold')
                ax.set_xlabel("Beta Score", fontsize=18)
                ax.set_ylabel("Density", fontsize=18)
                ax.tick_params(axis='both', which='major', labelsize=14)
                ax.legend(loc='upper left', fontsize=12)

            for j in range(i + 1, n_rows * cols_per_row): fig.delaxes(axes.flatten()[j])
            plt.tight_layout()
            dist_name = f"QC_Distributions_{exp_tag}.png"
            plt.savefig(os.path.join(root, dist_name), dpi=300, bbox_inches='tight')
            plt.close()

            # --- PLOT 2: GOLD STANDARD TRACKING (Big Text) ---
            subset = df[df['Gene'].astype(str).isin(gold_std_ids)].copy()
            if not subset.empty:
                melted = subset.melt(id_vars='Gene', value_vars=beta_cols, var_name='Condition', value_name='Beta')
                melted['Condition'] = melted['Condition'].str.replace('|beta', '', regex=False)

                plt.figure(figsize=(16, 9))
                plt.axhline(0, color='black', linestyle='--', lw=3, alpha=0.9, label="Zero Baseline")
                sns.lineplot(data=melted, x='Condition', y='Beta', hue='Gene', marker='o', alpha=0.7, markersize=10, lw=2)
                
                plt.ylim(min(melted['Beta'].min() - 0.2, -1.2), 0.6) 
                plt.title(f"Essential Gene Beta Tracking: {exp_tag}", fontsize=26, fontweight='bold')
                plt.xticks(rotation=30, ha='right', fontsize=16)
                plt.yticks(fontsize=16)
                plt.ylabel("Beta Score (Fitness Effect)", fontsize=20)
                plt.xlabel("Condition", fontsize=20)
                plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=11, title="Phatr3 ID", title_fontsize=13, ncol=1)
                plt.grid(axis='y', linestyle=':', alpha=0.6)
                
                plt.tight_layout()
                track_name = f"QC_Essential_Tracking_{exp_tag}.png"
                plt.savefig(os.path.join(root, track_name), dpi=300, bbox_inches='tight')
                plt.close()
                
            print(f"✅ Finished {exp_tag}: Large-text plots saved.")

def quantile_normalize_data(df, beta_cols):
    """
    Standardizes distributions across conditions.
    
    How it works:
    1. It sorts the Beta scores for each condition (ranking).
    2. It calculates the mean value for each rank across all selected columns.
    3. It reassigns these mean values back to the genes based on their original rank.
    Result: All conditions now have the same mean and variance.
    """
    sub_df = df[beta_cols].copy()
    # Find the average value for each rank across all conditions
    rank_mean = sub_df.stack().groupby(sub_df.rank(method='first').stack().astype(int)).mean()
    # Map those average values back to the original ranks
    normalized = sub_df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return normalized

def calculate_fdr_quantile_matching(delta_betas):
    """
    Estimates FDR using the Quantile Matching protocol.
    
    How it works:
    1. Delta-Beta: The difference between treatment and control (T - C).
    2. Sigma (σ): We find the 68.2 percentile of absolute delta-betas. 
       In a normal distribution, this is exactly 1 standard deviation.
    3. Z-score: delta_beta / sigma.
    4. P-value: Derived from the Z-score using a standard normal distribution.
    5. FDR: Benjamini-Hochberg correction for multiple testing.
    """
    # 68th percentile of absolute values is a robust estimator for Sigma
    sigma = np.quantile(np.abs(delta_betas), 0.68)
    
    # Avoid division by zero if sigma is somehow 0
    if sigma == 0:
        return np.ones(len(delta_betas))
        
    z_scores = delta_betas / sigma
    # Two-tailed p-value
    p_values = 2 * (1 - stats.norm.cdf(np.abs(z_scores)))
    
    # Benjamini-Hochberg FDR correction
    _, fdr, _, _ = multipletests(p_values, method='fdr_bh')
    return fdr

def run_mle_quantile_norm(mle_dir, norm_dir, selector_csv):
    """
    Coordinates Normalization using high-level category matching.
    """
    all_gene_files = []
    for root, dirs, files in os.walk(mle_dir):
        for f in files:
            if f.endswith(".gene_summary.txt"):
                all_gene_files.append(os.path.join(root, f))
    
    if not all_gene_files:
        print(f"❌ No gene_summary.txt files found in {mle_dir}")
        return

    # 1. Regenerate CSV if Experiment column is wrong
    if not os.path.exists(selector_csv):
        print(f"📄 Generating clean comparison template: {selector_csv}")
        template_rows = []
        for file_path in all_gene_files:
            path_parts = os.path.normpath(file_path).split(os.sep)
            # path_parts[-3] is 'Metals' or 'Nitrogen'
            exp_cat = path_parts[-3] 
            
            df_temp = pd.read_table(file_path, nrows=1)
            beta_cols = [c for c in df_temp.columns if '|beta' in c]
            for c1 in beta_cols:
                for c2 in beta_cols:
                    if c1 != c2: template_rows.append([exp_cat, c1, c2, 0])
        
        pd.DataFrame(template_rows, columns=['Experiment', 'Treatment_Col', 'Control_Col', 'Analyze_Flag']).to_csv(selector_csv, index=False)
        print("🛑 Template ready. Set Analyze_Flag=1 and run again.")
        return

    # 2. Process pairs
    selector = pd.read_csv(selector_csv)
    active_pairs = selector[selector['Analyze_Flag'] == 1]
    
    for file_path in all_gene_files:
        path_parts = os.path.normpath(file_path).split(os.sep)
        exp_cat = path_parts[-3] 
        perm_folder = path_parts[-2]
        
        # Exact string matching
        exp_pairs = active_pairs[active_pairs['Experiment'] == exp_cat]
        
        if exp_pairs.empty:
            print(f"⏭️ Skipping {perm_folder}: Category '{exp_cat}' not active in CSV.")
            continue
        
        print(f"🧪 Processing {exp_cat} -> {perm_folder}")
        df = pd.read_table(file_path)
        involved_cols = list(set(exp_pairs['Treatment_Col'].tolist() + exp_pairs['Control_Col'].tolist()))
        
        # Math: Norm and FDR
        df[involved_cols] = quantile_normalize_data(df, involved_cols)
        
        for _, row in exp_pairs.iterrows():
            treat, ctrl = row['Treatment_Col'], row['Control_Col']
            label = f"{treat.split('|')[0]}_vs_{ctrl.split('|')[0]}"
            delta = df[treat] - df[ctrl]
            df[f"{label}_delta_beta"] = delta
            df[f"{label}_fdr"] = calculate_fdr_quantile_matching(delta.values)

        # 3. Mirrored Output
        exp_out_dir = os.path.join(norm_dir, exp_cat, perm_folder)
        os.makedirs(exp_out_dir, exist_ok=True)
        out_name = os.path.basename(file_path).replace(".gene_summary.txt", "_normalized_hits.csv")
        df.to_csv(os.path.join(exp_out_dir, out_name), index=False)
        print(f"✅ Saved to: {exp_out_dir}/{out_name}")

def plot_gene_of_interest_volcano(file_path, genes_dict, out_dir, hard_db, hard_fdr, len_db, len_fdr):
    """
    Highlights specific genes with distinct markers and labels.
    Legend is placed horizontally on top, and title is removed.
    """
    df = pd.read_csv(file_path)
    os.makedirs(out_dir, exist_ok=True)
    
    db_cols = [c for c in df.columns if '_delta_beta' in c]
    
    # Distinct markers and colors for your targets
    markers = ['o', 's', 'D', '^', 'v', 'p', '*', 'h']
    colors = sns.color_palette("husl", len(genes_dict))
    
    for db_col in db_cols:
        prefix = db_col.replace('_delta_beta', '')
        fdr_col = f"{prefix}_fdr"
        df['neg_log_fdr'] = -np.log10(df[fdr_col] + 1e-12)
        
        plt.figure(figsize=(20, 14))
        
        # 1. Background (Increased visibility: s=80, alpha=0.3)
        sns.scatterplot(data=df, x=db_col, y='neg_log_fdr', 
                        color='lightgray', alpha=0.3, s=80, label='Other', zorder=1)
        
        # 2. Draw Significance Lines (Increased weight for visibility)
        plt.axvline(hard_db, color='red', linestyle='--', alpha=0.8, lw=3)
        plt.axhline(-np.log10(hard_fdr), color='red', linestyle='--', alpha=0.8, lw=3)
        
        plt.axvline(len_db, color='orange', linestyle=':', alpha=0.8, lw=3)
        plt.axhline(-np.log10(len_fdr), color='orange', linestyle=':', alpha=0.8, lw=3)

        # 3. Highlight Specific Genes
        for i, (gene_id, common_name) in enumerate(genes_dict.items()):
            gene_row = df[df['Gene'] == gene_id]
            
            if not gene_row.empty:
                x = gene_row[db_col].values[0]
                y = gene_row['neg_log_fdr'].values[0]
                
                marker = markers[i % len(markers)]
                color = colors[i]
                
                plt.scatter(x, y, color=color, edgecolor='black', s=600, 
                            marker=marker, zorder=10, label=f"{common_name}", alpha=0.9)

        # 4. Styling: Remove Title, Dynamic X-label, Massive Fonts
        plt.xlabel(f"$\Delta$Beta ({prefix})", fontsize=36, fontweight='bold', labelpad=20)
        plt.ylabel("-log10(FDR)", fontsize=36, fontweight='bold', labelpad=20)
        
        plt.grid(linestyle=':', alpha=0.4)
        plt.tick_params(axis='both', which='major', labelsize=28, length=12, width=4)
       
        # 5. Horizontal Legend on top
        # bbox_to_anchor places it just above the plot; ncol spreads it horizontally
        plt.legend(bbox_to_anchor=(0.5, 1.02), loc='lower center', 
                   ncol=6, fontsize=24, frameon=False, title_fontsize=20, labelspacing=1.5)
        
        save_name = f"Volcano_Targets_{prefix}.png"
        plt.savefig(os.path.join(out_dir, save_name), dpi=300, bbox_inches='tight')
        plt.close()

def generate_hit_reports_and_plots(norm_dir, hard_db, hard_fdr, len_db, len_fdr, top_n, genes_dict):
    """
    Generates Tiered CSVs and High-Vis Volcano Plots.
    Summary and Legend are placed horizontally outside/on top of the graph.
    """
    subfolder_master_data = {}

    for root, dirs, files in os.walk(norm_dir):
        for file in files:
            if not file.endswith("_normalized_hits.csv"): continue
            
            file_path = os.path.join(root, file)
            target_profile_dir = os.path.join(root, "__genes_profile")
            plot_gene_of_interest_volcano(file_path, genes_dict, target_profile_dir, hard_db, hard_fdr, len_db, len_fdr)

            df = pd.read_csv(file_path)
            path_parts = os.path.normpath(file_path).split(os.sep)
            perm_tag = path_parts[-2] 
            
            db_cols = [c for c in df.columns if '_delta_beta' in c]
            
            for db_col in db_cols:
                prefix = db_col.replace('_delta_beta', '')
                fdr_col = f"{prefix}_fdr"
                
                def get_significance_label(row):
                    if row[db_col] <= hard_db and row[fdr_col] <= hard_fdr:
                        return "Tier 1 (Hard)"
                    elif row[db_col] <= len_db and row[fdr_col] <= len_fdr:
                        return "Tier 2 (Lenient)"
                    return "None"

                df['Significance_Tier'] = df.apply(get_significance_label, axis=1)
                
                n_hard = len(df[df['Significance_Tier'] == "Tier 1 (Hard)"])
                n_lenient = len(df[df['Significance_Tier'] == "Tier 2 (Lenient)"])
                
                # --- Master Data Storage Logic ---
                current_hits = df[df['Significance_Tier'] != "None"].copy()
                if not current_hits.empty:
                    master_entry = current_hits.copy()
                    master_entry['Comparison'] = prefix
                    master_entry = master_entry.rename(columns={'Gene': 'Gene_ID', db_col: 'Delta_Beta', fdr_col: 'FDR'})
                    if root not in subfolder_master_data:
                        subfolder_master_data[root] = []
                    subfolder_master_data[root].append(master_entry[['Gene_ID', 'Comparison', 'Delta_Beta', 'FDR', 'Significance_Tier']])

                # --- High-Vis Plotting ---
                fig, ax = plt.subplots(figsize=(20, 14))
                df['neg_log_fdr'] = -np.log10(df[fdr_col] + 1e-12)
                
                sns.scatterplot(data=df[df['Significance_Tier'] == "None"], 
                                x=db_col, y='neg_log_fdr', color='lightgray', alpha=0.3, s=80, label='Not Significant')
                sns.scatterplot(data=df[df['Significance_Tier'] == "Tier 2 (Lenient)"], 
                                x=db_col, y='neg_log_fdr', color='orange', alpha=0.7, s=250, label=f'Lenient (FDR < {len_fdr})')
                sns.scatterplot(data=df[df['Significance_Tier'] == "Tier 1 (Hard)"], 
                                x=db_col, y='neg_log_fdr', color='red', alpha=0.9, s=450, 
                                label=f'Hard (FDR < {hard_fdr})', edgecolor='black')

                plt.axvline(len_db, color='black', linestyle='--', alpha=0.5, lw=3)
                plt.axhline(-np.log10(len_fdr), color='black', linestyle='--', alpha=0.5, lw=3)

                # 1. HORIZONTAL SELECTION SUMMARY (Top Left)
                stats_text = f"[Hard:{n_hard}] | [Lenient:{n_lenient}]"
                ax.text(0.0, 1.05, stats_text, transform=ax.transAxes, fontsize=22,
                        verticalalignment='top', fontweight='bold', color='darkred')

                # 2. HORIZONTAL LEGEND (Top Right)
                # loc='lower right' relative to the bbox anchor makes it sit on top
                plt.legend(bbox_to_anchor=(1.0, 1.0), loc='lower right', 
                           ncol=3, fontsize=18, frameon=False, title_fontsize=20)

                plt.xlabel(f"$\Delta$Beta ({prefix}, {perm_tag})", fontsize=30, fontweight='bold', labelpad=20)
                plt.ylabel("-log10(FDR)", fontsize=30, fontweight='bold', labelpad=20)
                ax.tick_params(axis='both', which='major', labelsize=28, length=12, width=4)
                plt.grid(axis='both', linestyle=':', alpha=0.4, lw=2)
                
                plt.savefig(os.path.join(root, f"Volcano_{prefix}.png"), dpi=300, bbox_inches='tight')
                plt.close()

    # --- Save Master Summaries ---
    for folder_path, hits_list in subfolder_master_data.items():
        master_df = pd.concat(hits_list, ignore_index=True)
        master_df = master_df.sort_values(by=['Significance_Tier', 'Delta_Beta'])
        exp_name = os.path.basename(os.path.dirname(folder_path))
        master_path = os.path.join(folder_path, f"Master_Hits_{exp_name}_Summary.csv")
        master_df.to_csv(master_path, index=False)
        print(f"✅ Created Subfolder Master List: {master_path}")

def prepare_koala_database(koala_paths, output_path):
    """
    Step 11: Consolidates KEGG annotations into a global reference folder.
    """
    if os.path.exists(output_path):
        print(f"⏩ Step 11: Global KEGG file already exists. Skipping.")
        return

    # Create the 'processed' folder if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    print("🚀 Step 11: Merging KEGG/KOALA files into Global Reference...")
    all_dfs = []

    for fpath in koala_paths:
        if not os.path.exists(fpath):
            print(f"⚠️ Warning: Missing: {os.path.basename(fpath)}")
            continue
            
        # Load col 0: GeneID, col 1: K-number
        df = pd.read_csv(fpath, sep='\t', header=None, names=['Gene_ID', 'KEGG_ID'], 
                         on_bad_lines='skip', low_memory=False)
        
        # Clean suffixes (.p1, .t1)
        df['Gene_ID'] = df['Gene_ID'].astype(str).str.split('.').str[0]
        df = df.dropna(subset=['KEGG_ID'])
        all_dfs.append(df)

    if not all_dfs:
        print("❌ Error: No data found.")
        return

    # Combine and deduplicate
    master_koala = pd.concat(all_dfs, ignore_index=True).drop_duplicates()
    master_koala.to_csv(output_path, sep='\t', index=False)
    
    print(f"✨ Step 11 Complete: Saved to {output_path}")
    print(f"📊 Total unique entries: {len(master_koala)}")

def clean_eggnog_file(input_path, output_path):
    """
    Cleans the raw eggNOG file, keeping ALL columns for maximum data retention.
    """
    if os.path.exists(output_path):
        print(f"⏩ Cleaned eggNOG already exists. Skipping.")
        return

    print(f"🧹 Cleaning eggNOG file (All Columns): {os.path.basename(input_path)}")
    
    try:
        # 1. Find the header row dynamically
        header_line = 0
        with open(input_path, 'r') as f:
            for i, line in enumerate(f):
                if line.startswith('#query'):
                    header_line = i
                    break
        
        # 2. Read the full file
        df = pd.read_csv(
            input_path, 
            sep='\t', 
            skiprows=header_line,
            engine='python'
        )
        
        # 3. Rename '#query' to 'Gene_ID'
        df.rename(columns={df.columns[0]: 'Gene_ID'}, inplace=True)
        
        # 4. Clean Gene IDs: 'Phatr3_draftJ1000.t1' -> 'Phatr3_draftJ1000'
        df['Gene_ID'] = df['Gene_ID'].astype(str).str.split('.').str[0]
        
        # 5. Handle potential empty rows/trailers
        df = df.dropna(subset=['Gene_ID'])
        
        # 6. Save as a comprehensive TSV
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        df.to_csv(output_path, sep='\t', index=False)
        
        print(f"✨ Success! Full database saved to: {output_path}")
        print(f"📊 Columns preserved: {len(df.columns)}")

    except Exception as e:
        print(f"❌ Error processing eggNOG file: {e}")

def normalize_id(gene_id):
    """
    Standardizes Phatr3 IDs by stripping 'draft' and underscores.
    Example: 
    'Phatr3_Jdraft278' -> 'phatr3j278'
    'Phatr3_draftJ913' -> 'phatr3j913'
    'Phatr3_EG01223'   -> 'phatr3eg01223'
    """
    if pd.isna(gene_id):
        return ""
    
    # 1. Lowercase and strip whitespace
    name = str(gene_id).lower().strip()
    
    # 2. Remove the word 'draft' entirely
    name = name.replace('draft', '')
    
    # 3. Remove underscores to handle 'phatr3_j' vs 'phatr3j'
    name = name.replace('_', '')
    
    # 4. Strip protein versioning (everything after the first dot)
    name = name.split('.')[0]
    
    return name

def annotate_all_experiment_hits(norm_dir, eggnog_path):
    """
    Step 12: Annotates hits while preserving the Normalized ID for verification.
    """
    if not os.path.exists(eggnog_path):
        print(f"❌ Error: eggNOG database not found at {eggnog_path}")
        return

    print(f"📖 Loading Global eggNOG Database...")
    eggnog_df = pd.read_table(eggnog_path)
    
    # Create the normalized key for the database side
    eggnog_df['Normalized_ID'] = eggnog_df['Gene_ID'].apply(normalize_id)
    # Remove original Gene_ID from the eggnog side to avoid name collisions (Gene_ID_x/y)
    eggnog_subset = eggnog_df.drop(columns=['Gene_ID'])

    for root, dirs, files in os.walk(norm_dir):
        # Find every Master Hits summary file
        target_files = [f for f in files if f.startswith("Master_Hits_") and f.endswith("_Summary.csv")]
        
        for file in target_files:
            master_path = os.path.join(root, file)
            out_path = master_path.replace(".csv", "_ANNOTATED.csv")
            
            print(f"🧬 Processing and Annotating: {file}")
            hits_df = pd.read_csv(master_path)
            
            # Create the normalized key on the CRISPR hits side
            hits_df['Normalized_ID'] = hits_df['Gene_ID'].apply(normalize_id)
            
            # Perform the merge. 
            # This keeps your original 'Gene_ID' and adds 'Normalized_ID' + eggNOG columns.
            annotated_df = pd.merge(
                hits_df, 
                eggnog_subset, 
                on='Normalized_ID', 
                how='left'
            )
            
            # Fill missing annotations with N/A for clean viewing
            annotated_df = annotated_df.fillna("N/A")
            
            # Reorder columns slightly so IDs are at the very beginning
            cols = annotated_df.columns.tolist()
            # Move Normalized_ID to be the second column right after the original Gene_ID
            if 'Normalized_ID' in cols:
                cols.insert(1, cols.pop(cols.index('Normalized_ID')))
                annotated_df = annotated_df[cols]

            annotated_df.to_csv(out_path, index=False)
            print(f"   ✅ Saved Annotated file with ID verification to: {os.path.basename(out_path)}")

    print("\n✨ Step 12 Complete: You can now compare 'Gene_ID' and 'Normalized_ID' in your CSVs.")

def run_functional_plotting(norm_dir, plot_dir, cog_mapping):
    """
    Step 13: Generates combined COG plots. 
    Fixes IndexError by ensuring all Categories have both Tier 1 and Tier 2 rows.
    """
    print("\n🚀 Step 13: Generating Combined Tier COG Plots...")
    
    for root, _, files in os.walk(norm_dir):
        for file in [f for f in files if f.endswith("_ANNOTATED.csv")]:
            csv_path = os.path.join(root, file)
            df = pd.read_csv(csv_path)
            
            exp_name = file.replace("_ANNOTATED.csv", "")
            exp_root_dir = os.path.join(plot_dir, exp_name)
            os.makedirs(exp_root_dir, exist_ok=True)

            for comp in df['Comparison'].unique():
                comp_df = df[df['Comparison'] == comp].copy()
                target_tiers = ['Tier 1 (Hard)', 'Tier 2 (Lenient)']
                plot_df = comp_df[comp_df['Significance_Tier'].isin(target_tiers)].copy()
                
                if plot_df.empty:
                    continue

                # 1. Map COG and Wrap
                plot_df['COG_Full'] = plot_df['COG_category'].astype(str).str[0].str.upper().map(cog_mapping)
                plot_df['COG_Full'] = plot_df['COG_Full'].fillna('Unclassified / No COG')
                
                def manual_wrap(text, words_per_line=3):
                    words = str(text).split()
                    return "\n".join([" ".join(words[i:i+words_per_line]) for i in range(0, len(words), words_per_line)])
                
                plot_df['Category'] = plot_df['COG_Full'].apply(manual_wrap)

                # 2. Aggregate counts
                counts = plot_df.groupby(['Category', 'Significance_Tier']).size().reset_index(name='Count')

                # --- THE FIX: Reindex to ensure every Category has both Tiers ---
                all_cats = counts['Category'].unique()
                mux = pd.MultiIndex.from_product([all_cats, target_tiers], names=['Category', 'Significance_Tier'])
                counts = counts.set_index(['Category', 'Significance_Tier']).reindex(mux, fill_value=0).reset_index()
                # ---------------------------------------------------------------

                # Sort order based on Tier 1
                sort_order = counts[counts['Significance_Tier'] == 'Tier 1 (Hard)'].sort_values('Count', ascending=False)['Category']
                
                # Dynamic Spacing
                num_categories = len(all_cats)
                plt.figure(figsize=(20, max(10, num_categories * 1.2)))

                ax = sns.barplot(
                    data=counts, 
                    x='Count', 
                    y='Category', 
                    hue='Significance_Tier', 
                    order=sort_order,
                    palette={'Tier 1 (Hard)': '#2c7fb8', 'Tier 2 (Lenient)': '#7fcdbb'}
                )
                
                plt.title(f"Functional Enrichment: {comp}\nExperiment: {exp_name}", fontsize=28, pad=40, fontweight='bold')
                plt.xlabel("Number of Genes", fontsize=24, labelpad=20)
                plt.ylabel("", labelpad=30)
                plt.xticks(fontsize=20)
                plt.yticks(fontsize=18)
                plt.legend(title="Significance Tier", title_fontsize='22', fontsize='20', loc='lower right')
                
                for container in ax.containers:
                    ax.bar_label(container, padding=10, fontsize=16, fontweight='bold')
                
                plt.savefig(os.path.join(exp_root_dir, f"{comp}_Combined_COG.png"), bbox_inches='tight', dpi=300)
                plt.close()
            
            print(f"   ✅ Successfully generated plots for {exp_name}")

def plot_top_hits_by_beta(norm_dir, plot_dir, cog_mapping):
    """
    Step 13: Top 20 essential gene plots. 
    Maintains consistent colors for COG categories across all experiments
    and uses higher intensity colors for better slide visibility.
    """
    print("\n🏆 Step 13: Generating Top 20 Essential Gene Plots (Consistent Colors)...")
    
    # 1. CREATE GLOBAL COLOR MAP
    # This ensures "Category A" is always "Color X" across every plot you generate.
    all_cog_names = sorted(list(set(cog_mapping.values())) + ['Unclassified / No COG'])
    # Using 'Set2' or 'Dark2' for higher intensity than pastels
    colors = sns.color_palette("Set2", len(all_cog_names))
    cog_color_dict = dict(zip(all_cog_names, colors))

    for root, _, files in os.walk(norm_dir):
        for file in [f for f in files if f.endswith("_ANNOTATED.csv")]:
            csv_path = os.path.join(root, file)
            df = pd.read_csv(csv_path)
            
            exp_name = file.replace("_ANNOTATED.csv", "")
            exp_root_dir = os.path.join(plot_dir, exp_name)
            os.makedirs(exp_root_dir, exist_ok=True)

            for comp in df['Comparison'].unique():
                comp_df = df[df['Comparison'] == comp].copy()
                
                target_tiers = ['Tier 1 (Hard)', 'Tier 2 (Lenient)']
                hits_df = comp_df[comp_df['Significance_Tier'].isin(target_tiers)].copy()
                
                if hits_df.empty:
                    continue

                # 2. Sort and Take Top 20
                top_20 = hits_df.sort_values('Delta_Beta', ascending=True).head(20).copy()

                # 3. Cleanup Data
                top_20['COG_Full'] = top_20['COG_category'].astype(str).str[0].str.upper().map(cog_mapping)
                top_20['COG_Full'] = top_20['COG_Full'].fillna('Unclassified / No COG')
                top_20['Description'] = top_20['Description'].fillna('-')
                
                id_to_desc = dict(zip(top_20['Gene_ID'], top_20['Description']))

                # 4. Setup Plot
                plt.figure(figsize=(22, 14))
                
                # Use the fixed cog_color_dict to ensure consistency across files
                ax = sns.barplot(
                    data=top_20,
                    x='Delta_Beta',
                    y='Gene_ID',
                    hue='COG_Full',
                    dodge=False, 
                    palette=cog_color_dict 
                )

                # 5. Styling (Dynamic X-axis, No Title)
                plt.xlabel(f"Delta Beta ({comp})", fontsize=24, labelpad=20, fontweight='bold')
                plt.ylabel("Gene ID", fontsize=22, labelpad=10)
                plt.xticks(fontsize=18)
                plt.yticks(fontsize=18)
                
                # 6. Legend: 3 columns, positioned at top
                plt.legend(title="COG Functional Category", title_fontsize='22', fontsize='16', 
                           loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=3, frameon=False)

                # 7. Text placement: Description inside the bar
                ytick_labels = ax.get_yticklabels()
                for i, label_obj in enumerate(ytick_labels):
                    gene_id = label_obj.get_text()
                    desc = str(id_to_desc.get(gene_id, "-"))
                    display_text = (desc[:75] + '...') if len(desc) > 75 else desc
                    
                    ax.text(-0.02, i, f"{display_text}  ", 
                            color='black', va='center', ha='right', 
                            fontsize=15, fontweight='bold')

                # 8. Save
                plot_filename = f"{comp}_Top20_Essential_Genes.png"
                plt.savefig(os.path.join(exp_root_dir, plot_filename), bbox_inches='tight', dpi=300)
                plt.close()

            print(f"   ✅ Consistent Top 20 plots saved for: {exp_name}")

