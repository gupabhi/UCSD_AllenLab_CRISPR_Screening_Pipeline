import os

# =============================================================================
# CRISPR MAGeCK Pipeline — Configuration
# Round 2: 2026_02_CRISPR_Pool2  |  AVITI SE150  |  TSCC
# =============================================================================

# --- 1. EXPERIMENT IDENTIFIER ---
EXP_NAME = "2026_02_CRISPR_Pool2"

# --- 2. PIPELINE STEP CONTROL ---
# Set each step True/False to run or skip.
# Typical first run: Steps 0-5, then 6 (pause to fill design matrix), then 7-13.
STEP_0_CHECK_ENV          = True
STEP_1_INITIALIZE         = True
STEP_1b_PREFLIGHT_CHECK   = True   # Sanity-check anchor presence in raw reads (fast)
STEP_2_METADATA_TEMPLATE  = True  # Pre-filled CSV is provided; keep False
STEP_3_VALIDATE_AND_SUMMARY = True
STEP_4_EXTRACT_GUIDES     = True
STEP_5_COUNT              = False
STEP_6_DESIGN_TEMPLATE    = False
STEP_7_MLE                = False
STEP_8_MLE_QC             = False
STEP_9_QUANTILE_NORM      = False
STEP_10_VISUALIZE_HITS    = False
STEP_11_PREPARE_DB        = False
STEP_12_ANNOTATE_HITS     = False
STEP_13_PLOT_COG          = False

# --- 3. SEQUENCING MODE ---
# SE150 from AVITI — reads come from one direction only; disable reverse-complement search.
SINGLE_END = True

# Set True when running under SLURM (sbatch) to skip interactive y/n prompts.
# Set False when running interactively to allow manual inspection at Step 3.
NON_INTERACTIVE = True

# --- 4. COMPUTE ---
N_CORES = 32   # Must match -c in the .sh SBATCH header

# --- 5. DIRECTORY STRUCTURE ---
EXP_INPUT_DIR  = f"input/{EXP_NAME}"
EXP_OUTPUT_DIR = f"output/{EXP_NAME}"

# Input sub-directories (relative to pipeline root)
LIB_DIR = os.path.join(EXP_INPUT_DIR, "common/library")
ANN_DIR = os.path.join(EXP_INPUT_DIR, "common/PT_protein_anno")
GLOBAL_DB_DIR = os.path.join(ANN_DIR, "processed")

# Raw sequencing data — absolute TSCC path, no file copying needed.
# Files are named:  {SampleID}_R1.fastq.gz  (SE150, R1 only)
RAW_DATA_DIR = (
    "/tscc/projects/ps-allenlab/archdata/zfussy/"
    "allen_lab_ampliconseq/CRISPR-pool2/AAHZ_CRIPR_POOL2"
)

# --- 6. INPUT FILES ---
META_CSV              = os.path.join(EXP_INPUT_DIR, f"metadata_{EXP_NAME}.csv")
LIBRARY_CSV           = os.path.join(LIB_DIR, "JonScreenAllGuides_mageck_format.csv")
NEG_CONTROL_TXT       = os.path.join(LIB_DIR, "NegativeControl_mageck_format.txt")
POS_CONTROL_CSV       = os.path.join(LIB_DIR, "PositiveControl_20EssentialGenes.csv")
ESSENTIAL_GENES_TXT   = os.path.join(LIB_DIR, "essential_genes_annots.txt")
COMPARISON_SELECTOR_CSV = os.path.join(EXP_INPUT_DIR, "mle_comparison_pairs.csv")
RAW_EGGNOG            = os.path.join(ANN_DIR, "Phatr3.emapper.annotations")
GO_ASSOC_FILE         = os.path.join(ANN_DIR, "Phatr_goatools-association.tsv")
KOALA_FILES = [
    os.path.join(ANN_DIR, "2026_03_06_GhostKOALA_Ptricornutum_annot.txt"),
    os.path.join(ANN_DIR, "2026_03_06_BlastKOALA_Ptricornutum_annot_1.txt"),
    os.path.join(ANN_DIR, "2026_03_06_BlastKOALA_Ptricornutum_annot_2.txt"),
]

# --- 7. OUTPUT DIRECTORIES ---
PREPROCESS_DIR    = os.path.join(EXP_OUTPUT_DIR, "1_preprocess_reads")
COUNT_DIR         = os.path.join(EXP_OUTPUT_DIR, "2_mageck_count")
DESIGN_BASE_DIR   = os.path.join(EXP_OUTPUT_DIR, "3_mageck_design_matrix")
MLE_DIR           = os.path.join(EXP_OUTPUT_DIR, "4_mageck_mle")
NORM_DIR          = os.path.join(EXP_OUTPUT_DIR, "5_mageck_mle_norm")
COG_PLOT_DIR      = os.path.join(EXP_OUTPUT_DIR, "6_Plots_top_hits", "COG_distribution")
BETA_PLOT_DIR     = os.path.join(EXP_OUTPUT_DIR, "6_Plots_top_hits", "BETA_TopHits_Description")
COG_SPECIFIC_PLOT_DIR = os.path.join(EXP_OUTPUT_DIR, "6_Plots_top_hits", "COG_Specific_Hits")

# --- 8. OUTPUT FILES ---
PREPROCESS_SUMMARY_CSV = os.path.join(PREPROCESS_DIR, "extraction_summary.csv")
UNIFIED_KOALA_MASTER   = os.path.join(GLOBAL_DB_DIR, "Phatr3_Unified_KEGG_KOALA.tsv")
CLEAN_EGGNOG           = os.path.join(GLOBAL_DB_DIR, "Phatr3_Cleaned_eggNOG.tsv")
CONTROL_SGRNA_FILE     = NEG_CONTROL_TXT

# --- 9. PREPROCESSING CONSTANTS ---
FWD_ANCHOR = "CAAAAAACACCTTCAAAGTC"   # U6 promoter sequence (upstream of guide)
REV_ANCHOR = "GCTATTTCTAGCTCTAAAAC"   # Scaffold sequence (only used if SINGLE_END=False)
GUIDE_LEN  = 20

# Preflight check: sample this many reads from the first FASTQ to verify anchor presence.
PREFLIGHT_N_READS    = 10_000
PREFLIGHT_MIN_RATE   = 0.10   # Warn (but don't abort) if < 10 % of sampled reads hit the anchor

# --- 10. MLE PARAMETERS ---
PERMUTATION_ROUND        = 10
MAX_SGRNA_PERMUTATION    = 15

# --- 11. HIT CRITERIA ---
HARD_DB_CUTOFF   = -1.5   # Delta-Beta: strong depletion
HARD_FDR_CUTOFF  = 0.01   # 1 % FDR
LENIENT_DB_CUTOFF  = -1.0
LENIENT_FDR_CUTOFF = 0.05
TOP_HIT_COUNT = 50

# --- 12. GENES OF INTEREST ---
GENES_DICT = {
    "Phatr3_J54983":  "NR",
    "Phatr3_J42577":  "bZIP18",
    "Phatr3_J42514":  "HSF1g",
    "Phatr3_J43051":  "HSF1a",
    "Phatr3_J50624":  "J50624",
    "Phatr3_EG01412": "EG01412",
    "Phatr3_J11337":  "J11337",
    "Phatr3_J52260":  "J52260",
}

COG_MAP = {
    "J": "[J] Translation, ribosomal structure and biogenesis",
    "A": "[A] RNA processing and modification",
    "K": "[K] Transcription",
    "L": "[L] Replication, recombination and repair",
    "B": "[B] Chromatin structure and dynamics",
    "D": "[D] Cell cycle control, cell division, chromosome partitioning",
    "Y": "[Y] Nuclear structure",
    "V": "[V] Defense mechanisms",
    "T": "[T] Signal transduction mechanisms",
    "M": "[M] Cell wall/membrane/envelope biogenesis",
    "N": "[N] Cell motility",
    "Z": "[Z] Cytoskeleton",
    "W": "[W] Extracellular structures",
    "U": "[U] Intracellular trafficking, secretion, and vesicular transport",
    "O": "[O] Posttranslational modification, protein turnover, chaperones",
    "C": "[C] Energy production and conversion",
    "G": "[G] Carbohydrate transport and metabolism",
    "E": "[E] Amino acid transport and metabolism",
    "F": "[F] Nucleotide transport and metabolism",
    "H": "[H] Coenzyme transport and metabolism",
    "I": "[I] Lipid transport and metabolism",
    "P": "[P] Inorganic ion transport and metabolism",
    "Q": "[Q] Secondary metabolites biosynthesis, transport and catabolism",
    "R": "[R] General function prediction only",
    "S": "[S] Function unknown",
}
COG_interests = ["A", "K"]

# --- 13. OPTIONAL: PCR CYCLE COMPARISON ---
# When True, the pipeline runs a separate MLE for PCR_Cycle_Compare experiment,
# comparing 12-cycle vs 24-cycle library prep for matched samples.
# Relevant 12c files: T0_NHNOA12C1/2, T1_NO3B2_12C1/2, T1_NFRB_12C1/2
ENABLE_PCR_CYCLE_COMPARISON = False
