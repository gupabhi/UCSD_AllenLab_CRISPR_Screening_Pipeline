import os

# --- 1. EXPERIMENT IDENTIFIER ---
EXP_NAME = "2025_08_CRISPR_Pool1"

# --- 2. PIPELINE CONTROL ---
STEP_0_CHECK_ENV = False
STEP_1_INITIALIZE = False
STEP_2_METADATA_TEMPLATE = False
STEP_3_VALIDATE_AND_SUMMARY = False 
STEP_4_EXTRACT_GUIDES = False
STEP_5_COUNT = False
STEP_6_DESIGN_TEMPLATE = False  
STEP_7_MLE = False
STEP_8_MLE_QC = False
STEP_9_QUANTILE_NORM = False
STEP_10_VISUALIZE_HITS = False
STEP_11_PREPARE_DB = False
STEP_12_ANNOTATE_HITS = False
STEP_13_PLOT_COG = True

# Number of CPU cores assigned in your .sh file
N_CORES = 32

# --- 3. DIRECTORY STRUCTURE ---
EXP_INPUT_DIR = f"input/{EXP_NAME}"
EXP_OUTPUT_DIR = f"output/{EXP_NAME}"

# Input Data Folders
LIB_DIR = os.path.join(EXP_INPUT_DIR, "1_sgRNAs")
SEQ_PARENT_DIR = os.path.join(EXP_INPUT_DIR, "2_sequencing")
ANN_DIR = os.path.join(EXP_INPUT_DIR, "3_annotation")
GLOBAL_DB_DIR = os.path.join(ANN_DIR, "processed")

# Input Data Files
SEQ_SUBFOLDER = "251223_LH00444_0455_A22YNNMLT3"
RAW_DATA_DIR = os.path.join(SEQ_PARENT_DIR, SEQ_SUBFOLDER)

META_CSV = os.path.join(EXP_INPUT_DIR, f"metadata_{EXP_NAME}.csv")
LIBRARY_CSV = os.path.join(LIB_DIR, "JonScreenAllGuides_mageck_format.csv")
NEG_CONTROL_TXT = os.path.join(LIB_DIR, "NegativeControl_mageck_format.txt")
POS_CONTROL_CSV = os.path.join(LIB_DIR,  "PositiveControl_20EssentialGenes.csv")
ESSENTIAL_GENES_TXT = os.path.join(LIB_DIR, "essential_genes_annots.txt")
COMPARISON_SELECTOR_CSV = os.path.join(EXP_INPUT_DIR, "mle_comparison_pairs.csv")
RAW_EGGNOG = os.path.join(ANN_DIR, "Phatr3.emapper.annotations")
GO_ASSOC_FILE = os.path.join(ANN_DIR, "Phatr_goatools-association.tsv")
KOALA_FILES = [
    os.path.join(ANN_DIR, "2026_03_06_GhostKOALA_Ptricornutum_annot.txt"),
    os.path.join(ANN_DIR, "2026_03_06_BlastKOALA_Ptricornutum_annot_1.txt"),
    os.path.join(ANN_DIR, "2026_03_06_BlastKOALA_Ptricornutum_annot_2.txt")
]

# Output Data Folders
PREPROCESS_DIR = os.path.join(EXP_OUTPUT_DIR, "1_preprocess_reads")
COUNT_DIR = os.path.join(EXP_OUTPUT_DIR, "2_mageck_count")
DESIGN_BASE_DIR = os.path.join(EXP_OUTPUT_DIR, "3_mageck_design_matrix")
MLE_DIR = os.path.join(EXP_OUTPUT_DIR, "4_mageck_mle")
NORM_DIR = os.path.join(EXP_OUTPUT_DIR, "5_mageck_mle_norm")
COG_PLOT_DIR = os.path.join(EXP_OUTPUT_DIR, "6_Plots_top_hits", "COG_distribution")
BETA_PLOT_DIR = os.path.join(EXP_OUTPUT_DIR, "6_Plots_top_hits", "BETA_TopHits_Description")
COG_SPECIFIC_PLOT_DIR = os.path.join(EXP_OUTPUT_DIR, "6_Plots_top_hits", "COG_Specific_Hits")

# Output Data Files
PREPROCESS_SUMMARY_CSV =  os.path.join(PREPROCESS_DIR, "extraction_summary.csv")
UNIFIED_KOALA_MASTER = os.path.join(GLOBAL_DB_DIR, "Phatr3_Unified_KEGG_KOALA.tsv")        
CLEAN_EGGNOG = os.path.join(GLOBAL_DB_DIR, "Phatr3_Cleaned_eggNOG.tsv")

# --- 4. PREPROCESSING CONSTANTS ---
FWD_ANCHOR = "CAAAAAACACCTTCAAAGTC"  #U6 Promoter Sequence
REV_ANCHOR = "GCTATTTCTAGCTCTAAAAC" # scaffold sequence
GUIDE_LEN = 20

# --- 5. MLE PARAMETERS ---
PERMUTATION_ROUND = 10 
MAX_SGRNA_PERMUTATION = 15
CONTROL_SGRNA_FILE = NEG_CONTROL_TXT # Using your existing NEG_CONTROL_TXT variable

# MLE QC Analysis
CONTROL_SGRNA_FILE = NEG_CONTROL_TXT 
ESSENTIAL_GENES_FILE = POS_CONTROL_CSV

# --- 7. HIT CRITERIA ---
HARD_DB_CUTOFF = -1.5     # Delta-Beta: Strong depletion
HARD_FDR_CUTOFF = 0.01    # High confidence (1% FDR)

LENIENT_DB_CUTOFF = -1.0  # Moderate depletion
LENIENT_FDR_CUTOFF = 0.05 # Standard confidence (5% FDR)

TOP_HIT_COUNT = 50        # Number of genes to export per pair

# Your specific genes of interest mapping ID -> Common Name
GENES_DICT = {
    'Phatr3_J54983': 'NR',
    'Phatr3_J42577': 'bZIP18',
    'Phatr3_J42514': 'HSF1g',
    'Phatr3_J43051': 'HSF1a',
    'Phatr3_J50624': 'J50624',
    'Phatr3_EG01412': 'EG01412',
    'Phatr3_J11337': 'J11337',
    'Phatr3_J52260': 'J52260'
}

COG_MAP = {
    # INFORMATION STORAGE AND PROCESSING
    'J': '[J] Translation, ribosomal structure and biogenesis',
    'A': '[A] RNA processing and modification',
    'K': '[K] Transcription',
    'L': '[L] Replication, recombination and repair',
    'B': '[B] Chromatin structure and dynamics',

    # CELLULAR PROCESSES AND SIGNALING
    'D': '[D] Cell cycle control, cell division, chromosome partitioning',
    'Y': '[Y] Nuclear structure',
    'V': '[V] Defense mechanisms',
    'T': '[T] Signal transduction mechanisms',
    'M': '[M] Cell wall/membrane/envelope biogenesis',
    'N': '[N] Cell motility',
    'Z': '[Z] Cytoskeleton',
    'W': '[W] Extracellular structures',
    'U': '[U] Intracellular trafficking, secretion, and vesicular transport',
    'O': '[O] Posttranslational modification, protein turnover, chaperones',

    # METABOLISM
    'C': '[C] Energy production and conversion',
    'G': '[G] Carbohydrate transport and metabolism',
    'E': '[E] Amino acid transport and metabolism',
    'F': '[F] Nucleotide transport and metabolism',
    'H': '[H] Coenzyme transport and metabolism',
    'I': '[I] Lipid transport and metabolism',
    'P': '[P] Inorganic ion transport and metabolism',
    'Q': '[Q] Secondary metabolites biosynthesis, transport and catabolism',

    # POORLY CHARACTERIZED
    'R': '[R] General function prediction only',
    'S': '[S] Function unknown'
}

COG_interests  = ['A', 'K']