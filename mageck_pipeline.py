"""
mageck_pipeline.py  —  AllenLab CRISPR MAGeCK Pipeline

Run all steps by setting the STEP_* flags in mageck_config.py.
Designed for both interactive use and non-interactive SLURM submission.
"""

import os
import sys
import pandas as pd
import mageck_utils as utils
import mageck_config as config


def main():
    print("\n" + "=" * 70)
    print(f"  AllenLab CRISPR MAGeCK Pipeline")
    print(f"  Experiment : {config.EXP_NAME}")
    print(f"  Mode       : {'Single-end (SE)' if config.SINGLE_END else 'Paired-end (PE)'}")
    print(f"  Interactive: {not config.NON_INTERACTIVE}")
    print("=" * 70 + "\n")

    # -----------------------------------------------------------------------
    # Step 0: Check environment
    # -----------------------------------------------------------------------
    if config.STEP_0_CHECK_ENV:
        print("--- [Step 0] Checking environment ---")
        missing = utils.check_env()
        if missing:
            print(f"❌ Missing requirements: {missing}")
            print("   Install them, then re-run.")
            sys.exit(1)
        print("✅ All requirements satisfied.\n")

    # -----------------------------------------------------------------------
    # Step 1: Initialise directory skeleton
    # -----------------------------------------------------------------------
    if config.STEP_1_INITIALIZE:
        print("--- [Step 1] Initialising workspace ---")
        utils.initialize_workspace(
            config.EXP_INPUT_DIR,
            config.EXP_OUTPUT_DIR,
            config.LIB_DIR,
            config.ANN_DIR,
        )
        print()

    # -----------------------------------------------------------------------
    # Step 1b: Preflight anchor check (fast sanity check on raw reads)
    # -----------------------------------------------------------------------
    if config.STEP_1b_PREFLIGHT_CHECK:
        print("--- [Step 1b] Preflight anchor check ---")
        if not os.path.isdir(config.RAW_DATA_DIR):
            print(f"❌ RAW_DATA_DIR not found: {config.RAW_DATA_DIR}")
            sys.exit(1)
        utils.preflight_anchor_check(
            raw_dir    = config.RAW_DATA_DIR,
            fwd_anchor = config.FWD_ANCHOR,
            n_reads    = config.PREFLIGHT_N_READS,
            min_hit_rate = config.PREFLIGHT_MIN_RATE,
            single_end = config.SINGLE_END,
            rev_anchor = config.REV_ANCHOR,
        )

    # -----------------------------------------------------------------------
    # Step 2: Generate metadata template (skip if pre-built CSV provided)
    # -----------------------------------------------------------------------
    if config.STEP_2_METADATA_TEMPLATE:
        print("--- [Step 2] Generating metadata template ---")
        utils.generate_metadata_template(config.META_CSV, config.RAW_DATA_DIR, auto_parse=config.AUTO_PARSE_METADATA)
        print()

    # -----------------------------------------------------------------------
    # Step 3: Validate metadata and show experimental summary
    # -----------------------------------------------------------------------
    if config.STEP_3_VALIDATE_AND_SUMMARY:
        print("--- [Step 3] Validating metadata ---")
        is_valid, msg = utils.validate_metadata(
            config.META_CSV,
            non_interactive=config.NON_INTERACTIVE,
            single_end=config.SINGLE_END,
        )
        if not is_valid:
            print(f"❌ Metadata error: {msg}")
            sys.exit(1)
        print(f"✅ {msg}\n")

    # -----------------------------------------------------------------------
    # Step 4: Guide extraction (preprocessing)
    # -----------------------------------------------------------------------
    if config.STEP_4_EXTRACT_GUIDES:
        print("--- [Step 4] Extracting guides ---")
        meta = pd.read_csv(config.META_CSV)
        utils.run_guide_extraction(
            meta_df    = meta,
            raw_dir    = config.RAW_DATA_DIR,
            out_dir    = config.PREPROCESS_DIR,
            fwd        = config.FWD_ANCHOR,
            rev        = config.REV_ANCHOR,
            length     = config.GUIDE_LEN,
            n_cores    = config.N_CORES,
            single_end = config.SINGLE_END,
        )
        print()

    # -----------------------------------------------------------------------
    # Step 5: MAGeCK count
    # -----------------------------------------------------------------------
    if config.STEP_5_COUNT:
        print("--- [Step 5] MAGeCK count ---")
        utils.run_mageck_count(
            preprocess_dir = config.PREPROCESS_DIR,
            count_dir      = config.COUNT_DIR,
            library_csv    = config.LIBRARY_CSV,
            summary_csv    = config.PREPROCESS_SUMMARY_CSV,
            neg_path       = config.NEG_CONTROL_TXT,
            pos_path       = config.POS_CONTROL_CSV,
        )
        print()

    # -----------------------------------------------------------------------
    # Step 6: Design matrix management
    # -----------------------------------------------------------------------
    if config.STEP_6_DESIGN_TEMPLATE:
        print("--- [Step 6] Design matrix ---")
        ready_for_mle = utils.manage_all_design_matrices(
            config.PREPROCESS_SUMMARY_CSV,
            config.DESIGN_BASE_DIR,
        )
        if not ready_for_mle:
            print(
                "\n🛑 Pipeline paused at Step 6.\n"
                "   Open the design_matrix_*.txt file(s) in:\n"
                f"   {config.DESIGN_BASE_DIR}\n"
                "   Fill each row with 0 or 1:\n"
                "     - baseline column: always 1\n"
                "     - condition column: 1 for samples in that condition, 0 otherwise\n"
                "     - t0 rows: 0 in all condition columns (they ARE the baseline)\n"
                "   Then re-run with STEP_6_DESIGN_TEMPLATE=True to validate."
            )
            sys.exit(0)
        print()

    # -----------------------------------------------------------------------
    # Step 7: MAGeCK MLE
    # -----------------------------------------------------------------------
    if config.STEP_7_MLE:
        print("--- [Step 7] MAGeCK MLE ---")
        utils.run_all_mageck_mle(
            count_dir      = config.COUNT_DIR,
            mle_dir        = config.MLE_DIR,
            design_base_dir = config.DESIGN_BASE_DIR,
            n_cores        = config.N_CORES,
            perm_round     = config.PERMUTATION_ROUND,
            max_sgrna_perm = config.MAX_SGRNA_PERMUTATION,
            control_file   = config.CONTROL_SGRNA_FILE,
        )
        print()

    # -----------------------------------------------------------------------
    # Step 8: MLE QC plots + sanity tests
    # -----------------------------------------------------------------------
    if config.STEP_8_MLE_QC:
        print(f"--- [Step 8] MLE QC plots ---")
        utils.run_mle_qc_plots(
            mle_dir         = config.MLE_DIR,
            essential_txt   = config.ESSENTIAL_GENES_TXT,
            pos_control_csv = config.POS_CONTROL_CSV,
            count_dir       = config.COUNT_DIR,
        )
        print()

    # -----------------------------------------------------------------------
    # Step 9: Quantile normalisation & comparison
    # -----------------------------------------------------------------------
    if config.STEP_9_QUANTILE_NORM:
        print("--- [Step 9] Quantile normalisation ---")
        utils.run_mle_quantile_norm(
            mle_dir      = config.MLE_DIR,
            norm_dir     = config.NORM_DIR,
            selector_csv = config.COMPARISON_SELECTOR_CSV,
        )
        print()

    # -----------------------------------------------------------------------
    # Step 10: Hit reports & volcano plots
    # -----------------------------------------------------------------------
    if config.STEP_10_VISUALIZE_HITS:
        print("--- [Step 10] Hit reports & volcano plots ---")
        utils.generate_hit_reports_and_plots(
            norm_dir = config.NORM_DIR,
            hard_db  = config.HARD_DB_CUTOFF,
            hard_fdr = config.HARD_FDR_CUTOFF,
            len_db   = config.LENIENT_DB_CUTOFF,
            len_fdr  = config.LENIENT_FDR_CUTOFF,
            genes_dict = config.GENES_DICT,
        )
        utils.plot_condition_scatter(norm_dir = config.NORM_DIR)
        print()

    # -----------------------------------------------------------------------
    # Step 11: Annotation databases
    # -----------------------------------------------------------------------
    if config.STEP_11_PREPARE_DB:
        print("--- [Step 11] Preparing annotation databases ---")
        utils.prepare_koala_database(config.KOALA_FILES, config.UNIFIED_KOALA_MASTER)
        utils.clean_eggnog_file(config.RAW_EGGNOG, config.CLEAN_EGGNOG)
        print()

    # -----------------------------------------------------------------------
    # Step 12: Annotate hits
    # -----------------------------------------------------------------------
    if config.STEP_12_ANNOTATE_HITS:
        print("--- [Step 12] Annotating hits ---")
        utils.annotate_all_experiment_hits(config.NORM_DIR, config.CLEAN_EGGNOG)
        print()

    # -----------------------------------------------------------------------
    # Step 13: Functional plots
    # -----------------------------------------------------------------------
    if config.STEP_13_PLOT_COG:
        print("--- [Step 13] Functional / COG plots ---")
        utils.run_functional_plotting(
            norm_dir    = config.NORM_DIR,
            plot_dir    = config.COG_PLOT_DIR,
            cog_mapping = config.COG_MAP,
        )
        utils.plot_top_hits_by_beta(
            norm_dir    = config.NORM_DIR,
            plot_dir    = config.BETA_PLOT_DIR,
            cog_mapping = config.COG_MAP,
            top_n       = config.TOP_HIT_COUNT,
        )
        print()

    print("=" * 70)
    print("  ✅ PIPELINE COMPLETE")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
