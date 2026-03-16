import os
import sys
import pandas as pd
import mageck_utils as utils
import mageck_config as config

def main():
    # Step 0: Check Environment
    if config.STEP_0_CHECK_ENV:
        missing = utils.check_env()
        if missing:
            print(f"❌ Requirements missing: {missing}")
            sys.exit(1)

    # Step 1: Initialize Folders
    if config.STEP_1_INITIALIZE:
        is_new = utils.initialize_workspace(
            config.EXP_INPUT_DIR, config.EXP_OUTPUT_DIR, 
            config.SEQ_PARENT_DIR, config.SEQ_SUBFOLDER, config.LIB_DIR
        )
        if is_new:
            print("\n🛠️ Folders created. Please upload your files to the subfolder.")
            return

    # Step 2: Generate Metadata Template
    if config.STEP_2_METADATA_TEMPLATE:
        # This now looks specifically in the subfolder defined in config
        utils.generate_metadata_template(config.META_CSV, config.RAW_DATA_DIR)

    # Step 3: Validating Metadata and Generating Summary
    if config.STEP_3_VALIDATE_AND_SUMMARY:
        is_valid, msg = utils.validate_metadata(config.META_CSV)
        if not is_valid:
            print(f"❌ Metadata Error: {msg}")
            sys.exit(1)
        
        # Adding a manual pause so you can inspect the summary table
        user_input = input("Does the experimental summary above look correct? (y/n): ")
        if user_input.lower() != 'y':
            print("Operation cancelled by user. Please correct the metadata CSV.")
            sys.exit(0)
            
        print("✅ Metadata confirmed. Ready for next step.")

    # Step 4: Extracting Guides (Preprocessing)
    if config.STEP_4_EXTRACT_GUIDES:
        meta = pd.read_csv(config.META_CSV)
        
        utils.run_guide_extraction(
            meta, config.RAW_DATA_DIR, config.PREPROCESS_DIR, 
            config.FWD_ANCHOR, config.REV_ANCHOR, config.GUIDE_LEN,
            config.N_CORES # Passed from config
        )

    # Step 5: MAGeCK Count
    if config.STEP_5_COUNT:
  
        utils.run_mageck_count(
            config.PREPROCESS_DIR,
            config.COUNT_DIR,
            config.LIBRARY_CSV,
            config.PREPROCESS_SUMMARY_CSV,
            config.NEG_CONTROL_TXT,
            config.POS_CONTROL_CSV,
            config.GENES_DICT
        )
    
    # Step 6: Design Matrix Management
    if config.STEP_6_DESIGN_TEMPLATE:
        # Validates all experiments (Nitrogen, Metals, etc.) in one go
        ready_for_mle = utils.manage_all_design_matrices(
            config.PREPROCESS_SUMMARY_CSV, 
            config.DESIGN_BASE_DIR
        )
        
        if not ready_for_mle:
            print("\n🛑 Pipeline paused. Please fill the binary design matrices for each experiment.")
            sys.exit(0)
    
    # Step 7: Execute MLE
    if config.STEP_7_MLE:
        utils.run_all_mageck_mle(
            count_dir=config.COUNT_DIR,
            mle_dir=config.MLE_DIR,
            design_base_dir=config.DESIGN_BASE_DIR,
            n_cores=config.N_CORES,
            perm_round=config.PERMUTATION_ROUND,
            max_sgrna_perm=config.MAX_SGRNA_PERMUTATION,
            control_file=config.CONTROL_SGRNA_FILE
        )

    # Step 8: MLE Quality Control Plots
    if config.STEP_8_MLE_QC:
        print(f"\n--- Running Step 8 QC for {config.EXP_NAME} ---")
        utils.run_mle_qc_plots(
            mle_dir=config.MLE_DIR,
            essential_txt=config.ESSENTIAL_GENES_TXT,
            pos_control_csv=config.POS_CONTROL_CSV
        )
        
    if config.STEP_9_QUANTILE_NORM:
        print("\n--- [RUNNING] Step 9: Quantile Normalization & Comparative Analysis ---")
        utils.run_mle_quantile_norm(
            mle_dir=config.MLE_DIR,
            norm_dir=config.NORM_DIR,
            selector_csv=config.COMPARISON_SELECTOR_CSV
        )

    if config.STEP_10_VISUALIZE_HITS:
        print("\n--- [RUNNING] Step 10: Hit Reports & Volcano Plots ---")
        utils.generate_hit_reports_and_plots(
            norm_dir=config.NORM_DIR,
            hard_db=config.HARD_DB_CUTOFF,
            hard_fdr=config.HARD_FDR_CUTOFF,
            len_db=config.LENIENT_DB_CUTOFF,
            len_fdr=config.LENIENT_FDR_CUTOFF,
            genes_dict = config.GENES_DICT
        )

    if config.STEP_11_PREPARE_DB:
        print("\n--- [RUNNING] Step 11a: Preparing KOALA Annotation Database ---")
        utils.prepare_koala_database(
            koala_paths=config.KOALA_FILES,
            output_path=config.UNIFIED_KOALA_MASTER
        )
        print("\n--- [RUNNING] Step 11b: Preparing EGGNOG Annotation Database ---")
        utils.clean_eggnog_file(config.RAW_EGGNOG, config.CLEAN_EGGNOG)

    if config.STEP_12_ANNOTATE_HITS:
        print("\n--- [RUNNING] Step 12: Annotating CRISPR Hits ---")
        utils.annotate_all_experiment_hits(config.NORM_DIR, config.CLEAN_EGGNOG)

    if config.STEP_13_PLOT_COG:
        utils.run_functional_plotting(
            norm_dir=config.NORM_DIR,
            plot_dir=config.COG_PLOT_DIR,
            cog_mapping=config.COG_MAP
        )

        utils.plot_top_hits_by_beta(
            norm_dir=config.NORM_DIR, 
            plot_dir=config.BETA_PLOT_DIR, 
            cog_mapping=config.COG_MAP
        )

        utils.plot_cog_specific_hits(
            norm_dir=config.NORM_DIR, 
            plot_dir=config.COG_SPECIFIC_PLOT_DIR, 
            cog_mapping=config.COG_MAP,
            target_cogs=config.COG_interests,
            top_n=20  # Show top 20 genes instead of 10
        )


if __name__ == "__main__":
    main()