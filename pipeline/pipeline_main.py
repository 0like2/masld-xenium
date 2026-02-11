# pipeline/pipeline_main.py
# ------------------------------------------------------------
# MASLD Xenium Pipeline Orchestrator
# ------------------------------------------------------------

import os
import sys
import argparse
import yaml
import logging
from pathlib import Path

# Add current directory to path to allow imports if run from root
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

# Import Steps
import xenium_step0_formatting as step0
import xenium_step1_dataset_exploration as step1
import xenium_step2_segmentation_free_analysis as step2
import xenium_step3_resegmentation as step3
import xenium_step4_techniques_comparison as step4
import xenium_step5_optimal_expansion as step5
import xenium_step6_segmentation_benchmark as step6
import xenium_step7_simulation as step7

# Configure Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_config(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def get_config_for_step(base_config, step_output_dir):
    """Returns a copy of config with updated output_dir"""
    new_config = base_config.copy()
    new_config["output_dir"] = step_output_dir
    return new_config

def main():
    print("="*60)
    print("MASLD Xenium Pipeline Started")
    print("="*60)


    # --- 1. Load Config ---
    parser = argparse.ArgumentParser(description="MASLD Xenium Pipeline")
    parser.add_argument("--config", default=os.path.join(current_dir, "config.yaml"),
                        help="Path to config YAML file")
    args = parser.parse_args()
    config_path = args.config
    if not os.path.exists(config_path):
        logger.error(f"Config file not found at {config_path}")
        return

    config = load_config(config_path)
    print(f"Loaded configuration from {config_path}")
    
    # --- 1.1 Determine Output Structure ---
    # Input Path: .../data/SampleName_outs
    # Output Structure: xenium-output/SampleName_outs/stepX_name
    input_path = config["input_path"]
    
    # Attempt to extract a meaningful sample name from the input path
    # If input path ends with '/', strip it first
    input_path_clean = input_path.rstrip(os.sep)
    sample_name = os.path.basename(input_path_clean)
    
    # Base output directory
    base_output_root = config.get("output_dir", "xenium-output") # Default if not set
    sample_output_dir = os.path.join(base_output_root, sample_name)
    
    print(f"Sample Name Derived: {sample_name}")
    print(f"Sample Output Directory: {sample_output_dir}")
    os.makedirs(sample_output_dir, exist_ok=True)
    
    # Need to keep strict track of where data is for subsequent steps
    # We will update these as we go
    step0_adata_path = None
    step1_adata_path = None
    transcripts_csv_path = None # Will find in Step 0 or Input
    step3_adata_path = None
    step3_transcripts_path = None
    
    sample_tag = config["sample_tag"]


    # ---  2. Run Step 0: Formatting ---
    print("\n" + "-"*40)
    print("Step 0: Formatting")
    print("-"*40)
    
    step0_dir = os.path.join(sample_output_dir, "step0_formatting")
    os.makedirs(step0_dir, exist_ok=True)
    
    step0_config = get_config_for_step(config, step0_dir)
    step0_output_file = os.path.join(step0_dir, f"{sample_tag}.h5ad")
    
    if os.path.exists(step0_output_file):
         print(f"Skipping Step 0: Output found at {step0_output_file}")
    else:
        try:
            adata = step0.run_step0(step0_config)
            print("Step 0 completed successfully.")
        except Exception as e:
            logger.error(f"Step 0 failed: {e}")
            return # Critical failure
            
    step0_adata_path = step0_output_file
    
    # Search for transcripts.csv (or .parquet if configured)
    if os.path.exists(os.path.join(input_path, "transcripts.csv")):
        transcripts_csv_path = os.path.join(input_path, "transcripts.csv")
    elif os.path.exists(os.path.join(step0_dir, "transcripts.csv")):
        transcripts_csv_path = os.path.join(step0_dir, "transcripts.csv")

    if not transcripts_csv_path and config.get('use_parquet', False):
        parquet_path = os.path.join(input_path, "transcripts.parquet")
        if os.path.exists(parquet_path):
            transcripts_csv_path = parquet_path
            print(f"  > Using parquet transcripts: {parquet_path}")

    # --- 3. Run Step 1: Dataset Exploration ---
    print("\n" + "-"*40)
    print("Step 1: Dataset Exploration")
    print("-"*40)
    
    step1_dir = os.path.join(sample_output_dir, "step1_exploration")
    os.makedirs(step1_dir, exist_ok=True)
    step1_config = get_config_for_step(config, step1_dir)
    step1_config['previous_step_adata_path'] = step0_adata_path
    
    step1_output_file = os.path.join(step1_dir, f"{sample_tag}_step1_exploration.h5ad")
    
    if os.path.exists(step1_output_file):
         print(f"Skipping Step 1: Output found at {step1_output_file}")
    else:
        try:
            step1.run_step1(step1_config)
            print("Step 1 (Exploration) completed successfully.")
        except Exception as e:
            logger.error(f"Step 1 failed: {e}")
            return

    step1_adata_path = step1_output_file

    # --- 4. Run Step 2: Segmentation Free Analysis (Points2Regions) ---
    print("\n" + "-"*40)
    print("Step 2: Segmentation Free Analysis (Points2Regions)")
    print("-"*40)
    
    step2_dir = os.path.join(sample_output_dir, "step2_segmentation_free")
    os.makedirs(step2_dir, exist_ok=True)
    step2_config = get_config_for_step(config, step2_dir)
    step2_config['previous_step_adata_path'] = step1_adata_path
    
    step2_output_file = os.path.join(step2_dir, f"{sample_tag}_step2_points2regions.h5ad")
    
    if os.path.exists(step2_output_file):
         print(f"Skipping Step 2: Output found at {step2_output_file}")
    else:
        try:
            step2.run_step2(step2_config)
            print("Step 2 (Segmentation Free Analysis) completed successfully.")
        except Exception as e:
            logger.error(f"Step 2 failed: {e}")
            return


    # --- 5. Run Step 3: Resegmentation (Cellpose) ---
    # Formerly Step 4
    print("\n" + "-"*40)
    print("Step 3: Resegmentation (Cellpose)")
    print("-"*40)
    
    step3_dir = os.path.join(sample_output_dir, "step3_resegmentation")
    os.makedirs(step3_dir, exist_ok=True)
    step3_config = get_config_for_step(config, step3_dir)
    
    step3_marker = os.path.join(step3_dir, "step3_done.txt")
    
    step3_output_adata = os.path.join(step3_dir, f"{sample_tag}_step3_resegmented.h5ad")
    step3_output_transcripts = os.path.join(step3_dir, f"{sample_tag}_step3_transcripts_resegmented.csv")
    step3_output_masks = os.path.join(step3_dir, f"{sample_tag}_step3_resegmented_masks.tif")
    
    step3_mask_path = None

    if os.path.exists(step3_marker):
        print("Step 3 (Resegmentation) output already exists. Skipping...")
        step3_adata_path = step3_output_adata
        step3_transcripts_path = step3_output_transcripts
        if os.path.exists(step3_output_masks):
            step3_mask_path = step3_output_masks
    else:
        print("Running Step 3 (Resegmentation)...")
        
        # Find DAPI
        dapi_path = os.path.join(input_path, "morphology_focus.ome.tif")
        if not os.path.exists(dapi_path):
            dapi_path = os.path.join(input_path, "morphology_focus", "morphology_focus_0000.ome.tif")
        if not os.path.exists(dapi_path):
             dapi_path = os.path.join(input_path, "DAPI.tif")
        
        # Transcripts path (original)
        curr_transcripts = transcripts_csv_path if transcripts_csv_path else os.path.join(input_path, "transcripts.csv")

        # Check for generated domain map from Step 2
        potential_domain_map = os.path.join(step2_dir, "points2regions", "domain_polygons.json")
        if os.path.exists(potential_domain_map):
            print(f"  > Using automatically generated domain map: {potential_domain_map}")
            step3_config.setdefault('resegmentation', {}).setdefault('domain_assignment', {})['domain_map_path'] = potential_domain_map

        step3.run_step3(step3_config, dapi_path, curr_transcripts, step3_dir)

        if os.path.exists(step3_output_adata):
            with open(step3_marker, "w") as f:
                f.write("done")
            step3_adata_path = step3_output_adata
            step3_transcripts_path = step3_output_transcripts
            if os.path.exists(step3_output_masks):
                step3_mask_path = step3_output_masks
        else:
            print("  > Step 3 did not produce output (disabled or failed).")

    # --- 6. Run Step 4: Comparison & Validation (Metrics) ---
    # Formerly Step 3
    print("\n" + "-"*40)
    print("Step 4: Comparison & Validation")
    print("-"*40)
    
    step4_dir = os.path.join(sample_output_dir, "step4_techniques_comparison")
    os.makedirs(step4_dir, exist_ok=True)
    step4_config = get_config_for_step(config, step4_dir)
    
    step4_marker = os.path.join(step4_dir, "step4_done.txt")
    
    if os.path.exists(step4_marker):
        print("Step 4 (Comparison) output already exists. Skipping...")
    else:
        print("Running Step 4 (Comparison)...")
        # Needs Resegmented Data from Step 3
        if step3_adata_path and os.path.exists(step3_adata_path):
            print(f"  > Using resegmented data from Step 3: {step3_adata_path}")
            current_adata_path = step3_adata_path
            current_transcripts_path = step3_transcripts_path
        else:
            print("  > Step 3 output not found. Falling back to Step 1 data (Comparison only, no Validation).")
            current_adata_path = step1_adata_path
            current_transcripts_path = transcripts_csv_path # Raw transcripts
            
        step4.run_step4(step4_config, current_adata_path, current_transcripts_path, step4_dir, original_adata_path=step1_adata_path, original_transcripts_path=transcripts_csv_path)
        with open(step4_marker, "w") as f:
            f.write("done")


    # --- 7. Run Step 5: Optimal Expansion ---
    print("\n" + "-"*40)
    print("Step 5: Optimal Expansion")
    print("-"*40)
    
    step5_dir = os.path.join(sample_output_dir, "step5_optimal_expansion")
    os.makedirs(step5_dir, exist_ok=True)
    step5_config = get_config_for_step(config, step5_dir)
    
    # Prefer Step 3 (Resegmented) data for expansion if available
    if step3_adata_path and os.path.exists(step3_adata_path):
        step5_config['previous_step_adata_path'] = step3_adata_path
    else:
        step5_config['previous_step_adata_path'] = step1_adata_path
        
    step5_marker = os.path.join(step5_dir, f"{sample_tag}_step5_done.txt")
    if os.path.exists(step5_marker):
         print(f"Skipping Step 5: Output marker found at {step5_marker}")
    else:
        try:
            step5.run_step5(step5_config)
            step5_csv = os.path.join(step5_dir, f"{sample_tag}_step5_expanded_transcripts.csv")
            if os.path.exists(step5_csv):
                with open(step5_marker, "w") as f:
                    f.write("done")
                print("Step 5 (Optimal Expansion) completed successfully.")
            else:
                print("  > Step 5 did not produce output (disabled or failed).")
        except Exception as e:
            logger.error(f"Step 5 failed: {e}")

            
    # --- 8. Run Step 6: Segmentation Benchmark (with Optional Baysor) ---
    print("\n" + "-"*40)
    print("Step 6: Segmentation Benchmark (incl. Baysor)")
    print("-"*40)
    
    step6_dir = os.path.join(sample_output_dir, "step6_benchmark")
    os.makedirs(step6_dir, exist_ok=True)
    step6_config = get_config_for_step(config, step6_dir)
    
    # Inject input_dir for Baysor raw transcript access
    step6_config['input_dir'] = input_path
    
    # Inject inputs for benchmark
    # 1. Nuclei (Step 0)
    if step0_adata_path:
        step6_config.setdefault('benchmark', {})['input_nuclei'] = step0_adata_path
        
    # 2. Advanced (Cellpose from Step 3)
    if step3_adata_path:
        step6_config.setdefault('benchmark', {})['input_advanced'] = step3_adata_path
        
    # 3. Prior for Baysor (TIF Mask from Step 3)
    if step3_mask_path:
        step6_config.setdefault('baysor', {})['prior_segmentation_tif'] = step3_mask_path
        
    # 4. Optimal Expansion (Step 5)
    # Step 5 outputs a CSV: {sample_tag}_step5_expanded_transcripts.csv
    # Note: Currently Step 5 logic assigns DOMAINS not CELL IDs.
    step5_csv = os.path.join(step5_dir, f"{sample_tag}_step5_expanded_transcripts.csv")
    if os.path.exists(step5_csv):
        step6_config.setdefault('benchmark', {})['input_expansion'] = step5_csv
        
    step6_marker = os.path.join(step6_dir, "step6_done.txt")
    
    if os.path.exists(step6_marker):
        print("Step 6 output already exists. Skipping...")
    else:
        try:
            step6.run_step6(step6_config)
            with open(step6_marker, "w") as f:
                f.write("done")
            print("Step 6 completed successfully.")
        except Exception as e:
            logger.error(f"Step 6 failed: {e}")

    # --- 9. Run Step 7: Simulation ---
    print("\n" + "-"*40)
    print("Step 7: Simulation & Benchmarking")
    print("-"*40)
    
    step7_dir = os.path.join(sample_output_dir, "step7_simulation")
    os.makedirs(step7_dir, exist_ok=True)
    step7_config = get_config_for_step(config, step7_dir)
    
    step7_marker = os.path.join(step7_dir, f"{sample_tag}_step7_done.txt")
    if os.path.exists(step7_marker):
        print(f"Skipping Step 7: Output marker found at {step7_marker}")
    else:
        try:
            step7.run_step7(step7_config)
            with open(step7_marker, "w") as f:
                f.write("done")
            print("Step 7 (Simulation) completed successfully.")
        except Exception as e:
            logger.error(f"Step 7 failed: {e}")


    print("\n" + "="*60)
    print("Pipeline Execution Finished")
    print("="*60)

if __name__ == "__main__":
    main()
