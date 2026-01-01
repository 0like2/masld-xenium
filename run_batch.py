import sys
import os
import argparse
from xenium_pipeline_main import XeniumPipeline, load_config, setup_logging

def run_from_config(config_path):
    print(f"Loading configuration from: {config_path}")
    config = load_config(config_path)
    
    # Extract parameters from config
    input_path = config.get('input_anndata_path')
    output_dir = config.get('output_dir', './output')
    sample_tag = config.get('sample_tag', 'sample')
    steps = config.get('pipeline_steps', [0, 1, 2, 4, 6])
    
    # Tissue type for GT
    tissue_type = config.get('tissue_type', 'brain') # Default to brain if not specified
    
    print(f"\n{'='*60}")
    print(f"Processing Dataset: {sample_tag}")
    print(f"Input: {input_path}")
    print(f"Steps: {steps}")
    print(f"{'='*60}")

    # Prepare specific kwargs for pipeline steps
    pipeline_kwargs = {
        'generate_ground_truth': True, # Always try to generate if missing
        'tissue_type': tissue_type,
        'download_dir': 'data/reference',
        'simulation_kwargs': config.get('simulation_kwargs', {}),
        'preprocess_kwargs': config.get('preprocess_kwargs', {})
    }

    try:
        pipeline = XeniumPipeline(
            xenium_input_path=input_path,
            output_path=output_dir,
            sample_tag=sample_tag,
            load_preprocessed=True # Default behavior
        )
        
        pipeline.run_full_pipeline(steps=steps, **pipeline_kwargs)
        print(f"✓ {sample_tag} processing complete.")
        
    except Exception as e:
        print(f"✗ Failed to process {sample_tag}: {e}")
        raise e

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Xenium Pipeline from Config")
    parser.add_argument("--config", type=str, required=True, help="Path to config.yaml")
    args = parser.parse_args()
    
    # Setup logging
    BASE_CONFIG = load_config("config.yaml") # Load base config for logging settings
    setup_logging(BASE_CONFIG)
    
    run_from_config(args.config)
