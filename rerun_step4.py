
from xenium_pipeline_main import XeniumPipeline, setup_logging, load_config
import os

config = load_config("config_human_brain.yaml")
setup_logging(config)

# Initialize pipeline with existing Step 1 H5AD
pipeline = XeniumPipeline(
    xenium_input_path="./output/human_brain/human_brain_step1_preprocessed.h5ad", 
    sample_tag="human_brain",
    load_preprocessed=True
)

# Run Step 4 (This should be fast)
print("Running Step 4 to generate optimization curve data...")
pipeline.step4_optimal_expansion()
print("Done!")
