
import visualize_results as vis
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)

# Config
xenium_input = Path("/data1/project/20rak/masld_xenium/human_brain") # Raw input path
output_dir = Path("./output/human_brain")
sample_tag = "human_brain"

# Check if raw path exists
if not xenium_input.exists():
    # If using formatted data, we might not have raw input here?
    # Let's check where the user keeps raw data. 
    # From previous logs: xenium_input_path="./output/human_brain/human_brain_step1_preprocessed.h5ad" was used for pipeline init.
    # BUT Step 0 format_xenium takes raw input.
    # Where is the raw input?
    # User's other corpus: /data1/project/20rak/masld_xenium -> /data1/project/20rak/masld_xenium
    # Let's search for transcripts.parquet in the project dir to find the raw data.
    pass

# For now, let's try to assume a standard path or pass the formatted folder if it has copies?
# The formatted folder ./output/human_brain only has h5ad and csvs.
# The `xb/formatting.py` *copies* transcripts.csv if decompressed.
# Wait, XeniumPipeline init usually takes the RAW folder. 
# I will search for the raw folder.

vis.visualize_segmentation_overlay(xenium_input, output_dir, sample_tag)
