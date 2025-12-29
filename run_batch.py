
import sys
import os
from xenium_pipeline_main import XeniumPipeline, load_config, setup_logging

def run_batch():
    datasets = [
        ("human_brain", "./data/unprocessed_adata/human_brain.h5ad"),
        ("h_breast_1", "./data/unprocessed_adata/h_breast_1.h5ad")
    ]
    
    # Steps: [0, 1, 2, 4, 6] (권장 전체 분석)
    steps = [0, 1, 2, 4, 6]

    for tag, path in datasets:
        print(f"\n{'='*60}")
        print(f"Processing Dataset: {tag}")
        print(f"{'='*60}")
        
        try:
            pipeline = XeniumPipeline(
                xenium_input_path=path,
                output_path="./output",
                sample_tag=tag,
                load_preprocessed=True
            )
            
            pipeline.run_full_pipeline(steps=steps)
            print(f"✓ {tag} processing complete.")
            
        except Exception as e:
            print(f"✗ Failed to process {tag}: {e}")

if __name__ == "__main__":
    # 로깅 설정 초기화 (한 번만)
    config = load_config("config.yaml")
    setup_logging(config)
    
    run_batch()
