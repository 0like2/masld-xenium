#!/bin/bash
# run_all.sh
# Runs Xenium pipeline on 3 datasets in order of size (smallest first).

echo "Starting Batch Execution..."
date

# 1. h_breast_2 (~650MB)
echo "----------------------------------------"
echo "1/3. Skipping h_breast_2 (Completed)..."
echo "----------------------------------------"
# python run_batch.py --config config_h_breast_2.yaml

# 2. human_brain (~1.2GB)
echo "----------------------------------------"
echo "2/3. Running human_brain..."
echo "----------------------------------------"
python run_batch.py --config config_human_brain.yaml

# 3. h_breast_1 (~5GB)
echo "----------------------------------------"
echo "3/3. Running h_breast_1..."
echo "----------------------------------------"
python run_batch.py --config config_h_breast_1.yaml

echo "All datasets completed."
date
