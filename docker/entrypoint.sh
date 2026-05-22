#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# MASLD Xenium Pipeline - Docker Entrypoint
# ============================================================
# Overrides config.yaml paths via environment variables,
# then runs the pipeline or a diagnostic command.
#
# Environment variables:
#   INPUT_PATH   - Path to Xenium input data (mounted volume)
#   OUTPUT_DIR   - Path to write results (mounted volume)
#   SAMPLE_TAG   - Sample identifier (e.g. human_alzheimers)
#   RUN_MODE     - "run" (default), "test", or "shell"
# ============================================================

RUN_MODE="${RUN_MODE:-run}"

echo "========================================"
echo "MASLD Xenium Pipeline - Docker Container"
echo "========================================"
echo "Run mode  : ${RUN_MODE}"
echo "Input path: ${INPUT_PATH:-<not set>}"
echo "Output dir: ${OUTPUT_DIR:-<not set>}"
echo "Sample tag: ${SAMPLE_TAG:-<not set>}"
echo ""

# GPU / CUDA status
if command -v nvidia-smi &>/dev/null; then
    echo "--- GPU Info ---"
    nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader 2>/dev/null || echo "nvidia-smi failed"
    echo ""
fi

python -c "
import torch
print(f'PyTorch {torch.__version__}')
print(f'  CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'  CUDA device: {torch.cuda.get_device_name(0)}')
mps = hasattr(torch.backends, 'mps') and torch.backends.mps.is_available()
print(f'  MPS available:  {mps}')
" 2>/dev/null || true
echo ""

# --- Mode: shell (debugging) ---
if [ "${RUN_MODE}" = "shell" ]; then
    echo "Dropping into shell..."
    exec /bin/bash
fi

# --- Mode: test (import & GPU check) ---
if [ "${RUN_MODE}" = "test" ]; then
    echo "--- Running import tests ---"
    python -c "
import sys
modules = [
    'numpy', 'pandas', 'scipy', 'anndata', 'scanpy', 'squidpy',
    'torch', 'cellpose', 'shapely', 'geopandas', 'rasterio',
    'tifffile', 'h5py', 'yaml', 'xb',
]
failed = []
for m in modules:
    try:
        __import__(m)
        print(f'  [OK] {m}')
    except ImportError as e:
        print(f'  [FAIL] {m}: {e}')
        failed.append(m)

import torch
print(f'\ntorch.cuda.is_available() = {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'CUDA GPU: {torch.cuda.get_device_name(0)}')
mps = hasattr(torch.backends, 'mps') and torch.backends.mps.is_available()
print(f'torch.backends.mps.is_available() = {mps}')

if failed:
    print(f'\nFailed imports: {failed}')
    sys.exit(1)
print('\nAll imports OK.')
"
    exit $?
fi

# --- Mode: run (pipeline) ---
if [ "${RUN_MODE}" = "run" ]; then
    # Generate runtime config by patching paths
    RUNTIME_CONFIG="/tmp/runtime_config.yaml"
    ORIGINAL_CONFIG="/app/pipeline/config.yaml"

    python -c "
import yaml, sys, os

config_path = '${ORIGINAL_CONFIG}'
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)

# Override paths from environment variables
input_path = os.environ.get('INPUT_PATH')
output_dir = os.environ.get('OUTPUT_DIR')
sample_tag = os.environ.get('SAMPLE_TAG')

if input_path:
    config['input_path'] = input_path
if output_dir:
    config['output_dir'] = output_dir
if sample_tag:
    config['sample_tag'] = sample_tag

out_path = '${RUNTIME_CONFIG}'
with open(out_path, 'w') as f:
    yaml.dump(config, f, default_flow_style=False)

print(f'Runtime config written to {out_path}')
print(f'  input_path : {config.get(\"input_path\")}')
print(f'  output_dir : {config.get(\"output_dir\")}')
print(f'  sample_tag : {config.get(\"sample_tag\")}')
"

    echo ""
    echo "--- Starting pipeline ---"
    exec python /app/pipeline/pipeline_main.py --config "${RUNTIME_CONFIG}"
fi

echo "Unknown RUN_MODE: ${RUN_MODE}"
echo "Valid modes: run, test, shell"
exit 1
