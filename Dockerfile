# ============================================================
# MASLD Xenium Pipeline - Multi-stage Docker Build
# ============================================================
# Base: nvidia/cuda 12.4 + Python 3.11
# Stage 1 (builder): compile wheels with build dependencies
# Stage 2 (runtime): lean image with runtime libs only
# ============================================================

# ----------------------------------------------------------
# Stage 1: Builder
# ----------------------------------------------------------
FROM nvidia/cuda:12.4.1-devel-ubuntu22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=1

# System build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 python3.11-dev python3.11-venv python3-pip \
    build-essential cmake git pkg-config \
    libhdf5-dev libgdal-dev libproj-dev libspatialindex-dev libgeos-dev \
    libffi-dev libssl-dev libjpeg-dev libpng-dev libtiff-dev \
    && rm -rf /var/lib/apt/lists/*

# Make python3.11 the default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 \
    && update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1

# Upgrade pip
RUN python -m pip install --no-cache-dir --upgrade pip setuptools wheel

# Build wheels from requirements
WORKDIR /build
COPY requirements.txt .
RUN pip wheel --no-cache-dir --wheel-dir=/build/wheels \
    --extra-index-url https://download.pytorch.org/whl/cu124 \
    -r requirements.txt

# ----------------------------------------------------------
# Stage 2: Runtime
# ----------------------------------------------------------
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04 AS runtime

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV MPLBACKEND=Agg

# Runtime system libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 python3.11-venv python3-pip \
    libhdf5-103 libgdal30 libproj22 libspatialindex6 libgeos3.10.2 \
    libffi7 libssl3 libjpeg8 libpng16-16 libtiff5 \
    git wget curl \
    && rm -rf /var/lib/apt/lists/*

# Make python3.11 the default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 \
    && update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1

RUN python -m pip install --no-cache-dir --upgrade pip setuptools

# Install Baysor pre-built binary
RUN apt-get update && apt-get install -y --no-install-recommends unzip \
    && rm -rf /var/lib/apt/lists/* && \
    wget --retry-connrefused --waitretry=5 --tries=3 -O /tmp/baysor.zip \
    https://github.com/kharchenkolab/Baysor/releases/download/v0.7.1/baysor-x86_x64-linux-v0.7.1_build.zip && \
    unzip -q /tmp/baysor.zip -d /tmp/baysor_extracted && \
    find /tmp/baysor_extracted -type f -name "baysor" -executable | head -1 | xargs -I {} cp {} /usr/local/bin/baysor && \
    chmod +x /usr/local/bin/baysor && \
    rm -rf /tmp/baysor.zip /tmp/baysor_extracted

# Install wheels from builder stage
COPY --from=builder /build/wheels /tmp/wheels
RUN pip install --no-cache-dir --no-index --find-links=/tmp/wheels /tmp/wheels/*.whl \
    && rm -rf /tmp/wheels

# Clone Banksy_py (submodule not available in .dockerignore'd .git context)
WORKDIR /app
RUN git clone --depth 1 https://github.com/prabhakarlab/Banksy_py.git /app/Banksy_py

# Copy pipeline code (self-contained, no xb/ or notebooks/ dependency)
COPY pipeline/ /app/pipeline/

# Copy entrypoint
COPY docker/entrypoint.sh /app/entrypoint.sh
RUN chmod +x /app/entrypoint.sh

# Add Banksy_py to Python path
ENV PYTHONPATH="/app:${PYTHONPATH:-}"

WORKDIR /app
ENTRYPOINT ["/app/entrypoint.sh"]
CMD []
