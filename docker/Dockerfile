# Dockerfile for Scanpy Single-Cell RNA-seq Pipeline
# Usage: docker build -t scanpy-pipeline:1.0.0 -f docker/Dockerfile .

FROM python:3.12-slim

LABEL maintainer="Xiaobin Zheng"
LABEL description="Scanpy single-cell RNA-seq analysis environment"
LABEL version="1.0.0"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    gfortran \
    libhdf5-dev \
    libopenblas-dev \
    liblapack-dev \
    pkg-config \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /workspace

# Copy requirements file
COPY docker/requirements.txt /tmp/requirements.txt

# Install Python packages
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

# Verify installation
RUN python -c "import scanpy; print(f'Scanpy version: {scanpy.__version__}')" && \
    python -c "import igraph; print(f'igraph version: {igraph.__version__}')" && \
    python -c "import leidenalg; print(f'leidenalg version: {leidenalg.__version__}')"

# Set matplotlib backend
ENV MPLBACKEND=Agg

# Default command
CMD ["/bin/bash"]
