# Environment Management Guide

## Current Implementation

⚠️ **Important**: The current pipeline (v1.0) assumes Python and scanpy are already installed in your system environment. This works but has limitations for reproducibility.

## Environment Options

### Option 1: System Python (Current - Basic)

**Pros:**
- ✅ Fast execution (no container overhead)
- ✅ Works immediately if deps are installed
- ✅ Easy to debug

**Cons:**
- ❌ Not fully reproducible
- ❌ Version conflicts possible
- ❌ Manual dependency management

**Setup:**
```bash
pip install scanpy igraph leidenalg openpyxl
```

### Option 2: Conda Environment (Recommended)

**Pros:**
- ✅ Isolated environment
- ✅ Version pinning
- ✅ Easy to share
- ✅ Cross-platform

**Setup:**
See `conda/environment.yml` (I'll create this)

### Option 3: Docker Container (Best for Reproducibility)

**Pros:**
- ✅ Fully reproducible
- ✅ No local installation needed
- ✅ Portable across systems
- ✅ Version locked

**Setup:**
See `docker/Dockerfile` (I'll create this)

### Option 4: Singularity (Best for HPC)

**Pros:**
- ✅ Works on HPC clusters
- ✅ No root required
- ✅ Fully reproducible
- ✅ Docker compatible

**Setup:**
Can convert from Docker image

## Detailed Setup Instructions Below

---

## Option 1: System Python (Quick Start)

### Requirements
- Python 3.10 or higher
- pip package manager

### Installation
```bash
pip install scanpy==1.10.1 igraph==0.11.4 leidenalg==0.10.2 openpyxl==3.1.2
```

### Run Pipeline
```bash
nextflow run main.nf --input /path/to/data --outdir results
# or
nextflow run main.nf -profile standard --input /path/to/data --outdir results
```

### Verification
```bash
python -c "import scanpy; print(f'Scanpy: {scanpy.__version__}')"
python -c "import igraph; print(f'igraph: {igraph.__version__}')"
python -c "import leidenalg; print(f'leidenalg: {leidenalg.__version__}')"
```

---

## Option 2: Conda Environment (Recommended for Most Users)

### Requirements
- Conda or Mamba installed
- Get conda: https://docs.conda.io/en/latest/miniconda.html

### Installation
```bash
# Create environment from file
cd /home/workspace/scanpytest/nextflow-scanpy
conda env create -f conda/environment.yml

# Activate environment
conda activate scanpy-pipeline

# Verify
python -c "import scanpy; print(scanpy.__version__)"
```

### Run Pipeline
```bash
# Option A: With activated environment
conda activate scanpy-pipeline
nextflow run main.nf --input /path/to/data --outdir results

# Option B: Let Nextflow manage conda (recommended)
nextflow run main.nf -profile conda --input /path/to/data --outdir results
```

### Advantages
- ✅ Isolated from system Python
- ✅ Reproducible across machines
- ✅ Easy to share (just share environment.yml)
- ✅ Works on Linux, macOS, Windows

### Update Environment
```bash
conda env update -f conda/environment.yml --prune
```

### Remove Environment
```bash
conda env remove -n scanpy-pipeline
```

---

## Option 3: Docker (Best for Reproducibility)

### Requirements
- Docker installed
- Get Docker: https://docs.docker.com/get-docker/

### Build Container
```bash
cd /home/workspace/scanpytest/nextflow-scanpy
bash docker/build.sh
```

### Or Pull Pre-built (if available)
```bash
docker pull scanpy-pipeline:1.0.0
```

### Run Pipeline
```bash
nextflow run main.nf -profile docker --input /path/to/data --outdir results
```

### Manual Docker Run (for testing)
```bash
# Test the container
docker run --rm scanpy-pipeline:1.0.0 python -c "import scanpy; print(scanpy.__version__)"

# Interactive session
docker run -it --rm -v $(pwd):/workspace scanpy-pipeline:1.0.0 bash

# Run a script
docker run --rm \
    -v /path/to/data:/data \
    -v /path/to/results:/results \
    scanpy-pipeline:1.0.0 \
    python process_script.py
```

### Advantages
- ✅ 100% reproducible
- ✅ Same environment everywhere
- ✅ No dependency conflicts
- ✅ Easy to share (push to Docker Hub)
- ✅ Version locked

### Push to Registry
```bash
# Tag for registry
docker tag scanpy-pipeline:1.0.0 youruser/scanpy-pipeline:1.0.0

# Push
docker push youruser/scanpy-pipeline:1.0.0

# Update nextflow.config to use your image
# process.container = 'youruser/scanpy-pipeline:1.0.0'
```

---

## Option 4: Singularity (HPC Clusters)

### Requirements
- Singularity/Apptainer installed (common on HPC)
- Or ability to convert Docker images

### Build from Docker
```bash
# Convert Docker image to Singularity
singularity build scanpy-pipeline.sif docker://scanpy-pipeline:1.0.0

# Or from Docker Hub
singularity build scanpy-pipeline.sif docker://youruser/scanpy-pipeline:1.0.0
```

### Run Pipeline
```bash
nextflow run main.nf -profile singularity --input /path/to/data --outdir results
```

### Advantages
- ✅ Works on HPC without root
- ✅ Fully reproducible
- ✅ Compatible with Docker images
- ✅ Better for shared clusters

### Cache Directory
```bash
# Set cache location (important for HPC)
export NXF_SINGULARITY_CACHEDIR="/path/to/cache"
```

---

## Comparison Table

| Feature | System Python | Conda | Docker | Singularity |
|---------|--------------|-------|--------|-------------|
| **Setup Speed** | ⚡ Fastest | ⚡ Fast | 🐌 Slow (first time) | 🐌 Slow (first time) |
| **Reproducibility** | ❌ Low | ✅ Good | ✅✅ Excellent | ✅✅ Excellent |
| **Portability** | ❌ Poor | ✅ Good | ✅✅ Excellent | ✅✅ Excellent |
| **HPC Compatible** | ✅ Yes | ✅ Yes | ❌ Often no | ✅✅ Yes |
| **Root Required** | ❌ No | ❌ No | ✅ Sometimes | ❌ No |
| **Disk Space** | ✅ Small | 🔶 Medium | 🔶 Medium | 🔶 Medium |
| **Isolation** | ❌ None | ✅ Process-level | ✅✅ Full | ✅✅ Full |
| **Best For** | Quick tests | Development | Production | HPC clusters |

---

## Recommended Workflow

### For Development/Testing
1. Start with **System Python** for quick iteration
2. Switch to **Conda** once analysis is working

### For Production/Publication
1. Use **Docker** for complete reproducibility
2. Share container via Docker Hub or registry
3. Include container hash in methods section

### For HPC Clusters
1. Use **Singularity** profile
2. Build container once, cache it
3. Share .sif file or Docker image reference

---

## Troubleshooting

### "Command not found: python"
```bash
# Ensure Python is installed
which python
python --version

# Or use python3
which python3
```

### "Module not found: scanpy"
```bash
# Install in current environment
pip install scanpy

# Or activate conda environment
conda activate scanpy-pipeline
```

### Docker: "Permission denied"
```bash
# Add user to docker group (Linux)
sudo usermod -aG docker $USER
# Log out and back in

# Or run with sudo (not recommended)
sudo nextflow run main.nf -profile docker ...
```

### Conda: "Environment not found"
```bash
# List environments
conda env list

# Recreate if needed
conda env create -f conda/environment.yml
```

### Nextflow: "Process terminated with exit code 1"
```bash
# Check specific process work directory
cat work/xx/xxxxxxxxxxxx/.command.err

# Test Python environment
python -c "import scanpy, igraph, leidenalg"
```

---

## Environment Variables

### For Conda
```bash
export CONDA_PKGS_DIRS=/path/to/conda/cache
export CONDA_ENVS_PATH=/path/to/conda/envs
```

### For Docker
```bash
export DOCKER_TMPDIR=/path/to/tmp
export DOCKER_HOST=unix:///var/run/docker.sock
```

### For Singularity
```bash
export SINGULARITY_CACHEDIR=/path/to/cache
export SINGULARITY_TMPDIR=/path/to/tmp
```

### For Nextflow
```bash
export NXF_CONDA_CACHEDIR=/path/to/conda/cache
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity/cache
```

---

## Best Practices

1. **Version Control**: Always specify exact versions in requirements
2. **Testing**: Test environment on small dataset before large runs
3. **Documentation**: Note which environment/profile was used
4. **Sharing**: Share environment.yml or Dockerfile with publications
5. **Archiving**: Keep container images with data for long-term storage

---

## Future Enhancements (Roadmap)

- **v1.1**: Pre-built Docker images on Docker Hub
- **v1.2**: Bioconda package for easy installation
- **v1.3**: Cloud-native containers (AWS, GCP, Azure)
- **v1.4**: GPU-accelerated containers for large datasets