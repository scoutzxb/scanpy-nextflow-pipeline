# Singularity Guide for HPC

## Overview

Singularity is the preferred containerization solution for HPC environments because:
- ✅ Works without root privileges (once container is built)
- ✅ Integrates with HPC schedulers (SLURM, PBS, etc.)
- ✅ No daemon required
- ✅ Better security model for shared systems
- ✅ Can convert Docker containers

## Quick Start

### Option 1: Build from Docker Image (Easiest)

```bash
# Pull and convert the Docker image to Singularity
singularity pull docker://python:3.12-slim

# Or build directly from our Dockerfile
singularity build scanpy-pipeline.sif docker-daemon://scanpy-pipeline:1.0.0
```

### Option 2: Build from Definition File (Recommended)

```bash
cd scanpy-nextflow-pipeline

# Build the container (requires sudo or fakeroot)
singularity build scanpy-pipeline_1.0.0.sif singularity/scanpy-pipeline.def

# Or use the build script
./singularity/build.sh
```

### Option 3: Use Pre-built Container (if available)

```bash
# Download from Sylabs Cloud (if uploaded)
singularity pull library://scoutzxb/default/scanpy-pipeline:1.0.0

# Or from a shared location on your HPC
cp /shared/containers/scanpy-pipeline_1.0.0.sif .
```

## Running the Pipeline with Singularity

### Basic Usage

```bash
nextflow run main.nf \
    -profile singularity \
    --input /path/to/10x_data \
    --outdir results \
    -resume
```

### With SLURM Scheduler

```bash
nextflow run main.nf \
    -profile slurm_singularity \
    --input /path/to/10x_data \
    --outdir results \
    -resume
```

### Custom Container Path

If you've built or downloaded the container to a specific location:

```bash
nextflow run main.nf \
    -profile singularity \
    --container /path/to/scanpy-pipeline_1.0.0.sif \
    --input /path/to/10x_data \
    --outdir results
```

## Building the Container on HPC

### Method 1: With sudo (if you have admin access)

```bash
sudo singularity build scanpy-pipeline.sif singularity/scanpy-pipeline.def
```

### Method 2: With fakeroot (most HPC systems)

```bash
# Check if fakeroot is available
singularity version | grep fakeroot

# Build with fakeroot
singularity build --fakeroot scanpy-pipeline.sif singularity/scanpy-pipeline.def
```

### Method 3: Remote Build (no local privileges needed)

```bash
# Requires Sylabs Cloud account (free)
# Sign up at: https://cloud.sylabs.io/

# Login
singularity remote login

# Build remotely
singularity build --remote scanpy-pipeline.sif singularity/scanpy-pipeline.def
```

### Method 4: Build from Docker (if Docker is available)

```bash
# First build Docker image
cd docker
docker build -t scanpy-pipeline:1.0.0 .

# Convert to Singularity
singularity build scanpy-pipeline.sif docker-daemon://scanpy-pipeline:1.0.0
```

### Method 5: Build on Personal Machine, Transfer to HPC

```bash
# On your local machine (with sudo)
singularity build scanpy-pipeline.sif singularity/scanpy-pipeline.def

# Transfer to HPC
scp scanpy-pipeline.sif username@hpc.cluster.edu:~/containers/
```

## Configuration

### Update nextflow.config (if needed)

The default configuration should work, but you can customize:

```groovy
profiles {
    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
        process.container = 'scanpy-pipeline_1.0.0.sif'
        singularity.cacheDir = "${HOME}/.nextflow/singularity"
    }
}
```

### Environment Variables

```bash
# Set Singularity cache directory
export SINGULARITY_CACHEDIR=$HOME/.singularity/cache

# Set Nextflow Singularity cache
export NXF_SINGULARITY_CACHEDIR=$HOME/.nextflow/singularity

# Bind additional paths (if needed)
export SINGULARITY_BINDPATH="/scratch,/data"
```

## Testing the Container

### Test Python and Scanpy

```bash
singularity exec scanpy-pipeline.sif python -c "import scanpy; print(scanpy.__version__)"
```

### Test All Dependencies

```bash
singularity exec scanpy-pipeline.sif python -c "
import scanpy as sc
import igraph
import leidenalg
import pandas as pd
import numpy as np
print('All packages imported successfully!')
print(f'Scanpy version: {sc.__version__}')
print(f'igraph version: {igraph.__version__}')
"
```

### Interactive Shell

```bash
# Enter container interactively
singularity shell scanpy-pipeline.sif

# Inside container
Singularity> python
>>> import scanpy
>>> scanpy.__version__
```

## Common Issues and Solutions

### Issue: "FATAL: container creation failed: mount error"

**Solution**: Check bind paths and permissions

```bash
# Explicitly bind required directories
export SINGULARITY_BINDPATH="/home,/scratch,/data"
```

### Issue: "WARNING: Skipping mount ... doesn't exist in container"

**Solution**: This is usually harmless, but you can silence it:

```bash
export SINGULARITY_BIND="/scratch:/scratch,/data:/data"
```

### Issue: Cannot build without root

**Solutions**:
1. Use `--fakeroot` flag (if available)
2. Use remote build with Sylabs Cloud
3. Build on your local machine and transfer
4. Ask your HPC admin to build it

### Issue: "library does not exist"

**Solution**: The container is looking for system libraries. Try:

```bash
# Bind system libraries
singularity exec --bind /lib:/lib,/usr/lib:/usr/lib scanpy-pipeline.sif python script.py
```

## Performance Tips

### 1. Use Local Storage for Cache

```bash
# Use node-local storage (faster)
export SINGULARITY_CACHEDIR=/tmp/$USER/singularity
export NXF_SINGULARITY_CACHEDIR=/tmp/$USER/nextflow_singularity
```

### 2. Pre-download Container Images

```bash
# Pull container before running pipeline
singularity pull scanpy-pipeline.sif docker://scanpy-pipeline:1.0.0

# Then use it
nextflow run main.nf -profile singularity --container scanpy-pipeline.sif
```

### 3. Shared Container Location

For team use, store container in shared location:

```bash
# Build once in shared location
SHARED_CONTAINER=/shared/containers/scanpy-pipeline_1.0.0.sif

# Everyone uses the same container
nextflow run main.nf \
    -profile singularity \
    --container $SHARED_CONTAINER \
    --input data --outdir results
```

## Example SLURM Script

```bash
#!/bin/bash
#SBATCH --job-name=scanpy_pipeline
#SBATCH --output=scanpy_%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# Load modules
module load singularity
module load nextflow

# Set cache directories
export SINGULARITY_CACHEDIR=$SCRATCH/.singularity
export NXF_SINGULARITY_CACHEDIR=$SCRATCH/.nextflow/singularity

# Run pipeline
nextflow run main.nf \
    -profile slurm_singularity \
    --input $SCRATCH/data/10x_data \
    --outdir $SCRATCH/results \
    --max_cpus 8 \
    --max_memory '32.GB' \
    -resume
```

## Best Practices

✅ **Build once, use many times**: Share container among team members  
✅ **Version your containers**: Include version in filename (scanpy-pipeline_1.0.0.sif)  
✅ **Test locally first**: Verify container works before large runs  
✅ **Use cache directories**: Set SINGULARITY_CACHEDIR for performance  
✅ **Document container source**: Note how/when container was built  
✅ **Archive containers**: Keep containers with published results for reproducibility

## Additional Resources

- [Singularity Documentation](https://sylabs.io/docs/)
- [Singularity on HPC](https://singularity-tutorial.github.io/)
- [Nextflow with Singularity](https://www.nextflow.io/docs/latest/singularity.html)
- [Sylabs Cloud Library](https://cloud.sylabs.io/library)
