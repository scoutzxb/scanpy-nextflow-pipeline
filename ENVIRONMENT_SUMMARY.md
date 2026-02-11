# Environment Management - Quick Summary

## How the Pipeline Handles Environment

The Nextflow pipeline supports **4 execution modes** via profiles:

```bash
# 1. Standard (default) - uses your current Python environment
nextflow run main.nf --input data --outdir results

# 2. Conda - Nextflow creates/manages conda environment automatically
nextflow run main.nf -profile conda --input data --outdir results

# 3. Docker - runs each process in a Docker container
nextflow run main.nf -profile docker --input data --outdir results

# 4. Singularity - uses Singularity containers (for HPC)
nextflow run main.nf -profile singularity --input data --outdir results
```

## What's Provided

### ✅ Conda Environment Definition
- `conda/environment.yml` - complete conda environment spec
- Pin
[truncated]
/build.sh` - automated build script

### ✅ Nextflow Profiles
- `nextflow.config` - 6 execution profiles configured:
  - `standard` - local Python
  - `conda` - conda-managed environment
  - `docker` - Docker containers
  - `singularity` - Singularity containers
  - `slurm_conda` - SLURM + conda
  - `slurm_singularity` - SLURM + singularity

### ✅ Comprehensive Documentation
- `ENVIRONMENT.md` - full 300+ line guide
- Setup instructions for all modes
- Troubleshooting tips
- Best practices

## Quick Start Options

### For Quick Testing (Current Setup)
```bash
# Already have scanpy installed
nextflow run main.nf --input /home/workspace/scanpytest --outdir results
```

### For Reproducible Analysis (Recommended)
```bash
# Let Nextflow handle everything
nextflow run main.nf -profile conda --input /home/workspace/scanpytest --outdir results
# First run creates conda environment (takes 5-10 min)
# Subsequent runs use cached environment (fast)
```

### For HPC Clusters
```bash
# Singularity is standard on most HPC
nextflow run main.nf -profile slurm_singularity \
    --input /data/scrna --outdir /results
```

### For Maximum Reproducibility
```bash
# Build Docker container once
bash docker/build.sh

# Run with Docker
nextflow run main.nf -profile docker --input data --outdir results
```

## Comparison at a Glance

| Mode | Setup Time | Reproducibility | Use Case |
|------|-----------|----------------|----------|
| **Standard** | 0 min (if installed) | ⭐⭐ | Quick tests, development |
| **Conda** | 5-10 min (first time) | ⭐⭐⭐⭐ | Most users, production |
| **Docker** | 10-15 min (build) | ⭐⭐⭐⭐⭐ | Publications, sharing |
| **Singularity** | 10-15 min (convert) | ⭐⭐⭐⭐⭐ | HPC clusters |

## Key Advantages

### Conda Profile (Recommended)
✅ Nextflow automatically creates and manages the environment  
✅ No manual conda commands needed  
✅ Environment cached between runs  
✅ Works across platforms  
✅ Easy to update (just update environment.yml)

### Docker/Singularity
✅ Completely isolated environment  
✅ Bit-for-bit reproducible  
✅ Share exact same environment  
✅ Include in publications  
✅ Long-term archival

## What Happens Behind the Scenes

### With `-profile conda`:
1. Nextflow reads `conda/environment.yml`
2. Creates conda environment (cached in `~/.nextflow/conda/`)
3. Runs each process in that environment
4. Environment persists for future runs

### With `-profile docker`:
1. Nextflow checks for `scanpy-pipeline:1.0.0` image
2. Pulls/builds if not present
3. Runs each process in container
4. Mounts input/output directories automatically

### With `-profile singularity`:
1. Converts Docker image to Singularity format (if needed)
2. Caches .sif file in `~/.nextflow/singularity/`
3. Runs each process in container
4. Works without root permissions

## Files Created

```
nextflow-scanpy/
├── conda/
│   └── environment.yml          # Conda environment specification
├── docker/
│   ├── Dockerfile              # Docker image definition
│   ├── requirements.txt        # Python packages
│   └── build.sh               # Build script
├── nextflow.config             # Includes all profiles
└── ENVIRONMENT.md              # Full documentation
```

## Bottom Line

**Your pipeline now has professional-grade environment management:**

1. **Current state**: Works with system Python (already functional)
2. **Best practice**: Add `-profile conda` for reproducibility
3. **Publication**: Use `-profile docker` for maximum reproducibility
4. **HPC**: Use `-profile slurm_singularity` for cluster computing

**No code changes needed** - just add the profile flag!
