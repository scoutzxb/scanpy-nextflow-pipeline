# Running on HPC Clusters

## Quick Start for HPC

### 1. Clone the repository on your HPC
```bash
git clone https://github.com/scoutzxb/scanpy-nextflow-pipeline.git
cd scanpy-nextflow-pipeline
```

### 2. Install Nextflow (if not already available)
```bash
# Download Nextflow
curl -s https://get.nextflow.io | bash

# Move to your bin directory
mkdir -p ~/bin
mv nextflow ~/bin/
export PATH=$HOME/bin:$PATH

# Verify installation
nextflow -version
```

### 3. Choose your execution method

#### Option A: Conda (Recommended for HPC)
```bash
# Load conda module (if needed)
module load conda  # or miniconda3, anaconda3, etc.

# Run with conda profile
nextflow run main.nf \
    -profile conda \
    --input /path/to/10x_data \
    --outdir results \
    -resume
```

#### Option B: Singularity (Best for shared HPC)
```bash
# Load singularity module
module load singularity

# Run with singularity profile
nextflow run main.nf \
    -profile singularity \
    --input /path/to/10x_data \
    --outdir results \
    -resume
```

#### Option C: SLURM with Conda
```bash
# Use SLURM scheduler
nextflow run main.nf \
    -profile slurm_conda \
    --input /path/to/10x_data \
    --outdir results \
    -resume
```

## HPC-Specific Configuration

### Custom SLURM Configuration

Create a custom config file `hpc.config`:

```groovy
// Custom HPC configuration
params {
    max_cpus = 16
    max_memory = '64.GB'
    max_time = '24.h'
}

process {
    executor = 'slurm'
    queue = 'normal'  // Change to your queue name
    
    // Account/partition settings
    clusterOptions = '--account=your_account'
    
    // Process-specific resources
    withName: 'LOAD_AND_QC' {
        cpus = 4
        memory = '16.GB'
        time = '2.h'
    }
    
    withName: 'MARKER_GENES' {
        cpus = 8
        memory = '32.GB'
        time = '4.h'
    }
}

singularity {
    enabled = true
    autoMounts = true
    cacheDir = '/path/to/singularity/cache'
}
```

Run with custom config:
```bash
nextflow run main.nf \
    -c hpc.config \
    --input /path/to/10x_data \
    --outdir results
```

## Common HPC Scenarios

### Scenario 1: Large Dataset (>50K cells)

```bash
nextflow run main.nf \
    -profile slurm_conda \
    --input /path/to/large_data \
    --outdir results \
    --max_cpus 32 \
    --max_memory '128.GB' \
    --max_time '12.h' \
    -resume
```

### Scenario 2: Multiple Samples (batch processing)

Create a sample sheet `samples.txt`:
```
sample1,/path/to/sample1_10x
sample2,/path/to/sample2_10x
sample3,/path/to/sample3_10x
```

Then process in a loop:
```bash
#!/bin/bash
while IFS=',' read -r sample path; do
    nextflow run main.nf \
        -profile slurm_conda \
        --input "$path" \
        --outdir "results_${sample}" \
        -resume
done < samples.txt
```

### Scenario 3: Interactive Session Testing

```bash
# Request interactive node
salloc -n 4 -t 2:00:00 --mem=16G

# Run pipeline
nextflow run main.nf \
    --input /path/to/test_data \
    --outdir test_results \
    -resume

# Exit when done
exit
```

## Monitoring Your Pipeline

### Check pipeline status
```bash
# View running jobs
squeue -u $USER

# Monitor Nextflow log
tail -f .nextflow.log

# Check specific process
ls -lh work/
```

### View resource usage
```bash
# After completion, check the report
firefox results/pipeline_info/execution_report.html
```

## Environment Setup

### Method 1: Pre-build Conda Environment (Faster startup)

```bash
# Create environment once
conda env create -f conda/environment.yml -n scanpy-pipeline

# Run pipeline using existing environment
conda activate scanpy-pipeline
nextflow run main.nf --input data --outdir results
```

### Method 2: Singularity Container (Most reproducible)

```bash
# Pull/build container once
singularity pull docker://python:3.12-slim
# Or build from Dockerfile:
singularity build scanpy-pipeline.sif docker/Dockerfile

# Run with container
nextflow run main.nf \
    -profile singularity \
    --input data \
    --outdir results
```

## Troubleshooting on HPC

### Issue: "Command not found: conda"
```bash
# Load conda module
module avail conda  # Find available conda modules
module load conda/latest  # Load conda
```

### Issue: "Singularity not found"
```bash
# Load singularity module
module avail singularity
module load singularity/latest
```

### Issue: Out of disk space in home directory
```bash
# Use scratch space for work directory
export NXF_WORK=/scratch/$USER/nextflow-work

# Run pipeline
nextflow run main.nf --input data --outdir results
```

### Issue: Jobs killed by scheduler
```bash
# Increase resource limits
nextflow run main.nf \
    --max_memory '64.GB' \
    --max_time '12.h' \
    --input data \
    --outdir results
```

### Issue: Conda environment creation is slow
```bash
# Use mamba (faster conda)
conda install mamba -c conda-forge

# Modify nextflow.config to use mamba:
# conda.useMamba = true
```

## Performance Tips

1. **Use `-resume`**: Always use `-resume` to continue from checkpoints
2. **Scratch space**: Run on fast scratch filesystem
3. **Cache directories**: Set cache dirs on fast storage
   ```bash
   export NXF_CONDA_CACHEDIR=/scratch/$USER/conda
   export NXF_SINGULARITY_CACHEDIR=/scratch/$USER/singularity
   ```
4. **Parallel jobs**: Nextflow automatically parallelizes independent tasks
5. **Monitor resources**: Check execution report to optimize future runs

## Example SLURM Submission Script

Create `run_pipeline.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=scanpy-pipeline
#SBATCH --output=scanpy-%j.out
#SBATCH --error=scanpy-%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=normal

# Load modules
module load conda
module load nextflow

# Set working directories
export NXF_WORK=/scratch/$USER/work
export NXF_CONDA_CACHEDIR=/scratch/$USER/conda

# Run pipeline
nextflow run main.nf \
    -profile slurm_conda \
    --input /path/to/10x_data \
    --outdir /path/to/results \
    --max_cpus 16 \
    --max_memory '64.GB' \
    -resume

echo "Pipeline complete!"
```

Submit with:
```bash
sbatch run_pipeline.sh
```

## Getting Help

1. Check Nextflow log: `.nextflow.log`
2. Check process logs: `work/*/command.log`
3. Check SLURM logs: `scanpy-*.out` and `scanpy-*.err`
4. Contact your HPC support for cluster-specific issues
5. Open an issue on GitHub for pipeline-related questions

## Best Practices for HPC

✅ Always use `-resume` to save time
✅ Use scratch/fast storage for work directory
✅ Pre-build conda environments for repeated use
✅ Monitor resource usage and adjust parameters
✅ Clean up work directory after successful runs
✅ Keep containers/environments in shared locations for team
✅ Document your HPC-specific configuration
