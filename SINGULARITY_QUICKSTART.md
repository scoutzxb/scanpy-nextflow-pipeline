# Singularity Quick Start - 3 Simple Steps

## Step 1: Build the Container

Choose ONE method:

### Method A: With Fakeroot (Most Common on HPC)
```bash
cd scanpy-nextflow-pipeline
singularity build --fakeroot scanpy-pipeline_1.0.0.sif singularity/scanpy-pipeline.def
```

### Method B: Remote Build (No Local Privileges)
```bash
singularity remote login  # First time only
singularity build --remote scanpy-pipeline_1.0.0.sif singularity/scanpy-pipeline.def
```

### Method C: Build on Local Machine, Transfer to HPC
```bash
# On your laptop/desktop (with sudo)
singularity build scanpy-pipeline_1.0.0.sif singularity/scanpy-pipeline.def

# Transfer to HPC
scp scanpy-pipeline_1.0.0.sif username@hpc.edu:~/scanpy-nextflow-pipeline/
```

## Step 2: Test the Container

```bash
# Quick test
singularity exec scanpy-pipeline_1.0.0.sif python -c "import scanpy; print(scanpy.__version__)"

# Should output: 1.10.1
```

## Step 3: Run the Pipeline

### Option A: Singularity Only (Local Execution)
```bash
nextflow run main.nf \
    -profile singularity \
    --input /path/to/10x_data \
    --outdir results \
    -resume
```

### Option B: Singularity + SLURM (HPC Cluster)
```bash
nextflow run main.nf \
    -profile slurm_singularity \
    --input /path/to/10x_data \
    --outdir results \
    -resume
```

### Option C: Custom Container Path
```bash
nextflow run main.nf \
    -profile singularity \
    --container /full/path/to/scanpy-pipeline_1.0.0.sif \
    --input /path/to/10x_data \
    --outdir results
```

---

## Common Scenarios

### If Your HPC Has Shared Containers
```bash
# Use the shared container
export CONTAINER=/shared/containers/scanpy-pipeline_1.0.0.sif

nextflow run main.nf \
    -profile singularity \
    --container $CONTAINER \
    --input data --outdir results
```

### If Building Fails
Try building from Docker instead:
```bash
# Pull Python base image and install manually
singularity pull docker://python:3.12-slim
# Then modify to add scanpy...

# Or ask your HPC admin to build it for you
```

### If You Get "Container Not Found"
Make sure you're in the pipeline directory:
```bash
cd scanpy-nextflow-pipeline
ls scanpy-pipeline_1.0.0.sif  # Should exist
```

Or specify full path:
```bash
--container /full/path/to/scanpy-pipeline_1.0.0.sif
```

---

## Complete Example

```bash
# 1. Clone repository
git clone https://github.com/scoutzxb/scanpy-nextflow-pipeline.git
cd scanpy-nextflow-pipeline

# 2. Build container
singularity build --fakeroot scanpy-pipeline_1.0.0.sif singularity/scanpy-pipeline.def

# 3. Test it
singularity exec scanpy-pipeline_1.0.0.sif python -c "import scanpy; print('Success!')"

# 4. Run pipeline
nextflow run main.nf \
    -profile slurm_singularity \
    --input /data/10x_pbmc \
    --outdir results \
    --max_cpus 8 \
    --max_memory '32.GB' \
    -resume

# 5. Check results
ls results/
```

---

## Troubleshooting

**"Permission denied" during build**
→ Use `--fakeroot` or `--remote` flag

**"Cannot find container"**
→ Use full path: `--container $(pwd)/scanpy-pipeline_1.0.0.sif`

**"Library does not exist"**
→ Container should be self-contained, this shouldn't happen. Rebuild if it does.

**Build takes forever**
→ Use remote build with Sylabs Cloud

---

For more details, see: **SINGULARITY_GUIDE.md**
