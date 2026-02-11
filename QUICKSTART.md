# Quick Start Guide

## Installation

1. **Install Nextflow**:
   ```bash
   curl -s https://get.nextflow.io | bash
   chmod +x nextflow
   sudo mv nextflow /usr/local/bin/
   ```

2. **Install Python dependencies**:
   ```bash
   pip install scanpy matplotlib pandas openpyxl igraph leidenalg
   ```

## Running the Pipeline

### Option 1: Use the example script

```bash
cd /home/workspace/scanpytest/nextflow-scanpy
./run_example.sh
```

### Option 2: Run manually

```bash
cd /home/workspace/scanpytest/nextflow-scanpy

nextflow run main.nf \
    --input /home/workspace/scanpytest \
    --outdir ./results \
    --min_genes 200 \
    --max_genes 2500 \
    --max_mt_percent 5 \
    --annotate_celltypes true
```

### Option 3: Resume a previous run

```bash
nextflow run main.nf \
    --input /home/workspace/scanpytest \
    --outdir ./results \
    -resume
```

## Checking Pipeline Status

While the pipeline is running, you can monitor progress:

```bash
# View real-time logs
tail -f .nextflow.log

# Check what's running
nextflow log

# View work directory
ls -lh work/
```

## Understanding Output

After completion, check:
- `results/analysis_summary.html` - Main report
- `results/07_visualization/` - All plots
- `results/pipeline_info/` - Execution reports

## Customizing for Your Data

1. **Adjust QC thresholds**:
   ```bash
   nextflow run main.nf --input data --max_mt_percent 10 --max_genes 3000
   ```

2. **Change clustering resolution** (higher = more clusters):
   ```bash
   nextflow run main.nf --input data --leiden_resolution 1.5
   ```

3. **Skip cell type annotation**:
   ```bash
   nextflow run main.nf --input data --annotate_celltypes false
   ```

## Common Issues

**Issue**: `igraph` not found  
**Solution**: `pip install igraph leidenalg`

**Issue**: Out of memory  
**Solution**: Increase memory: `--max_memory '32.GB'`

**Issue**: Pipeline fails midway  
**Solution**: Use `-resume` to continue from checkpoint

## Next Steps

1. Review the analysis summary report
2. Examine UMAPs and clustering results
3. Validate cell type annotations
4. Customize `modules/annotate.nf` for your specific cell types
5. Run differential expression analysis on identified cell types
