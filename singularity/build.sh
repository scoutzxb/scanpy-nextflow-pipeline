#!/bin/bash
# Build Singularity container for scanpy pipeline

set -e

CONTAINER_NAME="scanpy-pipeline_1.0.0.sif"
DEF_FILE="singularity/scanpy-pipeline.def"

echo "Building Singularity container..."
echo "This requires sudo/root privileges or --fakeroot flag"
echo ""

# Check if definition file exists
if [ ! -f "$DEF_FILE" ]; then
    echo "Error: Definition file not found: $DEF_FILE"
    exit 1
fi

# Build the container
echo "Building from: $DEF_FILE"
echo "Output: $CONTAINER_NAME"
echo ""

# Try different build methods
if command -v singularity &> /dev/null; then
    # Method 1: With sudo (if available)
    if sudo -n true 2>/dev/null; then
        echo "Building with sudo..."
        sudo singularity build "$CONTAINER_NAME" "$DEF_FILE"
    # Method 2: With fakeroot (if available)
    elif singularity version | grep -q "fakeroot"; then
        echo "Building with fakeroot..."
        singularity build --fakeroot "$CONTAINER_NAME" "$DEF_FILE"
    # Method 3: Remote build (requires Sylabs account)
    else
        echo "No sudo or fakeroot available. Trying remote build..."
        echo "This requires a Sylabs Cloud account."
        singularity build --remote "$CONTAINER_NAME" "$DEF_FILE"
    fi
else
    echo "Error: Singularity is not installed"
    exit 1
fi

echo ""
echo "Build complete! Container saved to: $CONTAINER_NAME"
echo ""
echo "Test the container:"
echo "  singularity exec $CONTAINER_NAME python -c 'import scanpy; print(scanpy.__version__)'"
echo ""
echo "Run pipeline with:"
echo "  nextflow run main.nf -profile singularity --input data --outdir results"
