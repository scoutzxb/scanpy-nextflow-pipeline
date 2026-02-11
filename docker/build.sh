#!/bin/bash
# Build Docker image for scanpy pipeline

set -e

IMAGE_NAME="scanpy-pipeline"
IMAGE_TAG="1.0.0"

echo "Building Docker image: ${IMAGE_NAME}:${IMAGE_TAG}"

cd "$(dirname "$0")/.."

docker build \
    -t ${IMAGE_NAME}:${IMAGE_TAG} \
    -t ${IMAGE_NAME}:latest \
    -f docker/Dockerfile \
    .

echo ""
echo "Build complete!"
echo "Image: ${IMAGE_NAME}:${IMAGE_TAG}"
echo ""
echo "Test with:"
echo "  docker run --rm ${IMAGE_NAME}:${IMAGE_TAG} python -c 'import scanpy; print(scanpy.__version__)'"
echo ""
echo "Run pipeline with:"
echo "  nextflow run main.nf -profile docker"
