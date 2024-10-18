#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate immunotyper-SR
echo immunotyper-SR --output_dir /output "$@"
exec immunotyper-SR --output_dir /output "$@"