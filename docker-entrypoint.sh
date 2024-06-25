#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
echo immunotyper-SR --output_dir /output "$@"
exec immunotyper-SR --output_dir /output "$@"