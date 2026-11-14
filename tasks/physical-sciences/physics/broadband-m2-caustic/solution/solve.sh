#!/bin/bash
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p /app/output
cd /solution
python /solution/m2_fit.py /app/data /app/output/results.json
