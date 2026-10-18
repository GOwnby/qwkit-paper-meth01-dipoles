#!/bin/bash
# Reproduce every METH-01 number and figure (about 10 min on 2 cores).
set -e
cd "$(dirname "$0")"
export PYTHONPATH=$(cd ../.. && pwd) OMP_NUM_THREADS=1 NPROC=${NPROC:-2}
python3 -W ignore scan.py > scan.log
python3 collect_numbers.py > /dev/null
python3 paper_figures.py
