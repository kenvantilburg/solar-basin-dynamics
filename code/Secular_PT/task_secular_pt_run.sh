#!/bin/bash

ic=$1

# Activate conda environment
source ~/.bashrc
conda activate ss_sim

echo "Starting secular PT run at core $ic ..."
python3 secular_pt_run.py
echo "Finished secular PT run at core $ic."