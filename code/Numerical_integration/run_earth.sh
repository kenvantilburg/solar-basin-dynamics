#!/bin/bash
#SBATCH --job-name=rkn12backwards
#SBATCH -N1 --exclusive --ntasks-per-node=128 -t 06-23:00:00
#SBATCH -C rome
#SBATCH --mem-per-cpu=5000
#SBATCH --mail-user=robertlasenby@gmail.com
#SBATCH --mail-type=ALL

module load modules-nix
module load nix/julia/1.5.3

parallel="./parallel --delay .2 -j 128 --joblog joblog_earth"
$parallel "julia driver4.jl EARTH {1} 4.5e9 1e3 outearth" ::: {1..128}
