#!/bin/bash

#SBATCH --job-name=NCPF_Calculation
#SBATCH --output=NPCF_Calculation.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=williamson.v@ufl.edu
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --account=zslepian
#SBATCH --qos=zslepian
module load python
pip install sarabande==1.0.0

python test_run.py -ell_max 4 -n_bins 20 -t_slice 500 -Ma 0.7 -Ms 0.7

