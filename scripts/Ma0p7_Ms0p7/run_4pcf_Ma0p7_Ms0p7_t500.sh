#!/bin/bash

#SBATCH --job-name=NCPF_Calculation_Ma0p7_Ms0p7_t500
#SBATCH --output=logs/NPCF_Calculation_Ma0p7_Ms0p7_t500.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=williamson.v@ufl.edu
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --account=zslepian
#SBATCH --qos=zslepian
module load python
pip install sarabande==1.0.1

python calc_4pcf.py -ell_max 1 -n_bins 20 -t_slice 500 -Ma 0.7 -Ms 0.7

