#!/bin/bash

#SBATCH --job-name=NCPF_Calculation
#SBATCH --output=logs/NPCF_Calculation.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=williamson.v@ufl.edu
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16gb
#SBATCH --account=zslepian
#SBATCH --qos=zslepian

echo imports
module load python
pip install sarabande

echo running 4pcf script
python calc_4pcf.py 
