#!/bin/bash

#SBATCH --job-name=NCPF_Calculation_0707
#SBATCH --output=logs/NPCF_Calculation_0707.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jamessunseri@berkeley.edu
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=32gb
#SBATCH --account=zslepian
#SBATCH --qos=zslepian

module load python
pip install tqdm
pip install sarabande

echo running 4pcf script
python calc_4PCF_script.py -M_s 0.7 -M_A 0.7 

