#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=williamson.v@ufl.edu
#SBATCH --time=80:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --account=zslepian
#SBATCH --qos=zslepian
#SBATCH --output=slurm_%a.out # stdout file
#SBATCH --error=slurm_%a.err  # stderr file
#SBATCH --array=500,550,600,650,700,750,800,850,900 

echo imports
module load python
pip install sarabande

echo environment
echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

echo running jobs
python calc_4pcf_array.py 
