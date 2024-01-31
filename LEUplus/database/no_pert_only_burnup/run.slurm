#!/bin/bash

#SBATCH -J serpentXSGen                         # Job name
#SBATCH -A fillInYourChargeCode                          # account to which job is charged, ex: GT-gburdell3
#SBATCH -N126 --ntasks-per-node=6                 # Number of nodes and cores per node
#SBATCH --cpus-per-task=4
#SBATCH --time=7:30:00                          # Walltime = 7.5h
#SBATCH -qembers                                # QOS Name (embers / inferno)
#SBATCH -oReport-%j.out                         # Combined output and error messages file

echo $SLURM_SUBMIT_DIR
cd   $SLURM_SUBMIT_DIR                            # change into directory from where job to be executed (where data is / script is located)

source ~/.bashrc                                # load environment
echo "Started on `/bin/hostname`"               # prints the name of the node job started on
srun python run.py
