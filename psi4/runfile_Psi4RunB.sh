#!/bin/bash
#SBATCH -t 14-00:00:00                                #Time for the job to run
#SBATCH --job-name=Psi4RunB          #Name of the job
#SBATCH -N 1                                            #Number of nodes required
#SBATCH -n 32                                           #Number of cores needed for the job
#SBATCH -p SKY32M192_L
#SBATCH --account=col_cmri235_uksr  #Name of account to run under

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user rdu230@uky.edu      #Where to send email
#SBATCH --error=SLURM_Psi4RunB_%j.err         #Name of error file
#SBATCH --output=SLURM_Psi4RunB_%j.out        #Name of output file

ulimit -u 32768
echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "
module load ccs/conda/psi4-1.3.2+ecbda83
. ~/.bashrc

python Psi4RunB.py
