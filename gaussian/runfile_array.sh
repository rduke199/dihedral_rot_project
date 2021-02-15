#!/bin/bash
#SBATCH -t 14-00:00:00   			      #Time for the job to run
#SBATCH --job-name=submitRotDiheds           #Name of the job
#SBATCH -N 1 	        			       	#Number of nodes required
#SBATCH -n 8				                #Number of cores needed for the job
#SBATCH -p SKY32M192_L
#SBATCH --account=col_cmri235_uksr  #Name of account to run under

#SBATCH --mail-type ALL			      	#Send email on start/end
#SBATCH --mail-user rdu230@uky.edu	#Where to send email
#SBATCH --error=SLURM_submitRotDiheds_%j.err		#Name of error file
#SBATCH --output=SLURM_submitRotDiheds_%j.out 	#Name of output file
#SBATCH --array=1-600%30

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "
echo "Job started at : $(date)"

module load ccs/gaussian/g16-A.03/g16-sandybridge
ulimit -s unlimited
ulimit -l unlimited
export OMP_NUM_THREADS=1

line=$(echo $SLURM_ARRAY_TASK_ID )
runfolder=$(sed -n "$line"p folders_to_run.txt)
echo $runfolder
cd $runfolder
time g16 *.gjf
echo "Job ended at : $(date)"
