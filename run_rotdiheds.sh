#!/bin/bash
#SBATCH -t 14-00:00:00   			      #Time for the job to run
#SBATCH --job-name=RunRotDiheds           #Name of the job
#SBATCH -N 1 	        			       	#Number of nodes required
#SBATCH -n 32				                #Number of cores needed for the job
#SBATCH -p SKY32M192_L
#SBATCH --account=col_cmri235_uksr  #Name of account to run under

#SBATCH --mail-type ALL			      	#Send email on start/end
#SBATCH --mail-user rdu230@uky.edu	#Where to send email
#SBATCH --error=SLURM_RunRotDiheds_%j.err		#Name of error file
#SBATCH --output=SLURM_RunRotDiheds_%j.out 	#Name of output file
#SBATCH --array=1-100%30

home=$(pwd)

cd rotated_dihed
for m in mols*/;
do
	if [ -d "${m}"]; then
		cd ${m}
		for d in mols*deg/;
		do
			cd ${d}
			for file in ./*.log
			do
				if [ -f "${file}" ]; then
				break
				else
				sbatch *.job
				fi
			done
			cd ..
		done
		cd ..
	fi
done
