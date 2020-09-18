#!/bin/bash
#SBATCH -t 14-00:00:00   			      #Time for the job to run
#SBATCH --job-name=submitRotDiheds           #Name of the job
#SBATCH -N 1 	        			       	#Number of nodes required
#SBATCH -n 32				                #Number of cores needed for the job
#SBATCH -p SKY32M192_L
#SBATCH --account=col_cmri235_uksr  #Name of account to run under

#SBATCH --mail-type ALL			      	#Send email on start/end
#SBATCH --mail-user rdu230@uky.edu	#Where to send email
#SBATCH --error=SLURM_submitRotDiheds_%j.err		#Name of error file
#SBATCH --output=SLURM_submitRotDiheds_%j.out 	#Name of output file
#SBATCH --array=1-100%30

module load ccs/gaussian/g16-A.03/g16-sandybridge
ulimit -u 32768
echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

home=$(pwd)

cd rotated_dihed
for m in mols*/
do
	if [ -d "${m}" ]; then
		cd ${m}
		for d in mols*deg/
		do
			cd ${d}
			for file in ./
			do
				if [[ ${file: -3} == "deg" ]]; then break 2; fi
			done
			fn="${d%?}.gjf"
			g16 ${fn}
			cd ..
		done
		cd ..
	fi
done
