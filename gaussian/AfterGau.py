import os
import re
import shutil
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput, GaussianOutput
from itertools import combinations_with_replacement, product

def FindEnergies():
    """
    Find final energies from all log files in cwd and copy them into table.txt
    """
    a = '0'
    main_dir = os.getcwd()
    for dir in os.listdir(main_dir):
        dir_path = os.path.join(main_dir,dir)
        os.chdir(dir_path)
        for i in os.listdir(dir_path):
            if i.endswith('.log'):
                mol = GaussianOutput(i)
                if mol.normal_termination:
                    Energy = mol.final_energy
                if normal == 1:
                    with open(main_dir+'/table.txt','a') as fout:
                        fout.write(i + '\t' + Energy + '\n')
                else:
                    print('Gaussian for {} did not terminate normally!'.format(i[:-4]))
        os.chdir(main_dir)

def FindBestConfromer(energy_file):
    """
    Find lowest energy confromer of each polymer from table.txt file. Make wtuning
    directory with best confromer files.
    """
    list1 = []
    #RECALL naming in table:
    #mols_<ring type>_<#repeat units>_<polymer#>_<subs>_<confromer#> Done:  E(RPM6) =  <energy> A.U. after   <#cycles>
    with open(energy_file,'r') as fn:
        for line in fn:
            try:
                mol = line.split()[0]
                ring_num = int(mol.split('_')[1])
                unit_num = int(mol.split('_')[2])
                polymer_num = int(mol.split('_')[3])
                confro_num = int(mol.split('_')[-1].split('.')[0])
                energy = float(line.split()[1])
                list1.append([ring_num,unit_num,polymer_num,confro_num,energy])
            except:
                pass
    biglist = []
    ring_type = [55,65,66]
    unit_nums = [1,3,5,7]
    for r in ring_type:
        for u in unit_nums:
            i = 0
            while i <= 100:
                energies = []
                for j in list1:
                    if j[0] == r and j[1] == u and j[2] == i:
                        energies.append(j[4])
                        m = min(energies)
                for j in list1:
                    if j[0] == r and j[1] == u and j[2] == i and j[4] == m:
                        biglist.append([j[0],j[1],j[2],j[3]])
                        #[ring_num, unit_num, polymer_num, confro_num]
                i+=1
    print('Done making list of lowest energy confromers')
    for i in range(len(biglist)):
        stronk = 'mols_{}_{}_{:02d}_*_{}.log'.format(biglist[i][0],biglist[i][1],biglist[i][2],biglist[i][3])
        subprocess.call('cp ../gau_output/{}/{} ../wtuning'.format(stronk[:-4],stronk),shell=True)

def WriteWTuning():
    """
    Write runfile for wtuning job
    """
    with open('wtuning.py','w+') as wt:
        wt.write("""from pymatgen.core.structure import Molecule
from sys import argv
from ocelot.task.wtuning import WtuningJob

fn = argv[1]
name = '_'.join(fn.split('_')[:-1])
pymol = Molecule.from_file(fn)
pymol.perturb(0.05)

job = WtuningJob(func='uLC-wPBE',basis='cc-pVDZ',name=name,nproc=8,mem=30,n_charge=0,n_spin=1,scheme='Jh')
route_parameters = job.route_params
route_parameters.update({'EmpiricalDispersion':'GD3'})
route_parameters.update({"scf":{"maxcyc=999":""}})
job.mol = pymol
job.geo_opt()
print('Done optimizing geometry with uLC-wPBE/cc-pVDZ.')
job.wtuning_cycle(max_cycles=0)
print('Done with wtuning.')""")

def RunJob(job_name,job):
    """
    Write runfile for and submit any SLURM job. job_name and job should each be a string.
    """
    runfile = os.path.join(os.getcwd(),"%s.job" %job_name)
    with open(runfile,"w+") as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH -t 14-00:00:00\n")
        fh.writelines("#SBATCH --job-name=%s\n" % job_name)
        fh.writelines("#SBATCH -N 1\n")
        fh.writelines("#SBATCH -n 8\n")
        fh.writelines("#SBATCH -p SKY32M192_L\n")
        fh.writelines("#SBATCH --account=col_cmri235_uksr\n")
        fh.writelines("#SBATCH --mail-type=ALL\n")
        fh.writelines("#SBATCH --mail-user=rdu230@uky.edu\n")
        fh.writelines("#SBATCH --error=SLURM_JOB_%j.err\n")
        fh.writelines("#SBATCH --output=SLURM_JOB_%j.out\n" )

        fh.writelines("ulimit -u 32768\n")
        fh.writelines("module load ccs/gaussian/g16-A.03/g16-sandybridge\n")
        fh.writelines('echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "\n')
        fh.writelines(". ~/.bashrc\n")
        fh.writelines("source $(conda info --base)/etc/profile.d/conda.sh\n")
        fh.writelines("conda activate myenv\n")

        fh.writelines(job)

    os.system("sbatch %s" %runfile)

def MoveFiles(destination,starts_with="",ends_with=""):
    """
    Moves all files that start/end with a certain string to a new destination
    """
    for i in os.listdir(os.getcwd()):
        if i.startswith(starts_with) and i.endswith(ends_with):
            shutil.move(i,destination+i)

def FuncAllFiles(func,starts_with="",ends_with=""):
    """
    Runs all files that start/end with a certain string through a function
    """
    for i in os.listdir(os.getcwd()):
        if i.startswith(starts_with) or i.endswith(ends_with):
            func(i)

home = os.getcwd()
os.chdir(home+'/gau_output')
FindEnergies()
print('Done finding energies')
FindBestConfromer('table.txt')
print('Done finding best confromer')

os.chdir(home+'/wtuning')
for f in os.listdir(os.getcwd()):
    file_name = '_'.join(f.split('_')[:-1])
    os.mkdir("{}/wtuning/{}/".format(home,file_name))
    MoveFiles("{}/wtuning/{}/".format(home,file_name),starts_with=str(f))
    os.chdir(file_name)
    WriteWTuning()
    RunJob('wtuning{}'.format(f[5:12]),"python wtuning.py "+str(f))
    os.chdir(home+'/wtuning')
print('Done Submitting wtuning')
