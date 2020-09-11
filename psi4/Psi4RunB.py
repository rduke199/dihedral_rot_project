import re
import os
import shutil
import subprocess

"""
After initital optimizations have been run on each confromer of all polymers,
this script finds the lowest energy confromer of each polymer. It then runs
ip fitting on the best confromer.
"""

def FindBestConfromer(energy_file='energy_table.txt'):
    """
    Find lowest energy confromer of each polymer from energy file. Moves
    best confromer '_geometry_optimization' files to the ip_fitting directory.
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
    for i in range(len(biglist)):
        final_geom = 'mols_{}_{}_{}_*_{}_geometry_final.xyz'.format(biglist[i][0],biglist[i][1],biglist[i][2],biglist[i][3])
        smile_file = 'mols_{}_{}_{}_*_{}.smi'.format(biglist[i][0],biglist[i][1],biglist[i][2],biglist[i][3])
        subprocess.call('cp ../psi4_output/{} ../ip_fitting/'.format(final_geom),shell=True)
        subprocess.call('cp ../smiles/{} ../rotated_dihed/'.format(smile_file),shell=True)

def WriteIPFit(script_name='ipfitting.py'):
    text="""import psi4
from psi4.driver.frac import ip_fitting
from sys import argv

mol_file = argv[1]
molecule_name = str(mol_file.split('/')[-1])[:-19]
with open(mol_file,'r') as mol:
    mol = psi4.core.Molecule.from_string(mol.read(), dtype='xyz')
    mol.set_molecular_charge(0)
    mol.set_multiplicity(1)

psi4.set_memory('2 GB')
psi4.set_num_threads(2)
psi4.set_output_file(molecule_name + '_ip_fitting.dat', False)
psi4.set_options({'reference': 'uks', 'basis': 'cc-pvdz'})
omega = ip_fitting('LRC-wPBE', 0.1, 2.0,molecule=mol)
with open("omega.txt",'a') as fn:
    fn.write(str(omega))"""

    with open(script_name,'w+') as wt:
        wt.write(text)

def RunJob(job_name,job):
    """
    Write runfile for and submit any psi4 SLURM job. job_name and job should
    each be a string.
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
        fh.writelines('echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "\n')
        fh.writelines(". ~/.bashrc\n")
        fh.writelines("conda deactivate\n")
        fh.writelines("module load ccs/conda/psi4-1.3.2+ecbda83\n")

        fh.writelines(job)

    os.system("sbatch %s" %runfile)

def MoveFiles(destination,starts_with="",ends_with=""):
    """
    Moves all files that start/end with a certain string to a new destination
    """
    for i in os.listdir(os.getcwd()):
        if i.startswith(starts_with) and i.endswith(ends_with):
            shutil.move(i,destination+i)

home = os.getcwd()

os.chdir(home+'/psi4_output')
FindBestConfromer()
print('Done finding best confromers')

os.chdir(home+'/ip_fitting')
for f in os.listdir(os.getcwd()):
    molecule_name = str(f.split('/')[-1])[:-19]
    os.mkdir("{}/{}/".format(os.getcwd(),molecule_name))
    shutil.move(f,"{}/{}/{}".format(os.getcwd(),molecule_name,f))
    os.chdir("{}/{}/".format(os.getcwd(),molecule_name))
    WriteIPFit()
    RunJob('ipfit_{}'.format(f[5:12]),'python ipfitting.py {}'.format(f))
    os.chdir(home+'/ip_fitting')
print("DONE ip fitting jobs submitted!")
