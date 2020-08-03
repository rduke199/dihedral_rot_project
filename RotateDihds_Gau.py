import os
import re
import numpy as np
import subprocess
import shutil

import pymatgen
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput

import ocelot
from ocelot.routines.conformerparser import pmgmol_to_rdmol
from ocelot.task.wtuning import WtuningJob
from ocelot.task.confgen import ConfGen

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolTransforms import SetDihedralDeg
from rdkit.Chem.rdmolfiles import MolToXYZBlock


def MakesGauInputs(mol_file,dihed_atoms_file,w,num_conform=19):
    """
    Creates Gaussian input files in its own folder for mol_file structues that
    freezes the central angle at iterations of 180/num_conform degrees. This
    function requires Gaussian output file of molecule in question and file
    with dihedral atoms to be in the current working directory.
    """

    molecule_name = mol_file.split('/')[-1][:-12]
    mol_pymat = Molecule.from_file(mol_file)
    molecule = pmgmol_to_rdmol(mol_pymat)[0]

    data = open(dihed_atoms_file,'r').readlines()
    dihedral1 = data[0].split(',')[0]
    dihedral2 = data[0].split(',')[1]

    d1atom1_num = int(dihedral1.split()[0])
    d1atom2_num = int(dihedral1.split()[1])
    d1atom3_num = int(dihedral1.split()[2])
    d1atom4_num = int(dihedral1.split()[3])

    phi = np.linspace(start=0, stop=180, num=num_conform)
    for P in phi:
        dir_name = "{}/{}_{:.0f}deg/".format(os.getcwd(),molecule_name,P)
        os.mkdir(dir_name)
        file_name = "{}/{}_{:.0f}deg".format(dir_name,molecule_name,P)

        SetDihedralDeg(molecule.GetConformer(),d1atom1_num, d1atom2_num, d1atom3_num, d1atom4_num, P)
        with open(dir_name+"temp.xyz",'w+') as fout:
            fout.write(str(MolToXYZBlock(molecule)))
        mol = Molecule.from_file(dir_name+'temp.xyz')
        os.remove(dir_name+'temp.xyz')
        gau = GaussianInput(mol=mol,charge=0,spin_multiplicity=1,functional='uLC-wPBE',basis_set='cc-pVDZ',route_parameters={"iop(3/107={}, 3/108={})".format(w,w):"","opt":"modredundant"},link0_parameters={'%mem':'5GB','%nproc':'1','%chk':'{}.chk'.format('{}.chk'.format(file_name.split('/')[-1]))})
        gjf_file = gau.write_file(dir_name+'temp.gjf')

        with open(dir_name+'temp.gjf') as temp:
            lines = temp.readlines()
        os.remove(dir_name+'temp.gjf')
        with open(file_name+'.gjf','w') as gjf:
            gjf.writelines([item for item in lines[:-2]])
            gjf.write("D * {} {} * F\n\n".format(d1atom2_num+1,d1atom3_num+1))

    print("Torsion Gaussian input files for {} finsihed!".format(molecule_name))

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
        fh.writelines("#SBATCH -n 32\n")
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
        if i.startswith(starts_with) and i.endswith(ends_with):
            func(i)

home = os.getcwd()

os.chdir(home+"/wtuning")
for i in os.listdir(os.getcwd()):
    i_path = os.path.join(home,"wtuning",i)
    if os.path.isdir(i_path):
        os.chdir(i_path)
        os.mkdir("{}/rotated_dihed/{}/".format(home,i))
        try:
            subprocess.call('cp *_opt_0.log output.log ../../rotated_dihed/{}'.format(i),shell=True)
            os.chdir(home+"/rotated_dihed/")
            MoveFiles("{}/rotated_dihed/{}/".format(home,i),starts_with=str(i[:-3]),ends_with='_dihed.txt')
            normal = 1
        except:
            print("opuput.log or opt_0.log file not found for {}. Omega tuning optimization may not have completed!".format(i))
            normal = 0

        if normal == 1:
            os.chdir("{}/rotated_dihed/{}/".format(home,i))
            with open('output.log','r') as fn:
                w_data = fn.readlines()[-2].split()[1]
                w = "0{}".format(w_data.split('.')[1])
            dihed_atoms_file = i+'_dihed.txt'
            print(dihed_atoms_file)
            MakesGauInputs('{}_opt_0.log'.format(i),dihed_atoms_file,w)
            print("Gaussian input files made for {}!".format(i))
        else:
            continue
        os.chdir(home+"/wtuning")

rd = os.path.join(home,'rotated_dihed')
os.chdir(rd)
mols = os.listdir(rd)
for m in mols:
    mpath = os.path.join(rd,m)
    if os.path.isdir(mpath):
    	os.chdir(mpath)
    	degs = os.listdir(mpath)
    	for d in degs:
            dpath = os.path.join(rd,m,d)
            if os.path.isdir(dpath):
                os.chdir(dpath)
                RunJob('gau{}'.format(d[5:]),"g16 {}.gjf".format(d))
                print("Gaussian jobs sumbitted for {}!".format(d))
            else:
                pass
            os.chdir(mpath)
    else:
        pass
    os.chdir(rd)
