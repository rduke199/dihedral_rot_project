import os
import shutil
import subprocess

"""
Runs psi4 energy calculations on each polymer (a smiles file) at varying
conformations where the central angle is rotated 180/num_conform degrees.
It uses the omega values from the ip fitting.
"""

def WritePsi4DihdsE(script_name='psi4_dihdsE.py'):
    text ="""import os
import numpy as np
import subprocess
import shutil

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolTransforms import SetDihedralDeg
from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.Chem.rdmolfiles import MolToXYZBlock

import psi4

def DihedEnergies(mol_file,dihed_atoms_file,w,num_conform=19):
    with open(dihed_atoms_file,'r') as dt:
        data = dt.readlines()
        d1 = data[0].split(',')[0]
        d2 = data[0].split(',')[1]

    #Note these atom numbers begin with atom 0. Rdkit begins with atom 0, but psi4 begins with atom 1.
    d1atom1_num = int(d1.split()[0])
    d1atom2_num = int(d1.split()[1])
    d1atom3_num = int(d1.split()[2])
    d1atom4_num = int(d1.split()[3])

    d2atom1_num = int(d2.split()[0])
    d2atom2_num = int(d2.split()[1])
    d2atom3_num = int(d2.split()[2])
    d2atom4_num = int(d2.split()[3])

    phi = np.linspace(start=0, stop=180, num=num_conform)
    psi4.set_memory('2 GB')

    for P in phi:
        molecule_name = mol_file.split('/')[-1][:-12]
        with open(mol_file,'r') as fn:
            rdmol = MolFromSmiles(str(fn.read()))

        dir_name = "{}/{}_{:.0f}deg/".format(os.getcwd(),molecule_name,P)
        os.mkdir(dir_name)
        file_name = "{}/{}_{:.0f}deg".format(dir_name,molecule_name,P)

        AllChem.EmbedMolecule(rdmol)
        SetDihedralDeg(rdmol.GetConformer(),d1atom1_num, d1atom2_num, d1atom3_num, d1atom4_num, P)
        with open(file_name+'.xyz','w+') as fout:
            fout.write(str(MolToXYZBlock(rdmol)))

        with open(file_name+'.xyz','r') as mol:
            molecule = psi4.core.Molecule.from_string(mol.read(), dtype='xyz')
            molecule.set_molecular_charge(0)
            molecule.set_multiplicity(1)

        #Remember that psi4 begins with atom 1.
        frozen_dihedral1 = "{} {} {} {}  {}".format(d1atom1_num+1, d1atom2_num+1, d1atom3_num+1, d1atom4_num+1, P)
        frozen_dihedral2 = "{} {} {} {}  {}".format(d2atom1_num+1, d2atom2_num+1, d2atom3_num+1, d2atom4_num+1, P)
        print("\nThe value of the frozen_dihedrals is: ", P)
        frozen_dihedral_total = frozen_dihedral1 + " " + frozen_dihedral2
        psi4.set_module_options('alpha',{'DFT_OMEGA':w})
        psi4.set_module_options('optking', {'ensure_bt_convergence': True})
        psi4.set_module_options('optking', {'fixed_dihedral': frozen_dihedral_total})  #set the fixed dihedral
        psi4.set_output_file(file_name+'_geometry_optimization.dat', False)
        psi4.set_options({'reference': 'uks', 'basis': 'cc-pvdz'})
        final_energy = psi4.optimize('LRC-wPBE', molecule=molecule)

        with open('energies_{}.txt'.format(molecule_name),'a') as fn:
            fn.write(str(P) + "\t" + str(final_energy) + "\n")

    print("Conformers of dihedral rodations for {} finsihed!".format(molecule_name))

mol_dir = os.getcwd().split('/')[-1]
dihed_atoms_file = '_'.join(str(mol_dir).split('_')[:-1])+'_dihed.txt'
mol_smi_file = '_'.join(str(mol_dir).split('_')[:-1])+'.smi'
with open('omega.txt','r') as fn:
    w = "0{}".format(fn.readlines()[0].split('.')[1])
DihedEnergies(mol_smi_file,dihed_atoms_file,w)
"""
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
ip_fitting = home+'/ip_fitting'
for mol_dir in os.listdir(ip_fitting):
    path = "{}/{}/".format(ip_fitting,mol_dir)
    if os.path.isdir(path)==True:
        os.chdir(path)
        molecule_name = str(mol_dir)[:-2]
        os.mkdir("{}/rotated_dihed/{}/".format(home,mol_dir))
        subprocess.call('cp *geometry_final.xyz omega.txt ../../rotated_dihed/{}'.format(mol_dir),shell=True)
        os.chdir(home+"/rotated_dihed/")
        dihed_atoms_file = molecule_name+'_dihed.txt'
        mol_smi_file = molecule_name+'.smi'
        MoveFiles("{}/rotated_dihed/{}/".format(home,mol_dir),starts_with=dihed_atoms_file)
        MoveFiles("{}/rotated_dihed/{}/".format(home,mol_dir),starts_with=mol_smi_file)
        os.chdir("{}/rotated_dihed/{}/".format(home,mol_dir))
        WritePsi4DihdsE()
        RunJob('dihd{}.job'.format(mol_dir[5:12]),'python psi4_dihdsE.py')
    os.chdir(ip_fitting)
print("Done finding torsion energies")
