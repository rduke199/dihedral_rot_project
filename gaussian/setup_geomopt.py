import os
import re
import json
import subprocess
from shutil import copyfile
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput
from pymatgen.io.gaussian import GaussianOutput

def MakesGauInputs(fn,w):
    """
    Convert an individually imputted xyz file to a gjf Gaussian input file
    """
    mol = Molecule.from_file(fn)
    mol_name = fn.split('/')[-1][:-12]
    gau = GaussianInput(mol=mol,charge=0,spin_multiplicity=1,functional='uLC-wPBE',basis_set='cc-pVDZ',  route_parameters={"iop(3/107={}, 3/108={})".format(w,w):"","opt":"modredundant"},   link0_parameters={'%mem':'5GB','%chk':'{}.chk'.format(mol_name)})
    gjf_file = gau.write_file('{}.gjf'.format(fn.split('.')[0]))
    return gjf_file

def SetupFiling(gather_dir,out_dir,angles_fn):
    with open(angles_fn) as f:
        angle_dict = json.load(f)
    for mol in angle_dict.keys():
        angle = angle_dict[mol]
        mol_fn = [i for i in os.listdir(gather_dir) if i.startswith(mol)][0]
        gather_mol_path = os.path.join(gather_dir,mol_fn)
        out_mol_path = os.path.join(out_dir,mol_fn)
        if os.path.isdir(gather_mol_path) and not os.path.isdir(out_mol_path):
            # set up rotated_dihed filing system
            gather_angle_path = os.path.join(gather_mol_path,[i for i in os.listdir(gather_mol_path) if i.endswith("{:.0f}deg".format(angle))][0])
            os.mkdir(out_mol_path)
            try:
                log_fn = [i for i in os.listdir(gather_angle_path) if i.endswith('.log')][0]
                copyfile(os.path.join(gather_mol_path,'output.log'), os.path.join(out_mol_path,'output.log'))
                copyfile(os.path.join(gather_angle_path,log_fn), os.path.join(out_mol_path,log_fn))
            except:
                print("Error copying files for {}".format(mol))
                continue

            # create gaussian input files for all degree rotations and run optimizations
            os.chdir(out_mol_path)
            try:
                with open('output.log','r') as fn:
                    w_data = fn.readlines()[-2].split()[1]
                    w = "0{}".format(w_data.split('.')[1])
                MakesGauInputs(log_fn,w)
                print("Gaussian input files made for {}!".format(mol))
            except:
                print("Error. Gaussian input not made for {}".format(mol))
        else:
            continue

def GetRunFolders(molecule_dir, out_dir, nflag='_'):
    mol_name = molecule_dir.split('/')[-1].split('.')[0]
    log_files = [x for x in os.listdir(molecule_dir) if x.endswith('.log') and x.startswith('mol') and not x.endswith('deg.log')]
    gjf_files = [x for x in os.listdir(molecule_dir) if x.endswith('.gjf')]
    if len(log_files) == 1:
        with open(os.path.join(molecule_dir,log_files[0])) as fn:
            log_fn = fn.readlines()
            last_line = log_fn[-1]
            if re.search('Normal termination', last_line):
                normal = 0
                for line in log_fn:
                    if re.search(' Optimized Parameters',line):
                        normal = 1
                        break
                    else:
                        continue
                if normal == 0:
                    print("Error. {} job had a normal termination but did not have optimized parameters.".format(mol_name))
            else:
                normal = 0
                print("Error. {} job did not have a normal termination.".format(mol_name))
    elif len(log_files) == 0:
        normal = 0
        print("Error. {} has no log file(s).".format(mol_name, len(log_files)))
    elif len(log_files) > 1:
        normal = 2
        print("Error. {} has {} log file(s).".format(mol_name, len(log_files)))
    if len(gjf_files) == 0:
        normal = 2
        print("Error. {} does not contain a gjf file.".format(mol_name))
    if normal == 0:
        txt_file = os.path.join(out_dir,'go_folders{}_to_run.txt'.format(nflag))
        with open(txt_file, "a+") as fn:
            fn.write(molecule_dir+"\n")

home = os.getcwd()
gather_dir = os.path.join(home,'rotated_dihed/')
dft_dir = os.path.join(home,'geomopt_dft/')
angles_fn = os.path.join(home, 'lowest_e_angles.txt')
run_folder_dir = os.path.join(home,'test/')

# SetupFiling(gather_dir,dft_dir,angles_fn)
mols = [m for m in os.listdir(dft_dir)]
for m in mols:
    mpath = os.path.join(dft_dir,m)
    if os.path.isdir(mpath):
        GetRunFolders(mpath, run_folder_dir, nflag=m[8])
