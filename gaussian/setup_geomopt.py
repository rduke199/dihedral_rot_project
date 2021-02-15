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

# full_redo_list = [  'mols_55_7_09_COC-COC',   'mols_65_5_25_n-cOC-N',     'mols_65_7_37_cOC-cOC-N',        'mols_66_5_13_c-n-cF-cOC',   'mols_66_7_13_c-n-cF-cOC',
#                     'mols_65_1_36_cOC-cOC-C',    'mols_65_5_37_cOC-cOC-N',  'mols_66_1_07_c-n-c-cOC',        'mols_66_5_14_n-n-cF-cOC',   'mols_66_7_14_n-n-cF-cOC',
#                     'mols_65_1_37_cOC-cOC-N',    'mols_65_7_13_c-cOC-N',   'mols_66_3_16_c-n-cOC-cOC',      'mols_66_5_16_c-n-cOC-cOC',  'mols_66_7_16_c-n-cOC-cOC',
#                     'mols_65_1_39_cOC-cOC-COC',  'mols_65_7_27_n-cOC-COC',  'mols_66_5_07_c-n-c-cOC',        'mols_66_5_17_n-n-cOC-cOC',  'mols_66_7_17_n-n-cOC-cOC',
#                     'mols_65_3_37_cOC-cOC-N',    'mols_65_7_33_cF-cOC-N',   'mols_66_5_07_c-n-c-cOC_opt_0',  'mols_66_7_07_c-n-c-cOC',
#                     'mols_65_3_39_cOC-cOC-COC',  'mols_65_7_36_cOC-cOC-C',  'mols_66_5_08_n-n-c-cOC',        'mols_66_7_08_n-n-c-cOC',
#                 ]
full_redo_list = [  'mols_66_3_00_c-c-c-', 'mols_66_3_02_n-n-c-', 'mols_66_3_08_n-n-c-cOC', 'mols_66_5_00_c-c-c-', 'mols_66_5_02_n-n-c-',
                    'mols_66_5_05_n-n-c-cF', 'mols_66_7_02_n-n-c-', 'mols_66_7_03_c-c-c-cF', 'mols_66_7_05_n-n-c-cF', 'mols_65_3_30_cF-cF-CF',
                    'mols_65_3_32_cF-cOC-C', 'mols_65_3_34_cF-cOC-CF', 'mols_65_3_37_cOC-cOC-N', 'mols_65_3_39_cOC-cOC-COC', 'mols_65_5_01_c-c-N',
                    'mols_65_5_06_c-n-CF', 'mols_65_5_07_c-n-COC', 'mols_65_5_10_c-cF-CF', 'mols_65_5_11_c-cF-COC', 'mols_65_5_18_n-n-CF',
                    'mols_65_5_19_n-n-COC', 'mols_65_5_32_cF-cOC-C', 'mols_65_5_37_cOC-cOC-N', 'mols_65_7_00_c-c-C', 'mols_65_7_01_c-c-N',
                    'mols_65_7_03_c-c-COC', 'mols_65_7_04_c-n-C', 'mols_65_7_05_c-n-N', 'mols_65_7_07_c-n-COC', 'mols_65_7_09_c-cF-N',
                    'mols_65_7_13_c-cOC-N', 'mols_65_7_14_c-cOC-CF', 'mols_65_7_18_n-n-CF', 'mols_65_7_28_cF-cF-C', 'mols_65_7_35_cF-cOC-COC',
                    'mols_65_7_38_cOC-cOC-CF', 'mols_66_1_05_n-n-c-cF', 'mols_66_1_08_n-n-c-cOC', 'mols_65_1_05_c-n-N', 'mols_65_1_12_c-cOC-C',
                    'mols_65_1_32_cF-cOC-C', 'mols_65_1_34_cF-cOC-CF', 'mols_65_1_35_cF-cOC-COC', 'mols_65_3_01_c-c-N', 'mols_65_3_05_c-n-N',
                    'mols_65_3_06_c-n-CF', 'mols_65_3_09_c-cF-N', 'mols_65_3_10_c-cF-CF', 'mols_65_3_11_c-cF-COC', 'mols_65_3_13_c-cOC-N',
                    'mols_65_3_15_c-cOC-COC'
                ]


def SetupFiling(gather_dir,out_dir,angles_fn, redo=False):
    with open(angles_fn) as f:
        angle_dict = json.load(f)
    for mol in angle_dict.keys():
        if redo == False or mol in full_redo_list:
            angle = angle_dict[mol]
            mol_fn = [i for i in os.listdir(gather_dir) if i.startswith(mol)][0]
            gather_mol_path = os.path.join(gather_dir,mol_fn)
            out_mol_path = os.path.join(out_dir,mol_fn)
            if os.path.isdir(gather_mol_path):
                # set up rotated_dihed filing system
                gather_angle_path = os.path.join(gather_mol_path,[i for i in os.listdir(gather_mol_path) if i.endswith("{:.0f}deg".format(angle))][0])
                if not os.path.isdir(out_mol_path): os.mkdir(out_mol_path)
                subprocess.call('rm {}/*deg.log'.format(out_mol_path),shell=True)
                subprocess.call('rm {}/*.gjf'.format(out_mol_path),shell=True)
                try:
                    log_fn = [i for i in os.listdir(gather_angle_path) if i.endswith('.log')][0]
                    copyfile(os.path.join(gather_mol_path,'output.log'), os.path.join(out_mol_path,'output.log'))
                    copyfile(os.path.join(gather_angle_path,log_fn), os.path.join(out_mol_path,log_fn))
                except:
                    print("Error copying files for {}".format(mol))
                    continue

                txt_file = os.path.join(out_dir,'go_folders_to_run.txt')
                with open(txt_file, "a+") as fn:
                    fn.write(out_mol_path+"\n")

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

SetupFiling(gather_dir,dft_dir,angles_fn,redo=True)
mols = [m for m in os.listdir(dft_dir)]
for m in mols:
    mpath = os.path.join(dft_dir,m)
    if os.path.isdir(mpath):
        GetRunFolders(mpath, run_folder_dir, nflag=m[8])
