import os
import json
import subprocess

home = os.getcwd()
xyz_dir = os.path.join(home, 'xyz/')
w_dir = os.path.join(home, 'wtuning/')

w_tuning_file = os.path.join(home, 'wtuning.py')
run_dirs_file = os.path.join(w_dir, 'folders_to_run.txt')
energy_data = {}

conformers = [c for c in os.listdir(xyz_dir) if c.endswith('.xyz')]
# collect energies
for conformer in conformers:
    try:
        file_path = os.path.join(xyz_dir, conformer)
        with open(file_path, 'r') as fn:
            data = fn.readlines()
        energy = data[1].split(' ')[2]
        mol = '_'.join(conformer.split('_')[:-1])
        end = conformer.split('_')[-1]
        conformer_num = end.split('.')[0]
        if mol in energy_data.keys():
            energy_data[mol][conformer_num] = energy
        else:
            energy_data[mol] = {}
    except:
        print('Energy not collected for {}'.format(mol))

# get min energy conformation for each molecule
for mol, mol_data in energy_data.items():
    low_e_conf = min(mol_data, key=mol_data.get)
    energy_data[mol]['lowest_energy'] = {low_e_conf: mol_data[low_e_conf]}
    low_e_xyz_file = "{}/{}_{}.xyz".format(xyz_dir, mol, low_e_conf)
    mol_w_dir = os.path.join(w_dir, mol)
    if not os.path.isdir(mol_w_dir):
        os.mkdir(mol_w_dir)
    subprocess.call('cp {} {}/{}.xyz'.format(low_e_xyz_file, mol_w_dir, mol), shell=True)
    subprocess.call('cp {} {}'.format(w_tuning_file, mol_w_dir), shell=True)
    # write w tuning directory to run folders file.
    with open(run_dirs_file, 'a+') as fn:
        fn.write(mol_w_dir + '\n')


with open('energy_data.json', 'w') as fn:
    json.dump(energy_data, fn, indent=2)
