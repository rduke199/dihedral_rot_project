import os
import re
import gc
import json
from pymatgen.core import Molecule
from pymatgen.core import structure
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.core.structure import IMolecule

def write_json(data, filename):
    with open(filename,'w') as f:
        json.dump(data, f, indent=2)

def MakeJSON(mol_dir,dihed_atoms_file,json_file):
    """
    For a molecule rotation directory, mol_dir, this finds the energy, xyz, and
    degree of rotation and places those into  json_file.
    """
    mol_name = str(mol_dir.split('/')[-1])
    try:
        log_fn = [x for x in os.listdir(mol_dir) if x.endswith('.log') and x.startswith('mol') and not x.endswith('deg.log')][0]
        log_path = os.path.join(mol_dir,log_fn)
        mol = GaussianOutput(log_path)
        if mol.properly_terminated:
            num_electrons = mol.electrons[0]
            eigens = list(mol.eigenvalues.values())[0]
            homo = eigens[num_electrons - 1] * 27.2114
            lumo = eigens[num_electrons] * 27.2114
            homo_lumo_gap = lumo - homo
            energy = mol.final_energy * 27.2114

            imol = IMolecule.from_file(log_path)
            with open(dihed_atoms_file,'r') as lines:
                data = lines.readlines()
                dihedral1 = data[0].split(',')[0]
            d1atom1_num = int(dihedral1.split()[0])
            d1atom2_num = int(dihedral1.split()[1])
            d1atom3_num = int(dihedral1.split()[2])
            d1atom4_num = int(dihedral1.split()[3])
            di_angle = imol.get_dihedral(d1atom1_num, d1atom2_num, d1atom3_num, d1atom4_num)

            normal = 1
        else:
            normal = 0
            print("Error. GeomOpt calculation did not properly terminate for {}".format(mol_name))
    except:
        normal = 0
        print("Error. Data NOT collected for {}. Energy calculations for dihedral rotations may not have finished!".format(mol_name))
    if normal == 1:
        json_data = {
        "molecule_name" : mol_name,
        "dihedral_angle" : di_angle,
        "energy" : energy,
        "homo" : homo,
        "lumo" : lumo,
        "homo_lumo_gap" : homo_lumo_gap
        }
        write_json(json_data,json_file)
        print("Data collected for {}".format(mol_name))
    else:
        pass

def WriteMaster(json_dir,master_json_file):
    data = {}
    data['molecules'] = []
    for f in os.listdir(json_dir):
        fpath = os.path.join(json_dir,f)
        with open(fpath) as fn:
            mol_info = json.load(fn)
        data['molecules'].append(mol_info)
    write_json(data, master_json_file)

cwd = os.getcwd()
gather_home = os.path.join(cwd,'geomopt_dft/')
dihedral_home = os.path.join(cwd, 'diheds/')
out_home = os.path.join(cwd,'geomopt_jsons/')
master_json_file = os.path.join(out_home,"master_geomopt.json")

for m in os.listdir(gather_home):
    mpath = os.path.join(gather_home,m)
    if os.path.isdir(mpath):
        json_file = "{}/{}.json".format(out_home,m)
        if os.path.isfile(json_file) == False:
            dihed_atoms_file = [x for x in os.listdir(dihedral_home) if x.startswith(m[:12])][0]
            dihpath = os.path.join(dihedral_home,dihed_atoms_file)
            MakeJSON(mpath,dihpath,json_file)
WriteMaster(out_home,master_json_file)
