import os
import re
import json
import pymatgen
from pymatgen.core import Molecule
from pymatgen.core import structure

def MakeJSON(mol_dir,json_file):
    """
    For a molecule rotation directory, mol_dir, this finds the energy, xyz, and
    degree of rotation and places those into  json_file.
    """
    mol_name = str(mol_dir.split('/')[-1])
    degree = str((mol_name.split('_')[-1])[:-3])
    energy = 0
    mol_xyz = 0
    try:
        for f in mol_dir:
            gen = [x for x in os.listdir(mol_dir) if x.endswith('.log')]
            for mol_file in gen:
                mol_file_path = os.path.join(mol_dir, mol_file)
                with open(mol_file_path) as fn:
                    normal = 0
                    for line in fn:
                        if re.search('Normal termination', line):
                            normal = 1
                        if re.search('SCF Done', line):
                            possible_E = line.split()[4]
                    if normal == 1:
                        energy = possible_E
                mol_pymat = Molecule.from_file(mol_file_path)
                mol_xyz = mol_pymat.to(fmt="xyz")
    except:
        print("Error! Data NOT collected for {}. Energy calculations for dihedral rotations may not have finished!".format(mol_name))
    if energy == 0 or mol_xyz == 0:
        print("Error! Data NOT collected for {}. Energy calculations for dihedral rotations may not have finished!".format(mol_name))
    else:
        print(degree)
        json_data = {
            "molecule_name" : mol_name,
            "degree" : degree,
            "xyz" : mol_xyz,
            "energy (Hartrees)" : energy
        }
        print(json_data)
        with open(json_file,'w+') as fn:
            json.dump(json_data, fn)
        print("Data collected for {}".format(d))


cwd = os.getcwd()

home = os.path.join(cwd,'rotated_dihed')
mols = [m for m in os.listdir(home) if os.path.isdir(os.path.join(home,m))]
for m in mols:
    mpath = os.path.join(home,m)
    degs = [d for d in os.listdir(mpath) if os.path.isdir(os.path.join(mpath,d))]
    for d in degs:
        dpath = os.path.join(home,m,d)
        json_file = "{}/jsons_nnff/{}.json".format(cwd,d)
        if not os.path.isfile(json_file):
            MakeJSON(dpath,json_file)
