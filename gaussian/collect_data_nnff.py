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
            gen = [x for x in os.listdir(os.getcwd()) if x.endswith('.log')]
            for mol_file in gen:
                with open(mol_file) as fn:
                    for line in fn:
                        if re.search(r'SCF Done',line):
                            energy = line.split()[4]
                mol_pymat = Molecule.from_file(mol_file)
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
            "energy" : energy
        }
        print(json_data)
        with open(json_file,'w+') as fn:
            json.dump(json_data, fn)
        print("Data collected for {}".format(d))


cwd = os.getcwd()

home = os.path.join(cwd,'rotated_dihed')
mols = os.listdir(home)
os.chdir(home)
for m in mols:
    if os.path.isdir(m):
        mpath = os.path.join(home,m)
        os.chdir(mpath)
        degs = os.listdir(mpath)
        for d in degs:
            if os.path.isdir(d):
                dpath = os.path.join(home,m,d)
                os.chdir(dpath)
                json_file = "{}/jsons/{}.json".format(cwd,d)
                MakeJSON(dpath,json_file)
            os.chdir(mpath)
    os.chdir(home)
