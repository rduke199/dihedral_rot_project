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


def check_normal_optimization(log_file):
    with open(log_file) as fn:
        last_line = fn[-1]
        if re.search('Normal termination', last_line):
            for line in fn:
                if re.search(' Optimized Parameters', line):
                    return True
                else:
                    raise UserWarning("Calculation terminated normally but did not sucessfully optimze parameters for {}".format(log_file))
        else:
            raise UserWarning("Termination error for {}".format(log_file))


def make_json(mol_dir, dihed_atoms_file, json_file):
    """
    For a molecule rotation directory, mol_dir, this finds the energy, xyz, and
    degree of rotation and places those into  json_file.
    """
    mol_name = str(mol_dir.split('/')[-1])
    try:
        log_fn = [x for x in os.listdir(dpath) if x.endswith('deg.log')][0]
        log_path = os.path.join(dpath, log_fn)
    except IndexError:
        print("Error. No log file for {} at {} degrees.".format(mol_name, degree))
        return None
    if check_normal_optimization(log_path):
        mol = GaussianOutput(log_path)
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
        normal = 0
        print("Error. GeomOpt calculation did not properly terminate for {}".format(mol_name))


def write_master_json(json_dir, master_json_file):
    data = {}
    for f in os.listdir(json_dir):
        fpath = os.path.join(json_dir, f)
        with open(fpath) as fn:
            mol_info = json.load(fn)
        mol_name = mol_info["molecule_name"]
        dihed_angle = mol_info["dihedral_angle"]
        energy = mol_info["energy"]
        data[mol_name] = (dihed_angle,energy)
    write_json(data, master_json_file)


cwd = os.getcwd()
gather_home = os.path.join(cwd,'geomopt_dft/')
dihedral_home = os.path.join(cwd, 'diheds/')
out_home = os.path.join(cwd,'geomopt_jsons/')
master_json_file =os.path.join(out_home,"master_geomopt.json")

# for m in os.listdir(gather_home):
#     mpath = os.path.join(gather_home,m)
#     if os.path.isdir(mpath):
#         json_file = "{}/{}.json".format(out_home,m)
#         if os.path.isfile(json_file) == False:
#             dihed_atoms_file = [x for x in os.listdir(dihedral_home) if x.startswith(m[:12])][0]
#             dihpath = os.path.join(dihedral_home,dihed_atoms_file)
#             MakeJSON(mpath,dihpath,json_file)
write_master_json(out_home,master_json_file)
