import os
import re
import json
import warnings
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.core.structure import IMolecule


def write_json(data, filename):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)


def check_normal_optimization(log_file):
    with open(log_file) as fn:
        data = fn.readlines()
    last_line = data[-1]
    if not re.search('Normal termination', last_line):
        warnings.warn("Termination error for {}".format(log_file))
    for line in data:
        if re.search(' Optimized Parameters', line):
            return True
    warnings.warn("Calculation terminated normally but did not successfully optimize parameters for "
                  "{}".format(log_file))
    return False


def make_json(mol_dir, dihed_atoms_file, json_file):
    """
    For a molecule rotation directory, mol_dir, this finds the energy, xyz, and
    degree of rotation and places those into  json_file.
    """
    mol_name = str(mol_dir.split('/')[-1])
    try:
        log_fn = [x for x in os.listdir(mol_dir) if x.endswith('deg.log')][0]
        log_path = os.path.join(mol_dir, log_fn)
    except IndexError:
        print("Error. No log file for {}.".format(mol_name))
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
        with open(dihed_atoms_file, 'r') as lines:
            data = lines.readlines()
            dihedral1 = data[0].split(',')[0]
        d1atom1_num = int(dihedral1.split()[0])
        d1atom2_num = int(dihedral1.split()[1])
        d1atom3_num = int(dihedral1.split()[2])
        d1atom4_num = int(dihedral1.split()[3])
        di_angle = imol.get_dihedral(d1atom1_num, d1atom2_num, d1atom3_num, d1atom4_num)

        json_data = {
            "molecule_name": mol_name,
            "dihedral_angle": di_angle,
            "energy": energy,
            "homo": homo,
            "lumo": lumo,
            "homo_lumo_gap": homo_lumo_gap
        }
        write_json(json_data, json_file)
        print("Data collected for {}".format(mol_name))
    else:
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
        data[mol_name] = (dihed_angle, energy)
    write_json(data, master_json_file)


cwd = os.getcwd()
gather_home = os.path.join(cwd, 'geomopt_dft/')
dihedral_home = os.path.join(cwd, 'diheds/')
out_home = os.path.join(cwd, 'jsons_geomopt/')
master_json_fn = os.path.join(out_home, "master_geomopt.json")

for m in os.listdir(gather_home):
    mpath = os.path.join(gather_home, m)
    if os.path.isdir(mpath):
        json_fn = "{}/{}.json".format(out_home, m)
        if not os.path.isfile(json_fn):
            dihed_atoms_fn = [x for x in os.listdir(dihedral_home) if x.startswith(m[:12])][0]
            dihpath = os.path.join(dihedral_home, dihed_atoms_fn)
            make_json(mpath, dihpath, json_fn)
write_master_json(out_home, master_json_fn)
