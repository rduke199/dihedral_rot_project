import os
import re
import json
import warnings
from pymatgen.core.structure import IMolecule
from pymatgen.io.gaussian import GaussianOutput
from ocelot.routines.conformerparser import pmgmol_to_rdmol


def write_json(data, filename):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)


def runtime_from_log(logfile):
    """
    Collects runtime in core hours from a logfile
    """
    time_patt = re.compile(r"\d+\.d+|\d+")
    time_data = []
    with open(logfile, "r") as f:
        line = f.readline()
        while line != "":
            if re.match("Job cpu time", line.strip()):
                time_data.extend(time_patt.findall(line))
            line = f.readline()
    if time_data:
        time_data = [float(time) for time in time_data]
        runtime = (time_data[0] * 86400 + time_data[1] * 3600 + time_data[2] * 60 + time_data[3]) / 3600
    else:
        runtime = 0
    return round(runtime, 3)


def check_normal_optimization(log_file):
    with open(log_file) as fn:
        data = fn.readlines()
    last_line = data[-1]
    if not re.search('Normal termination', last_line):
        warnings.warn("Termination error for {}".format(log_file))
    for line in data:
        if re.search(' Optimized Parameters', line):
            return True
    warnings.warn("Calculation terminated normally but did not sucessfully optimze parameters for {}".format(log_file))
    return False


def make_json(mol_dir, json_file, omega_file=None):
    """
    For a molecule directory, mol_dir, this finds the molecule_name, smiles, and tuned omega
    for the polymer. It then finds the optimized structures, energies, homos, lumos, and
    runtimes for each degree of rotation and places those into a json_file.
    """
    mol_name = mol_dir.split('/')[-1]
    xyz_dict, energy_dict, homo_dict, lumo_dict, runtime_dict = {}, {}, {}, {}, {}
    mol_smiles, omega = None, None
    if os.path.isfile(omega_file):
        with open(omega_file, 'r') as fn:
            try:
                omega_data = fn.readlines()[-1].split()[1]
                omega = "0{}".format(omega_data.split('.')[1])
            except IndexError:
                print('Error. wtuning for {} may not have finished'.format(mol_name))
                return None
    deg_folders = [deg for deg in os.listdir(mol_dir) if os.path.isdir(os.path.join(mol_dir, deg))]
    for deg in deg_folders:
        dpath = os.path.join(mol_dir, deg)
        if os.path.isdir(dpath):
            degree = str((deg.split('_')[-1])[:-3])
            try:
                log_fn = [x for x in os.listdir(dpath) if x.endswith('deg.log')][0]
                log_path = os.path.join(dpath, log_fn)
            except IndexError:
                print("Error. No log file for {} at {} degrees.".format(mol_name, degree))
                continue
            if check_normal_optimization(log_path):
                runtime = runtime_from_log(log_path)
                mol = IMolecule.from_file(log_path)
                mol_xyz = mol.to(fmt="xyz")

                out_mol = GaussianOutput(log_path)
                num_electrons = out_mol.electrons[0]
                eigens = list(out_mol.eigenvalues.values())[0]
                homo = eigens[num_electrons - 1] * 27.2114
                lumo = eigens[num_electrons] * 27.2114
                energy = out_mol.final_energy * 27.2114
                if mol_smiles is None:
                    mol_smiles = pmgmol_to_rdmol(mol)[1]
            else:
                print("Error. Optimization did not properly terminate for {} at {} degrees.".format(mol_name, degree))
                continue
            xyz_dict[degree] = mol_xyz
            energy_dict[degree] = energy
            homo_dict[degree] = homo
            lumo_dict[degree] = lumo
            runtime_dict[degree] = runtime
    json_data = {
        # Polymer data
        "molecule_name": mol_name,
        "smiles": mol_smiles,
        "tuned_omega": omega,
        # Rotation dictionaries
        "structures": xyz_dict,
        "energies": energy_dict,
        "homos": homo_dict,
        "lumos": lumo_dict,
        "runtimes": runtime_dict
    }

    write_json(json_data, json_file)
    print("Data collected for {}".format(mol_name))


def write_master_json(json_dir, master_json_file, prop=None):
    data = {}
    json_files = [f for f in os.listdir(json_dir) if f.startswith("mol")]
    for f in json_files:
        fpath = os.path.join(json_dir, f)
        with open(fpath) as fn:
            mol_info = json.load(fn)
        mol_name = mol_info["molecule_name"]
        if prop is None:
            data[mol_name] = mol_info
        else:
            prop_info = mol_info[prop]
            data[mol_name] = prop_info
    write_json(data, master_json_file)


def main():
    cwd = os.getcwd()
    gather_home = os.path.join(cwd, 'rotated_dihed')
    json_dir = os.path.join(cwd, 'jsons_rotdiheds')

    for mol in os.listdir(gather_home):
        mpath = os.path.join(gather_home, mol)
        json_file = "{}/{}.json".format(json_dir, mol)
        if os.path.isdir(mpath) and not os.path.isfile(json_file):
            omega_file = os.path.join(mpath, 'output.log')
            # try:
            make_json(mpath, json_file, omega_file=omega_file)
            # except:
            #     print("Error. Data not collected for {}".format(mol))

    master_energy_file = os.path.join(json_dir, "master_energy.json")
    write_master_json(json_dir, master_energy_file, prop="energies")

    master_homos_file = os.path.join(json_dir, "master_homos.json")
    write_master_json(json_dir, master_homos_file, prop="homos")

    master_lumos_file = os.path.join(json_dir, "master_lumos.json")
    write_master_json(json_dir, master_lumos_file, prop="lumos")

    master_smiles_file = os.path.join(json_dir, "master_smiles.json")
    write_master_json(json_dir, master_smiles_file, prop="smiles")

    master_omega_file = os.path.join(json_dir, "master_omega.json")
    write_master_json(json_dir, master_omega_file, prop="tuned_omega")


if __name__ == "__main__":
    main()
