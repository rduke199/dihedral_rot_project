import os
import re
import json


def write_json(data, filename):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)


def find_omega(omega_file):
    if os.path.isfile(omega_file):
        with open(omega_file, 'r') as fn:
            lines = fn.readlines()
            for line in reversed(lines):
                if re.search('new:', line):
                    omega_data = line.split()[1]
                    omega = "0{}".format(omega_data.split('.')[1])
                    return omega


def main():
    cwd = os.getcwd()
    gather_home = os.path.join(cwd, 'rotated_dihed')
    json_dir = os.path.join(cwd, 'jsons_rotdiheds')
    master_omega_file = os.path.join(json_dir, "master_omega_rd.json")

    omega_data = {}
    mols = [m for m in os.listdir(gather_home) if os.path.isdir(os.path.join(gather_home, m))]
    for mol in mols:
        mpath = os.path.join(gather_home, mol)
        mol_name = mpath.split('/')[-1]
        degs = [d for d in os.listdir(mpath) if os.path.isdir(os.path.join(mpath, d))]
        for deg in degs:
            gjf_file = os.path.join(mpath, deg, deg + ".gjf")
            try:
                with open(gjf_file, 'r') as fn:
                    omega_line = fn.readlines()[3]
                    omega = omega_line.split('3/108=')[1].split(')')[0]
                    omega_data[mol_name] = omega
                    print(omega)
                    break
            except:
                print('No omega gathered for {}'.format(mol_name))

        omega_file = os.path.join(mpath, 'output.log')
        omega = find_omega(omega_file)
        omega_data[mol_name] = omega

    write_json(omega_data, master_omega_file)


if __name__ == "__main__":
    main()
