import os
import numpy as np

import ocelot
import pymatgen
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput
from rdkit.Chem.rdmolfiles import MolToXYZBlock
from rdkit.Chem.rdMolTransforms import SetDihedralDeg
from ocelot.routines.conformerparser import pmgmol_to_rdmol

def MakesGauInputs(mol_file,dihed_atoms_file,w,num_conform=19):
    """
    Creates Gaussian input files in its own folder for mol_file structues that
    freezes the central angle at iterations of 180/num_conform degrees. This
    function requires Gaussian output file of molecule in question and file
    with dihedral atoms to be in the current working directory.
    """

    molecule_name = mol_file.split('/')[-1][:-12]
    molecule = pmgmol_to_rdmol(Molecule.from_file(mol_file))[0]

    with open(dihed_atoms_file,'r') as fn:
        data = fn.readlines()
        dihedral1 = data[0].split(',')[0]
        dihedral2 = data[0].split(',')[1]

    d1atom1_num = int(dihedral1.split()[0])
    d1atom2_num = int(dihedral1.split()[1])
    d1atom3_num = int(dihedral1.split()[2])
    d1atom4_num = int(dihedral1.split()[3])

    phi = np.linspace(start=0, stop=180, num=num_conform)
    for P in phi:
        dir_name = "{}/{}_{:.0f}deg/".format(os.getcwd(),molecule_name,P)
        file_name = "{}/{}_{:.0f}deg".format(dir_name,molecule_name,P)
        os.mkdir(dir_name)

        SetDihedralDeg(molecule.GetConformer(),d1atom1_num, d1atom2_num, d1atom3_num, d1atom4_num, P)
        with open(dir_name+"temp.xyz",'w+') as fn:
            fn.write(str(MolToXYZBlock(molecule)))
        mol = Molecule.from_file(dir_name+'temp.xyz')
        os.remove(dir_name+'temp.xyz')
        gau = GaussianInput(mol=mol,charge=0,spin_multiplicity=1,functional='uLC-wPBE',basis_set='cc-pVDZ',route_parameters={"iop(3/107={}, 3/108={})".format(w,w):"","opt":"modredundant"},link0_parameters={'%mem':'5GB','%chk':'{}.chk'.format('{}.chk'.format(file_name.split('/')[-1]))})
        gjf_file = gau.write_file(dir_name+'temp.gjf')
        (dir_name+'temp.gjf').close()

        with open(dir_name+'temp.gjf') as temp:
            lines = temp.readlines()
        os.remove(dir_name+'temp.gjf')
        with open(file_name+'.gjf','w') as gjf:
            gjf.writelines([item for item in lines[:-2]])
            gjf.write("D * {} {} * F\n\n".format(d1atom2_num+1,d1atom3_num+1))

    print("Torsion Gaussian input files for {} finsihed!".format(molecule_name))

#This part just goes into the correct file and runs MakesGauInputs.
home = os.getcwd()
os.chdir(home)
for i in os.listdir(os.getcwd()):
    if i.startswith("mols_"):
        os.chdir(i)
        for d in os.listdir(os.getcwd()):
            if os.path.isdir(d):
                break
            else:
                print(i)
                try:
                    with open('output.log','r') as fn:
                        w_data = fn.readlines()[-2].split()[1]
                        w = "0{}".format(w_data.split('.')[1])
                    dihed_atoms_file = i+'_dihed.txt'
                    MakesGauInputs('{}_opt_0.log'.format(i),dihed_atoms_file,w)
                    print("Gaussian input files made for {}!".format(i))
                    break
                except:
                    pass
        os.chdir(home)
