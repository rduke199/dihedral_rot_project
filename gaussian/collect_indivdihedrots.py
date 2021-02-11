import os
import re
from pymatgen.io.gaussian import GaussianOutput

def FindTorsionEnergies(rotdihed_dir, out_dir):
	"""
	Iterates through rotdihed_dir and finds energies for each torsion angle
	for the molecule. It then creates a energy cvs file	for the molecule.
	"""
	mols = [m for m in os.listdir(rotdihed_dir) if os.path.isdir(os.path.join(rotdihed_dir,m))]
	for m in mols:
		mol_name = m.split('.')[0]
		mpath = os.path.join(rotdihed_dir,m)
		log_file = [x for x in os.listdir(mpath) if x.endswith('.log')][0]
		log_path = os.path.join(mpath,log_file)
		mol = GaussianOutput(log_path)
		if mol.properly_terminated:
			normal = 0
			with open(log_path) as fn:
				log = fn.readlines()
				for line in log:
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
		if normal == 1:
			energy = mol.final_energy
			with open('{}/energies.cvs'.format(out_dir),'a') as fout:
				fout.write(mol_name + ',' + str(energy) + '\n')
				print("Dihedral rotations energies collected for {}!".format(mol_name))

cwd = os.getcwd()
dft_dir = os.path.join(cwd,'dft/')
if not os.path.isdir(energy_dir): os.mkdir(energy_dir)
FindTorsionEnergies(dft_dir,cwd)
print("Done gathering torsion energies!")
