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
		mpath = os.path.join(rotdihed_dir,m)
		degs = [d for d in os.listdir(mpath) if os.path.isdir(os.path.join(mpath,d))]
		for d in degs:
			try:
				dpath = os.path.join(mpath,d)
				log_file = [x for x in os.listdir(dpath) if x.endswith('.log')][0]
				log_path = os.path.join(dpath,log_file)
				mol = GaussianOutput(log_path)
				if mol.properly_terminated:
					with open(log_path) as fn:
						log = fn.readlines()
						normal = 0
						for line in log:
						    if re.search(' Optimized Parameters',line):
						        normal = 1
						        break
						    else:
						        continue
						if normal == 0:
						    print("Error. {} job had a normal termination but did not have optimized parameters.".format(d))
				else:
				    normal = 0
				    print("Error. {} job did not have a normal termination.".format(d))
				if normal == 1:
					energy = mol.final_energy
					with open('{}/energies_{}.cvs'.format(out_dir,m),'a') as fout:
						fout.write(d.split('_')[-1][:-3] + ',' + str(energy) + '\n')
						print("Dihedral rotations energies collected for {}!".format(d))
			except:
				print("Error finding energy for {}.".format(d))

cwd = os.getcwd()
dft_dir = os.path.join(cwd,'dft/')
energy_dir = os.path.join(cwd,'energies/')
if not os.path.isdir(energy_dir): os.mkdir(energy_dir)
FindTorsionEnergies(dft_dir,energy_dir)
print("Done gathering torsion energies!")
