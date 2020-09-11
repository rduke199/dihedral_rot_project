import os
import re
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

def FindTorsionEnergies():
	"""
	Iterates through dihedral directory and finds energies for each torsion angle
	for each molecule. It creates a energies directory with a energy cvs file
	for each molecule.
	"""
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
					try:
						gen = [x for x in os.listdir(os.getcwd()) if x.endswith('.log')]
						for log_file in gen:
							with open(log_file) as fn:
								for line in fn:
									if re.search(r'SCF Done',line):
										a = line.split()[4]
						with open('{}/energies/energies_{}.cvs'.format(cwd,m),'a') as fout:
							fout.write(d.split('_')[-1][:-3] + ',' + a + '\n')
						print("Dihedral rotations energies collected for {}!".format(d))
					except:
						print(".log file not found for {}. Energy calculations for dihedral rotations did not finish!".format(d))
						continue
				os.chdir(mpath)
		os.chdir(home)

def SortEnergies():
    """
    In directory with all energie .cvs files, this copies all energies files into
    a new "sorted_energies" directory. This contains directories for each unique
    monomer base which each contains energies files for all polymer lengths of
    this monomer base.
    """
    home = os.getcwd()
    try:
        os.mkdir(home+"/sorted_energies/")
    except:
        pass
    list_attrib = []
    #Recall naming  :  mols_<ring type>_<#repeat units>_<polymer#>_<subs>_<confromer#>
    for file in os.listdir(os.getcwd()):
        mol = str(file)
        try:
            ring_num = int(mol.split('_')[2])
            unit_num = int(mol.split('_')[3])
            polymer_num = int(mol.split('_')[4])
            substituents = mol.split('_')[5]
            list_attrib.append([ring_num,unit_num,polymer_num,substituents])
        except:
            pass
    ring_type = [55,65,66]
    unit_nums = [1,3,5,7]
    for r in ring_type:
        i=0
        while i <= 100:
            for j in list_attrib:
                if j[0] == r and j[2] == i:
                    mol_name = '{}/energies_mols_{}_{}_{:02d}_{}'.format(home,j[0],j[1],j[2],j[3])
                    dir_name = '{}/sorted_energies/allenergies_mols_{}_{:02d}_{}/'.format(home,j[0],j[2],j[3][:-4])
                    try:
                        os.mkdir(dir_name)
                    except:
                        pass
                    subprocess.call('cp {} {}'.format(mol_name,dir_name),shell=True)
            i+=1

def MakeTorsionPlot(energy_file):
	"""
	Plots the energy vs the torsion angle for a molecule given a torsion energy
	file. Saves the plot to a png file.
	"""
	mol_name = energy_file.strip('energies_ _0.cvs')
	edata = np.genfromtxt(fname=energy_file,delimiter=',', dtype='unicode')
	edata=edata.astype(np.float)

	energy_values = edata[:,1]
	phi = edata[:,0]

	energy_range = max(energy_values) - min(energy_values)

	plt.scatter(phi, energy_values, color='MediumVioletRed')
	plt.xlim(min(phi)-3, max(phi)+3)
	plt.xticks(np.linspace(start=0, stop=180, num=7))
	plt.ylim(top = max(energy_values) + energy_range*0.15,bottom = min(energy_values) - energy_range*0.15)
	plt.yticks(np.linspace(start=min(energy_values), stop=max(energy_values), num=5))
	plt.xlabel("dihedral angle in degrees")
	plt.ylabel("energy (Hartrees)")
	plt.title("{} torsion energy".format(mol_name))
	plt.savefig('torsionE_plt_{}.png'.format(mol_name),dpi=300)
	plt.close()

def OverlayPlts():
	"""
	Makes overlay plot of allenergies files in cwd.
	"""
	home = os.getcwd()
	mol_name = home.split('/')[-1].strip('allenergies_')
	plt.figure()
	gen = [x for x in os.listdir(os.getcwd()) if x.endswith('.cvs')]
	for energy_file in gen:
	    poly_len = str(energy_file).split('_')[3]
	    edata = np.genfromtxt(fname=energy_file,delimiter=',', dtype='unicode')
	    edata=edata.astype(np.float)

	    i=0
	    while i < edata.shape[0]:
	        if edata[i,0] == 0:
	            zero_energy = edata[i,1]
	        i+=1

	    energy_values = edata[:,1]
	    phi = edata[:,0]

	    plt.scatter(phi, energy_values-zero_energy, label=poly_len)

	plt.xlim(-3, 183)
	plt.xticks(np.linspace(start=0, stop=180, num=7))
	plt.xlabel("dihedral angle in degrees")
	plt.ylabel("energy (Hartrees)")
	plt.title("{} torsion energy".format(mol_name))
	plt.legend()
	plt.savefig('overlay_plt_{}.png'.format(mol_name),dpi=300)
	plt.close()

def MoveFiles(destination,starts_with="",ends_with=""):
	"""
	Moves all files that start/end with a certain string to a new destination
	"""
	for i in os.listdir(os.getcwd()):
		if i.startswith(starts_with) and i.endswith(ends_with):
			shutil.move(i,destination+i)

cwd = os.getcwd()
FindTorsionEnergies()
print("Done gathering torsion energies!")

os.chdir(cwd+'/energies')
for f in os.listdir(os.getcwd()):
	if f.startswith("energies"):
		try:
			MakeTorsionPlot(f)
		except:
			pass
os.mkdir('plots')
MoveFiles(cwd+'/energies/plots/',starts_with="torsionE")
print("Done making potential energy plots!")

os.chdir(cwd+"/energies")
SortEnergies()
os.mkdir(cwd+"/energies/overlay_plts")
os.chdir(cwd+"/energies/sorted_energies/")
for mol in os.listdir(os.getcwd()):
    if mol.startswith("allenergies_"):
        print(mol)
        os.chdir("{}/energies/sorted_energies/{}".format(cwd,mol))
        OverlayPlts()
        MoveFiles(cwd+"/energies/overlay_plts/",ends_with='.png')
    else:
        pass
print("Done makeing overlay plots!")
