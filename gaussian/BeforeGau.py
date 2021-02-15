from ocelot.task.confgen import ConfGen
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem
from itertools import combinations_with_replacement, product
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput

from ocelot.task.wtuning import WtuningJob
import subprocess
import re
import os
import shutil

def FindCntDihed(mol):
    """
    function for recording central dihedral atom numbers
    """
    #find all potential center bonds (single bonds between carbons)
    potential_cntbond = []
    for bond in mol.GetBonds():
        if str(bond.GetBondType()) == 'SINGLE':
            atom_a = bond.GetBeginAtom()
            atom_b = bond.GetEndAtom()
            bond_atoms = [atom_a, atom_b]
            if atom_a.GetAtomicNum() == 6 and atom_b.GetAtomicNum() == 6:
                potential_cntbond.append(bond_atoms)

    #find central bond
    num_cnt_bond = int((len(potential_cntbond)-1)/2)
    cnt_bond = potential_cntbond[num_cnt_bond]

    cntatom_a = cnt_bond[0]
    cntatom_b = cnt_bond[1]
    dihed1 = []
    dihed2 = []

    #assemble list of atoms in first dihedral angle
    for n in cntatom_a.GetNeighbors():
        if n.GetIdx() != cntatom_b.GetIdx():
            dihed1.append(n.GetIdx())
            dihed1.extend((cntatom_a.GetIdx(),cntatom_b.GetIdx()))
            break
    for n in cntatom_b.GetNeighbors():
        if n.GetIdx() != cntatom_a.GetIdx():
            dihed1.append(n.GetIdx())
            break
    #assemble list of atoms in second dihedral angle
    for n in cntatom_a.GetNeighbors():
        dihed2.append(n.GetIdx()) if n.GetIdx() not in dihed1 else dihed1
    dihed2.extend((cntatom_a.GetIdx(),cntatom_b.GetIdx()))
    for n in cntatom_b.GetNeighbors():
        dihed2.append(n.GetIdx()) if n.GetIdx() not in dihed1 else dihed1

    dihedral1 = str(dihed1[0])+" "+str(dihed1[1])+" "+str(dihed1[2])+" "+str(dihed1[3])+" "
    dihedral2 = str(dihed2[0])+" "+str(dihed2[1])+" "+str(dihed2[2])+" "+str(dihed2[3])+" "

    return dihedral1+','+dihedral2

def GenPolymers():
    """
    Generate polymers with substituent combinations and all their confromers
    """
    #These are the base SMILES for the different systems, 66 is biphenyl, 65 is phenyl thiophene, and 55 is bithiophene
    #I put in place holder characters so that they could be replaced later on, so W is NOT tungsten
    r_55 = ['C1=CX=C(S1)C2=YC=C(S2)','C1=CX=C(S1)C2=YC=C(S2)C1=CX=C(S1)C2=YC=C(S2)C1=CX=C(S1)C2=YC=C(S2)',
            'C1=CX=C(S1)C2=YC=C(S2)C1=CX=C(S1)C2=YC=C(S2)C1=CX=C(S1)C2=YC=C(S2)C1=CX=C(S1)C2=YC=C(S2)C1=CX=C(S1)C2=YC=C(S2)',
            'C1=CX=C(S1)C2=YC=C(S2)C1=CX=C(S1)C2=YC=C(S2)C1=CX=C(S1)C2=YC=C(S2)C1=CX=C(S1)C2=YC=C(S2)C1=CX=C(S1)C2=YC=C(S2)C1=CX=C(S1)C2=YC=C(S2)C1=CX=C(S1)C2=YC=C(S2)']

    r_65 = ['c1cXc(Yc1)C2=ZC=C(S2)','c1cXc(Yc1)C2=ZC=C(S2)c1cXc(Yc1)C2=ZC=C(S2)c1cXc(Yc1)C2=ZC=C(S2)',
            'c1cXc(Yc1)C2=ZC=C(S2)c1cXc(Yc1)C2=ZC=C(S2)c1cXc(Yc1)C2=ZC=C(S2)c1cXc(Yc1)C2=ZC=C(S2)c1cXc(Yc1)C2=ZC=C(S2)',
            'c1cXc(Yc1)C2=ZC=C(S2)c1cXc(Yc1)C2=ZC=C(S2)c1cXc(Yc1)C2=ZC=C(S2)c1cXc(Yc1)C2=ZC=C(S2)c1cXc(Yc1)C2=ZC=C(S2)c1cXc(Yc1)C2=ZC=C(S2)c1cXc(Yc1)C2=ZC=C(S2)']

    r_66 = ['c1cYc(Xc1)c2Wcc(cZ2)','c1cYc(Xc1)c2Wcc(cZ2)c1cYc(Xc1)c2Wcc(cZ2)c1cYc(Xc1)c2Wcc(cZ2)',
            'c1cYc(Xc1)c2Wcc(cZ2)c1cYc(Xc1)c2Wcc(cZ2)c1cYc(Xc1)c2Wcc(cZ2)c1cYc(Xc1)c2Wcc(cZ2)c1cYc(Xc1)c2Wcc(cZ2)',
            'c1cYc(Xc1)c2Wcc(cZ2)c1cYc(Xc1)c2Wcc(cZ2)c1cYc(Xc1)c2Wcc(cZ2)c1cYc(Xc1)c2Wcc(cZ2)c1cYc(Xc1)c2Wcc(cZ2)c1cYc(Xc1)c2Wcc(cZ2)c1cYc(Xc1)c2Wcc(cZ2)']

    #This part makes lists with the different combinations of subs, these can be plain lists as in 55, or lists of lists
    X_55 = list(combinations_with_replacement(['C','C(F)','N','C(OC)'],2))
    X_65 = list(product(list(combinations_with_replacement(['c','n','c(F)','c(OC)'],2)),['C','N','C(F)','C(OC)']))
    X_66 = list(product(list(combinations_with_replacement(['c','c(F)','c(OC)'],2)),list(combinations_with_replacement(['c','n'],2))))

    base_ring = [55,65,66]

    for r in base_ring:
        #Sup is how many repeat units there are (1,3,5,7)
        sup = 1
        for h in eval('r_'+str(r)):
            #Fn is the iteration of a polymer with a specific combination of subs
            fn=0
            for i in eval('X_'+str(r)):
                #The following 'smiles' variable is equal to the substituted SMILES string
                if r == 55:
                    smiles = h.replace('X',i[0]).replace('Y',i[1])
                    XYZ = '{}-{}'.format(i[0],i[1]).replace('(','').replace(')','')
                if r == 65:
                    smiles = h.replace('X',i[0][0]).replace('Y',i[0][1]).replace('Z',i[1])
                    XYZ = '{}-{}-{}'.format(i[0][0],i[0][1],i[1]).replace('(','').replace(')','')
                if r == 66:
                    smiles = h.replace('W',i[1][0]).replace('X',i[1][1]).replace('Y',i[0][0]).replace('Z',i[0][1])
                    XYZ = '{}-{}-{}-{}'.format(i[1][0],i[1][1],i[0][0],i[0][1]).replace('(','').replace(')','')
                with open(name+'.smi','a') as fout:
                    fout.write(str(smiles))
                mol = Chem.MolFromSmiles(smiles)
                molH = AllChem.AddHs(mol)
                #Name: 'mols_<ring type>_<#repeat units>_<polymer#>_<subs>'
                name = 'mols_{}_{}_{:02d}_{}'.format(r,sup,fn,XYZ)
                #recording the dihedral atom numbers
                FindCntDihed(molH)
                with open(name+'_dihed.txt','a') as fout:
                        fout.write(str(FindCntDihed(molH)) + '\n')
                #You might need to adjust this next line, I kind of just picked a 'prune_in_gen' value. This is what makes the confs
                c = ConfGen(m=molH,mol_name=name,prune_in_gen=0.5)
                c.genconfs(write_confs=True)
                fn+= 1
            sup+=2

def XYZtoGJF(fn):
    """
    Convert an individually imputted xyz file to a gjf Gaussian input file
    """
    mol = Molecule.from_file(fn)
    gau = GaussianInput(mol=mol,charge=0,spin_multiplicity=1,functional='PM6',basis_set='',route_parameters={'':''},link0_parameters={'%mem':'5GB','%nproc':'1','%chk':'{}.chk'.format(fn.split('.')[0])})
    gjf_file = gau.write_file('{}.gjf'.format(fn.split('.')[0]))
    return gjf_file

def RunJob(job_name,job):
    """
    Submit any SLURM job. job_name and job should each be a string.
    """
    runfile = os.path.join(os.getcwd(),"%s.job" %job_name)
    with open(runfile,"w+") as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH -t 14-00:00:00\n")
        fh.writelines("#SBATCH --job-name=%s\n" % job_name)
        fh.writelines("#SBATCH -N 1\n")
        fh.writelines("#SBATCH -n 8\n")
        fh.writelines("#SBATCH -p SKY32M192_L\n")
        fh.writelines("#SBATCH --account=col_cmri235_uksr\n")
        fh.writelines("#SBATCH --mail-type=ALL\n")
        fh.writelines("#SBATCH --mail-user=rdu230@uky.edu\n")
        fh.writelines("#SBATCH --error=SLURM_JOB_%j.err\n")
        fh.writelines("#SBATCH --output=SLURM_JOB_%j.out\n" )

        fh.writelines("ulimit -u 32768\n")
        fh.writelines("export OPENBLAS_NUM_THREADS=1\n")

        fh.writelines("module load ccs/gaussian/g16-A.03/g16-sandybridge\n")
        fh.writelines("source $(conda info --base)/etc/profile.d/conda.sh\n")
        fh.writelines("conda activate myenv\n")
        fh.writelines("echo 'Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST '\n")
        fh.writelines(". ~/.bashrc\n")

        fh.writelines(job)

    os.system("sbatch %s" %runfile)

def MoveFiles(destination,starts_with="",ends_with=""):
    """
    Moves all files that start/end with a certain string to a new destination
    """
    for i in os.listdir(os.getcwd()):
        if i.startswith(starts_with) and i.endswith(ends_with):
            shutil.move(i,destination+i)

def FuncAllFiles(func,starts_with="",ends_with=""):
    """
    Runs all files that start/end with a certain string through a function
    """
    for i in os.listdir(os.getcwd()):
        if i.startswith(starts_with) and i.endswith(ends_with):
            func(i)


home = os.getcwd()
os.mkdir('xyz')
os.mkdir('gjf')
os.mkdir('gau_output')
os.mkdir('wtuning')
os.mkdir('diheds')
os.mkdir('smiles')
os.mkdir('energies')


os.chdir(home+'/xyz')
GenPolymers()
MoveFiles(home+'/diheds/',ends_with='dihed.txt')
MoveFiles(home+'/smiles/',ends_with='.smi')
print("Done generating xyz")

FuncAllFiles(XYZtoGJF,ends_with='.xyz')
MoveFiles(home+"/gjf/",ends_with='.gjf')
print("Done converting xyz to gjf")

os.chdir(home+'/gjf')
for i in os.listdir(os.getcwd()):
    os.mkdir(home+'/gau_output/'+i[:-4])
    subprocess.call('cp {} {}/gau_output/{}/'.format(i,home,i[:-4]),shell=True)
    os.chdir(home+'/gau_output/'+i[:-4])
    # RunJob('gau{}'.format(i[5:12]),"g16 {}".format(i))
    os.chdir(home+'/gjf')
    print("Submitted {} gaussian job.".format(i))
print("Done submitting Gaussian job")
