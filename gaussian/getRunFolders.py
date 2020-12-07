import os
import re

cwd = os.getcwd()
os.chdir('rotated_dihed')
home = os.getcwd()
mols = [m for m in os.listdir(home) if os.path.isdir(m)]
for m in mols:
    mpath = os.path.join(home,m)
    os.chdir(mpath)
    degs = [d for d in os.listdir(mpath) if os.path.isdir(d)]
    for d in degs:
        dpath = os.path.join(home,m,d)
        os.chdir(dpath)
        log_files = [x for x in os.listdir(dpath) if x.endswith('.log')]
        gjf_files = [x for x in os.listdir(dpath) if x.endswith('.gjf')]
        if len(log_files) == 1:
            normal = 0
            with open(log_files[0]) as fn:
                last_line = fn.readlines()[-1]
                if re.search('Normal termination', last_line):
                    normal = 1
            if normal == 0:
                print("Error. {} dihedral rotation did not finish.".format(d))
                pass
            else:
                continue
        if len(gjf_files) == 0:
            print("Error. {} does not contain a gjf file.".format(d))
            continue
        else:
            unit_num = gjf_files[0][8]
            txt_file = '{}/test/folders_to_run{}.txt'.format(cwd,unit_num)
            with open(txt_file, "a+") as fn:
                fn.write(dpath+"\n")
            continue
        os.chdir(mpath)
        break
    os.chdir(home)
