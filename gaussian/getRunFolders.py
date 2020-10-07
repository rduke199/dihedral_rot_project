import os

cwd = os.getcwd()
os.chdir('rotated_dihed')
home = os.getcwd()
mols = os.listdir(home)
for m in mols:
    if os.path.isdir(m):
        mpath = os.path.join(home,m)
        os.chdir(mpath)
        degs = os.listdir(mpath)
        for d in degs:
            if os.path.isdir(d):
                dpath = os.path.join(home,m,d)
                os.chdir(dpath)
                files = os.listdir(dpath)
                for f in files:
                    fpath = os.path.join(home,m,d,f)
                    unit_num = f[8]
                    if f.endswith(".log"):
                        break
                    if f.endswith(".gjf"):
                        txt_file = '{}/folders_to_run{}.txt'.format(cwd,unit_num)
                        with open(txt_file, "a+") as fn:
                            fn.write(dpath+"\n")
                os.chdir(mpath)
        os.chdir(home)
