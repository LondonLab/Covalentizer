import shutil
import os
import sys
import PYMOLUtils
import ChemUtils

def cysteine_folders(pdb_file):
    with open('res.txt', 'r') as f:
        lines = [line.split() for line in f]
    cysteines = set([(line[1], line[2]) for line in lines])
    for i,cys in enumerate(cysteines):
        os.mkdir("CYS" + str(i))
        os.chdir("CYS" + str(i))
        with open('res.txt', 'w') as f:
            for line in lines:
                if line[1] == cys[0] and line[2] == cys[1]:
                    f.write('\t'.join(line) + '\n')
        for line in lines:
            if line[1] == cys[0] and line[2] == cys[1]:
                os.mkdir(line[0])
                os.chdir(line[0])
                PYMOLUtils.seperate_rec_lig('../../' + pdb_file, line[0], line[2])
                cysteine_rotamers('rec.pdb', line)
                os.chdir('../')
        os.chdir('../')


def cysteine_rotamers(rec, line):
    PYMOLUtils.pymol_mutate(rec, line[2], line[1])
    with open('prepare.sh', 'w') as f:
        f.write('python ' + os.environ["SCRIPTS"] + '/Prepare.py rec.pdb xtal-lig.pdb CYS ' + line[1] + ' HG')
    os.system('chmod 777 prepare.sh')
    ChemUtils.convert('xtal-lig.pdb', 'xtal-lig.mol2')
    recs = ['rec.pdb'] + [rec for rec in os.listdir('./') if 'rec_' in rec]
    for i,rec in enumerate(recs):
        fol = 'ROT' + str(i) + '/'
        os.mkdir(fol)
        if i == 0:
            shutil.copy(rec, fol)
        else:
            os.rename(rec, fol + 'rec.pdb')
        shutil.copy('xtal-lig.pdb', fol)
        shutil.copy('prepare.sh', fol)
