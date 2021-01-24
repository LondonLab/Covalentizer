import sys,os,math
import shutil
sys.path.append(os.environ["COVALENTIZER"])
from Code import *

def main(name, argv):
        if not len(argv) == 2:
                print_usage(name)
                return

        with open(argv[1], 'r') as f:
                line = [line.split() for line in f][0]
        print (argv[0], line[2], line[1])
        PYMOLUtils.pymol_mutate(argv[0], line[2], line[1])
        with open('prepare.sh', 'w') as f:
                f.write('python /home/danielza/CovaLib/Scripts/Prepare.py rec.pdb xtal-lig.pdb CYS ' + line[1] + ' HG')
        os.system('chmod 777 prepare.sh')
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

def print_usage(name):
        print "Usage : " + name + " <rec_file> <res_file>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
