import sys,os,math
sys.path.append(os.environ["COVALIB"])
from Code import *

def main(name, argv):
        if not len(argv) == 2:
                print_usage(name)
                return

        with open(argv[1], 'r') as f:
                lines = [line.split() for line in f]
        for line in lines:
                os.mkdir(line[0])
                os.chdir(line[0])
                PYMOLUtils.seperate_rec_lig(argv[0], line[0], line[2])
                os.chdir('../')

def print_usage(name):
        print "Usage : " + name + " <pdb> <res_file>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
