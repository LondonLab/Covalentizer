import sys,os,math
sys.path.append(os.environ["COVALENTIZER"])
from Code import *

def main(name, argv):
        if not len(argv) == 2:
                print_usage(name)
                return

        with open(argv[1], 'r') as f:
                lines = [line.split() for line in f]
        with open('filtered_res.txt', 'w') as f:
                for line in lines:
                        if PYMOLUtils.is_lig_single(argv[0], line[0], line[2]):
                                f.write('\t'.join(line) + '\n')

def print_usage(name):
        print "Usage : " + name + " <pdb> <res_file>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
