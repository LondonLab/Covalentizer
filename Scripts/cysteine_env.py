import sys,os,math
sys.path.append(os.environ["COVALENTIZER"])
from Code import *

def main(name, argv):
        if not len(argv) == 2:
                print_usage(name)
                return

        PYMOLUtils.env_cysteine(*argv)

def print_usage(name):
        print "Usage : " + name + " <pdb> <allowed_hets>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
