import sys,os,math
sys.path.append(os.environ["COVALENTIZER"])
from Code import *

def main(name, argv):
        if not len(argv) == 4:
                print_usage(name)
                return

        PYMOLUtils.scene_photo(argv[:3], argv[3])

def print_usage(name):
        print "Usage : " + name + " <rec> <lig> <covalentized> <picture_name>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
