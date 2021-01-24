import glob
import subprocess
from subprocess import Popen, PIPE
import sys
import os
import shutil
sys.path.append(os.environ["COVALENTIZER"])
from Code import *

def main(name, argv):
        if len(argv) != 2:
                print_usage(name)
                return

        cluster = Cluster.Cluster()
        cluster.runSingle("python " + os.environ["COVALENTIZER"] + "/Scripts/AutoProcess/process_pdb_small_linker.py " + argv[0] + " " + argv[1])

def print_usage(name):
        print "Usage : " + name + " <pdb_file> <lig_file>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
