import glob
import subprocess
from subprocess import Popen, PIPE
import sys
import os
import shutil
sys.path.append(os.environ["COVALIB"])
from Code import *
from Code.Covalentizer import *

def main(name, argv):
        if len(argv) != 1:
                print_usage(name)
                return

        cluster = Cluster.Cluster()
        cluster.runBatchJobs(argv[0], "python /home/covalentizer/CovaLib/Scripts/Covalentizer/AutoProcess/process_pdb_multiple_lig.py rec.pdb ligands.txt")

def print_usage(name):
        print "Usage : " + name + " <dirlist>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
