import subprocess
import sys
import os
sys.path.append(os.environ["COVALIB"])
from Code import *

def main(name, argv):
        if len(argv) != 1:
                print_usage(name)
                return

	cluster = Cluster.Cluster()
        line = 'python /home/danielza/CovaLib/Scripts/Combine.py cov_lib7/'

	cluster.runJobs(argv[0], line)

def print_usage(name):
        print "Usage : " + name + " <dirlist>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
