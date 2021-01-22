import subprocess
import sys
import os
sys.path.append(os.environ["COVALIB"])
from Code import *

def main(name, argv):
        if len(argv) != 2:
                print_usage(name)
                return

	cluster = Cluster.Cluster()
	#cluster.runJobs(argv[0], argv[1])
	cluster.runBatchJobs(argv[0], argv[1], mem='10000mb')

def print_usage(name):
        print "Usage : " + name + " <dirlist> <command>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
