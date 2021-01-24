import subprocess
import sys
import os
sys.path.append(os.environ["COVALENTIZER"])
from Code import *

def main(name, argv):
        if len(argv) != 1:
                print_usage(name)
                return

	cluster = Cluster.Cluster()
        with open(argv[0], 'r') as f:
                lines = [line[:-1] for line in f]

	#cluster.runCommands(lines)
	cluster.runBatchCommands(lines, mem='4000mb')

def print_usage(name):
        print "Usage : " + name + " <command list>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
