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
	lines = []

        with open()
	cluster.runCommands(lines)

def print_usage(name):
        print "Usage : " + name + " <List of commands>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
