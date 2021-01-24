import subprocess
import sys
import os
sys.path.append(os.environ["COVALENTIZER"])
from Code import *

def main(name, argv):
        if len(argv) == 0:
                print_usage(name)
                return
	PyUtils.create_folder(argv[0][:-5])
	i = 0 #-1
	#if len(argv) == 1:
	#	i = 0
	inside = False
        for line in open(argv[0], 'r'):
		if not(inside) and line[0] == '#':
			i += 1
			inside = True
			if len(argv) == 1 or str(i) in argv[1:]:
				poses_f = open(argv[0][:-5] + '/' + str(i) + '.mol2', 'w')
		if not line[0] == '#':
			inside = False
		if len(argv) == 1 or str(i) in argv[1:]:
			poses_f.write(line)

def print_usage(name):
        print "Usage : " + name + " <poses file name> <index(s),(optional)>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
