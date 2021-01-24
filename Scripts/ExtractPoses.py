import subprocess
import sys
import os
sys.path.append(os.environ["COVALENTIZER"])
from Code import *

def main(name, argv):
        if len(argv) != 2:
                print_usage(name)
                return
	PyUtils.create_folder('poses')
	i = 0
	with open(argv[1], 'r') as f:
		indices = [line.split()[0] for line in f.readlines()]
	inside = False
        for line in open(argv[0], 'r'):
		if not(inside) and line[0] == '#':
			i += 1
			inside = True
			if str(i) in indices:
				poses_f = open('poses/' + str(i) + '.mol2', 'w')
		if not line[0] == '#':
			inside = False
		if str(i) in indices:
			poses_f.write(line)

def print_usage(name):
        print "Usage : " + name + " <poses file name> <indices file>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
