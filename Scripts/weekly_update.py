import urllib
import datetime
from rdkit import Chem
date = datetime.datetime.now().strftime("%m_%d_%Y")
name = date + "_release.txt"
smiles_file = date + "_release.smi"
urllib.urlretrieve("http://www.wwpdb.org/files/new_release_structure_nonpolymer.tsv", name)
with open(name, 'r') as f:
    lines = [line.split() for line in f]
lines = [line for line in lines if len(line) == 3]
for i,line in enumerate(lines):
    lines[i].append(Chem.inchi.MolFromInchi(line[2]))
#with open(smiles_file, 'w') as f:
for line in lines:
    print Chem.MolToSmiles(Chem.MolFromSmiles('CCOC'))
        #f.write(Chem.MolToSmiles(line[3]) + '\t' + line[1] + '\t' + line[0] + '\n')
