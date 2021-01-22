import sys,os,shutil
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
sys.path.append(os.environ["COVALIB"])
from Code import *

with open('res.txt', 'r') as f:
    [mol, resi, chain] = f.readline().split()
with open('data.txt', 'r') as f:
    elecs = [line.split() for line in f]

text = "Ligand ID: " + mol + '\n'
text += "Cysteine Position: " + resi + '\n'
text += "Chain: " + chain + '\n'
smiles = []
i = 0
for e in elecs:
    if e[1] in smiles:
        continue
    else:
        i += 1
    f_name = e[1] + "_elec.smi"
    m_name = e[1] + '.mol2'
    shutil.copyfile(m_name, 'cand' + str(i) + '.mol2')
    with open(f_name, 'r') as f:
        smile = f.readline().split()[0]
    rmol = Chem.MolFromSmiles(smile)
    Draw.MolToFile(rmol, 'cand' + str(i) + '.png')
    PYMOLUtils.scene_photo('rec.pdb', 'xtal-lig.pdb', 'cand' + str(i) + '.mol2', 'scene' + str(i) + '.png')
    with open('cand' + str(i) + '.smi', 'w') as f:
        f.write(smile + '\n')
    text += "Electophilic Analog " + str(i) + ": " + smile + "\trmsd: " + str(round(float(e[2]),2)) + '\n'
with open('text.txt', 'w') as f:
    f.write(text)
