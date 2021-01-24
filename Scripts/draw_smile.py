import sys,os,shutil
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
sys.path.append(os.environ["COVALENTIZER"])
from Code import *

smile_file = sys.argv[1]
name = smile_file.split('.')[0]
png = name + '.png'
mol2 = name + '.mol2'
scene = "scene_" + png
with open(sys.argv[1], 'r') as f:
    smile = f.readline().split()[0]
rmol = Chem.MolFromSmiles(smile)
Draw.MolToFile(rmol, png)
PYMOLUtils.scene_photo('rec.pdb', 'xtal-lig.pdb', mol2, scene)
