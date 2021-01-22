import random
import sys
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import Recap
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdMolDescriptors
sys.path.append("/home/danielza/Work/Current/Covalentizer/scripts/Final_Scripts/")
from multi_elec import main_function

def main(name, argv):
    if len(argv) != 3:
        print_usage(name)
        return

    with open(argv[0], 'r') as f:
        smile = f.readline().split()[0]
    mol = Chem.MolFromSmiles(smile)
    if mol == None:
        return
    size = mol.GetNumHeavyAtoms()
    hierarch = Recap.RecapDecompose(mol)
    children = []
    for child in hierarch.GetAllChildren().keys() + [smile]:
        new_smiles = child.replace('[*]', '[H]')
        new = Chem.MolFromSmiles(new_smiles)
        if not new == None:
            new_size = new.GetNumHeavyAtoms()
            if new_size > 7 and new_size <= 25:
                if rdMolDescriptors.CalcNumRotatableBonds(new) <= 5:
                    children.append(Chem.MolToSmiles(new, isomericSmiles=True))
                    #children.append(new_smiles)
                core_smile = MurckoScaffold.MurckoScaffoldSmilesFromSmiles(new_smiles, includeChirality = True)
                core = Chem.MolFromSmiles(core_smile)
                if rdMolDescriptors.CalcNumRotatableBonds(core) <= 5 and core.GetNumHeavyAtoms() > 7:
                    children.append(core_smile)
    with open(argv[1], 'w') as f:
        i = 1
        for m in set(children):
            if len(m) > 0:
                f.write(m + '\t' + str(i) + '\n')
                i += 1
    
    main_function(["/home/danielza/Work/Current/Covalentizer/scripts/Final_Scripts/numbered_reaction.re", argv[1], argv[2]])
    
def print_usage(name):
    print "Usage : " + name + " <mol> <frags_file> <lib_file>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
