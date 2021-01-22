import random
import sys
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem

def main(name, argv):
    if len(argv) != 3:
        print_usage(name)
        return

    main_function(argv)

def main_function(argv):
    with open(argv[0], 'r') as f:
        reactions = [line.split() for line in f.readlines()]
    rxns = [[rdChemReactions.ReactionFromSmarts(r[0]), r[1]] for r in reactions]

    #read molport building blocks
    with open(argv[1], 'r') as f:
        lines_1 = [line.split()[0] for line in f.readlines()]

    #creating RDKit Mol objects
    molecules_1 = [Chem.MolFromSmiles(line) for line in lines_1]
    molecules_1 = [m for m in molecules_1 if m is not None]

    f = open(argv[2], 'w')
    write_mols = []
    # calculate reactions
    for i, rxn in enumerate(rxns):
        #find all valid molecules for a specific reaction
        reactents_smarts = rxn[0].GetReactants()
        has_patt_1 = [m.HasSubstructMatch(reactents_smarts[0]) for m in molecules_1]
        reactent_option_1 = [m for i,m in enumerate(molecules_1) if has_patt_1[i]]

        for j,mol1 in enumerate(reactent_option_1):
            #calculate the product for the reaction for these two compounds
            products = rxn[0].RunReactants((mol1,))
            # write to file
            for i, pro in enumerate(products):
                write_mols.append((Chem.MolToSmiles(pro[0], isomericSmiles=True), str(j+1) + '_' + str(i+1) + '_' + rxn[1]))
    uniq_smiles = []
    for m in write_mols:
        if not m[0] in uniq_smiles:
            f.write(m[0] + '\t' + m[1] + '\n')
            uniq_smiles.append(m[0])
    f.close()

def print_usage(name):
    print "Usage : " + name + " <reactions> <library> <output>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
