import random
import sys
import os
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import Recap
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def multi_react(argv):
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

def multi_linkers(argv, linker_smiles):
    with open(argv[0], 'r') as f:
        reactions = [line.split() for line in f.readlines()]
    rxns = [[rdChemReactions.ReactionFromSmarts(r[0]), r[1]] for r in reactions]

    #read molport building blocks                                                                                                                
    with open(argv[1], 'r') as f:
        lines_1 = [line.split()[0] for line in f.readlines()]
    with open(linker_smiles, 'r') as f:
        lines_2 = [line.split() for line in f.readlines()]

    #creating RDKit Mol objects                                                                                                                  
    molecules_1 = [Chem.MolFromSmiles(line) for line in lines_1]
    molecules_1 = [m for m in molecules_1 if m is not None]
    molecules_2 = [(Chem.MolFromSmiles(line[0]), line[1]) for line in lines_2]
    molecules_2 = [m for m in molecules_2 if m[0] is not None]

    f = open(argv[2], 'w')
    write_mols = []
    # calculate reactions
                                                                                                        
    for i, rxn in enumerate(rxns):
        #find all valid molecules for a specific reaction                                                                                        
        reactents_smarts = rxn[0].GetReactants()
        has_patt_1 = [m.HasSubstructMatch(reactents_smarts[0]) for m in molecules_1]
        reactent_option_1 = [m for i,m in enumerate(molecules_1) if has_patt_1[i]]
        has_patt_2 = [m[0].HasSubstructMatch(reactents_smarts[1]) for m in molecules_2]
        reactent_option_2 = [m for i,m in enumerate(molecules_2) if has_patt_2[i]]

        for j,mol1 in enumerate(reactent_option_1):
        #calculate the product for the reaction for these two compounds                                                                      
            for k,mol2 in enumerate(reactent_option_2):
            #calculate the product for the reaction for these two compounds                                                                     
                products = rxn[0].RunReactants((mol1,mol2[0]))
                # write to file                                                                                                                 
                for i, pro in enumerate(products):
                    write_mols.append((Chem.MolToSmiles(pro[0], isomericSmiles=True), str(j+1) + '_' + str(i+1) + '_' + mol2[1] + '_' + rxn[1]))

    uniq_smiles = []
    for m in write_mols:
        if not m[0] in uniq_smiles:
            f.write(m[0] + '\t' + m[1] + '\n')
            uniq_smiles.append(m[0])
    f.close()

def build_library(in_smile, frags, lib, rules = os.environ["COVALIB"] + "/Code/Covalentizer/numbered_reaction.re", linker_lib = False, linker_smiles = ''):
    argv = [in_smile, frags, lib]
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
            if new_size > 7:
                if new_size <= 25 and rdMolDescriptors.CalcNumRotatableBonds(new) <= 5:
                    children.append(Chem.MolToSmiles(new, isomericSmiles=True))
                core_smile = MurckoScaffold.MurckoScaffoldSmilesFromSmiles(new_smiles, includeChirality = True)
                core = Chem.MolFromSmiles(core_smile)
                if new_size <= 25 and rdMolDescriptors.CalcNumRotatableBonds(core) <= 5 and core.GetNumHeavyAtoms() > 7:
                    children.append(core_smile)
    with open(argv[1], 'w') as f:
        i = 1
        for m in set(children):
            if len(m) > 0:
                f.write(m + '\t' + str(i) + '\n')
                i += 1
    
    if not linker_lib:
        multi_react([rules, argv[1], argv[2]])
    else:
        multi_linkers([rules, argv[1], argv[2]], linker_smiles)

def get_MW(smile_f):
    with open(smile_f, 'r') as f:
        smile = f.readline().split()[0]
    mol = Chem.MolFromSmiles(smile)
    if mol == None:
        return None
    return Descriptors.MolWt(mol)

