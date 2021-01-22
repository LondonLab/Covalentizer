import random
import sys
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem import rdFMCS
import numpy as np

def main(name, argv):
    if len(argv) != 2:
        print_usage(name)
        return

    with open(argv[0], 'r') as f:
        mol1_f = [line.split() for line in f]
    end = mol1_f.index(['@<TRIPOS>BOND'])
    begin = mol1_f.index(['@<TRIPOS>ATOM'])
    mol1_f = mol1_f[begin + 1 : end]

    with open(argv[1], 'r') as f:
        mol2_f = [line.split() for line in f]
    end = mol2_f.index(['@<TRIPOS>BOND'])
    begin = mol2_f.index(['@<TRIPOS>ATOM'])
    mol2_f = mol2_f[begin + 1 : end]

    #with open(argv[2], 'r') as f:
    #    smiles = f.readline().split()[0]

    mol1 = Chem.MolFromMol2File(argv[0], sanitize=False)
    mol2 = Chem.MolFromMol2File(argv[1], sanitize=False)
    mols = [mol1, mol2]

    #print len(mol1.GetAtoms())
    #print len(mol2.GetAtoms())
    #mols = [Chem.RemoveHs(x) for x in mols]

    # Find the MCS:
    #mcs = rdFMCS.FindMCS(mols,threshold=0.8,completeRingsOnly=True,ringMatchesRingOnly=True)
    mcs = rdFMCS.FindMCS(mols, ringMatchesRingOnly=True, completeRingsOnly=True)
    patt = Chem.MolFromSmarts(mcs.smartsString)
    #patt = Chem.MolFromSmarts(smiles)
    refMol = mols[0]
    refMatch = refMol.GetSubstructMatch(patt)
    atoms_1 = [np.array([float(a) for a in mol1_f[i][2:5]]) for i in refMatch]
    #atomstag = [np.array([a for a in mol1_f[i]]) for i in refMatch]
    size = len(atoms_1)
    rmsVs = []
    for probeMol in mols[1:]:
        mv = probeMol.GetSubstructMatch(patt)
        atoms_2 = [np.array([float(a) for a in mol2_f[i][2:5]]) for i in mv]
        #atomstag2 = [np.array([a for a in mol2_f[i]]) for i in mv]
        rmsd = 0
        for i in range(size):
            rmsd += (np.linalg.norm(atoms_1[i] - atoms_2[i])**2)
        rmsd = np.sqrt(rmsd/size)
        print argv[1] + '\t' + str(rmsd)
        #rmsVs.append()
        #rms = AllChem.AlignMol(probeMol,refMol,atomMap=list(zip(mv,refMatch)))
        #rmsVs.append(rms)
    #print rmsVs


def print_usage(name):
    print "Usage : " + name + " <First Mol> <Second Mol>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
