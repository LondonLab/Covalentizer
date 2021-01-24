Covalentizer
This is the code for the Covalentizer protocol for creating covalent analogs
of reversible small molecules binders from the crystal structures.
The repository includes two folders, namely, Code and Scripts.

First, make an environment parameter for this folder called $COVALENTIZER

For using this repository you need to install the following programs and settings:

1. python 2.7.14.

2. Install RDKit, NumPy and pymol python packages. If you're using conda, use:
conda install -c rdkit rdkit
conda install -c anaconda numpy
conda install -c schrodinger pymol

3. Download and Install OpenBabel (http://openbabel.org/wiki/Main_Page)

4. For job scheduling our scripts rely on PBS (https://en.wikipedia.org/wiki/Portable_Batch_System)

5. Install DOCK and DOCKovalent, as well as all of their requirements (http://wiki.docking.org/index.php/Main_Page)

6. Set the following environment parameters: export OB=”/Path_to_OpenBabel/”, export COVALENTIZER=”/Path_to_the_git_folder/”, export PBS_HOME=”/Path_to_PBS_bin_folder/”, export DOCKBASE="/Path_to_DOCK_git/"

Usage:

The scripts to use Covalentizer are in Scripts/AutoProcess/.
There are two kinds of scripts there:
run*.py are scripts that send the job into the cluster directly.
process*.py are the jobs that the run*.py send. If you run them by themselves, they will run on your local node (not recommended).
The run*.py scripts have their appropriate process*.py scripts, usually with the same suffix. Here I will only elaborate on the run*.py scripts. The process*.py scripts take the same input.
lig_file has \n as delimiter.

python run_covalentizer_multiple_lig.py <pdb_file> <lig_file>

Explanation: This is the version run by the webserver. Pdb_file is the X-ray structure which includes the reversible ligand. Lig_file is a file with the names of the ligand residues (can be more than one) to consider, one in each line. You don’t need to specify the cysteine. It will detect it automatically.

python run_covalentizer.py <pdb_file> <lig_name>

Explanation: Obsolete version, which gets a single lig_name instead of a file.

python run_covalentizer_linker.py

Explanation: Not used by the webserver, but it’s useful nonetheless. This is for running an extended run with the di-amine linker library.

python run_covalentizer_small.py

Explanation: Also useful. The default of Covalentizer is to ignore ligand which are smaller than 300 D (to avoid crystallographic agents, etc.). When using this variant, it wouldn’t ignore them.

python run_covalentizer_small_linker.py

Explanation: A combination of the “linker” and “small” variants.

Inside the Results zip folder, there will be separate folder for each triad of cystein, ligand and rotamer.
For example: CYS432_E80_ROT3/
Inside each such folder, there will be the following files:
rec.pdb - the protein structure.
xtal-lig.pdb - the original reversible bound molecule.
chain.txt - the protein chain.
cys_position.txt - the position of the cysteine residue within the protein sequence.
Then, for each candidate, the following files (X is the candidate number ID):
X.mol2 - structure of the docked candidate.
X.rmsd - the RMSD of the maximum common substructure between the reversible molecule and the covalent candidate.
X.png - picture of the 2D representation of the covalent candidate.
scene_X.png - picture of the candidate bound to the cysteine.

Thanks for using Covalentizer.
