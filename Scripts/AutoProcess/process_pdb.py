import glob
import subprocess
from subprocess import Popen, PIPE
import sys
import os
import shutil
sys.path.append(os.environ["COVALIB"])
from Code import *
from Code.Covalentizer import *

def main(name, argv):
        if len(argv) != 2:
                print_usage(name)
                return

        f = open('log.txt', 'w')
        f.write('Covalentizer has started.\n')

        curr_dir = os.getcwd()
        cluster = Cluster.Cluster()
        with open('lig.name', 'w') as lig_f:
                lig_f.write(argv[1] + '\n')
        f.write('Looking for cysteine\'s tiol within 6A of any of the ligand atoms.\n')
        res = PYMOLUtils.env_cysteine(argv[0], 'lig.name')
        with open('res.txt', 'w') as f:
                for r in res:
                        f.write("\t".join(r) + '\n')
        if len(res) == 0:
                f.write('Did not find cysteines which are close to the ligand.\n')
                f.close()
                return

        f.write('Create folders for the different cysteines and rotamers.\n')
        CysUtils.cysteine_folders(argv[0])
        dirlist = glob.glob("CYS*/*/ROT*/")
        with open('dirlist', 'w') as dirf:
                for line in dirlist:
                        dirf.write(line + '\n')
        f.write('Start preparing the structures for docking.\n')
        prepare_jobs = cluster.runJobs('dirlist', './prepare.sh')
        os.mkdir('Ligands/')
        ligands = {}
        for r in res:
                if r[0] in ligands:
                        continue
                ligands[r[0]] = curr_dir + '/Ligands/' + r[0] + '/gz_files/'
                os.mkdir('Ligands/' + r[0])
                PYMOLUtils.seperate_rec_lig(argv[0], r[0], r[2])
                ChemUtils.convert('xtal-lig.pdb', 'Ligands/' + r[0] + '/' + 'smile.smi')
                os.remove('xtal-lig.pdb')
        os.chdir('Ligands')
        for lig in ligands:
                os.chdir(lig)
                CovUtils.build_library('smile.smi', 'frags.smi', 'lib.smi')
                p = Popen([os.environ["DOCKBASE"] + '/ligand/generate/fixed/build_covalent_lib.csh', 'lib.smi', '10', 'LIB'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
                lig_jobs = [line for line in p.communicate()[0].split('\n') if 'pbs' in line]
                os.chdir('../')
        os.chdir('../')
        Cluster.Cluster.wait(prepare_jobs + lig_jobs)
        dock_jobs = []
        for d in dirlist:
                lig = d.split('/')[1]
                os.chdir(d)
                DOCK_Prepare.DOCK_Prepare.changeNumSave(10)
                p = Popen(['python', os.environ["SCRIPTS"] + '/DOCKovalentTask.py', 'CovLib', ligands[lig], 'True'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
                dock_jobs += [line for line in p.communicate()[0].split('\n') if 'pbs' in line]
                os.chdir(curr_dir)
        Cluster.Cluster.wait(dock_jobs)
        f.close()
        os.system(os.environ["SCRIPTS"] + '/Covalentizer/Final_Scripts/combine.sh')
        os.mkdir('Results')
        for d in dirlist:
                web_folder = d + '/CovLib/web_files/'
                new_folder = 'Results/' + d.replace('/', '_')[:-1] + '/' 
                if os.path.isdir(web_folder):
                        shutil.move(web_folder, new_folder)

def print_usage(name):
        print "Usage : " + name + " <pdb_file> <lig_name>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
