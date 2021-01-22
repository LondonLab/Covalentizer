import glob
import subprocess
from subprocess import Popen, PIPE
import sys
import os
import shutil
import glob
sys.path.append(os.environ["COVALIB"])
from Code import *
from Code.Covalentizer import *

def main(name, argv):
        if len(argv) != 2:
                print_usage(name)
                return

        curr_dir = os.getcwd()
        cluster = Cluster.Cluster()
        f = open('log.txt', 'w')
        lines = []
        lig_dict = {}
        with open(argv[1], 'r') as lig_f:
                for lig in lig_f:
                        split_line = lig.split()
                        lig_dict[split_line[1]] = split_line[0]
        lig_names = 'ligands.names'
        with open(lig_names, 'w') as lig_f:
                for lig in lig_dict:
                        lig_f.write(lig + '\n')
        for pdb in open(argv[0], 'r'):
                lines.append("http://www.rcsb.org/pdb/files/" + pdb[:-1] + ".pdb")
        download_pdb = cluster.runCommandsArgs("wget", lines)
        Cluster.Cluster.wait(download_pdb)
        for i in glob.glob('un_sub_l.*'):
                os.remove(i)
        os.mkdir('PDB')
        for i in glob.glob('*.pdb'):
                new_fol = 'PDB/' + i.split('.')[0]
                os.mkdir(new_fol)
                os.rename(i, new_fol + '/' + i)
        os.chdir('PDB')
        ligands = {}
        for i in glob.glob('*/'):
                print i
                pdb_file = i.split('/')[0] + '.pdb'
                os.chdir(i)
                res = PYMOLUtils.env_cysteine(pdb_file, '../../' + lig_names)
                with open('res.txt', 'w') as f:
                        for r in res:
                                if not r[0] in ligands:
                                        ligands[r[0]] = 'Ligands/' + r[0] + '/gz_files/'
                                f.write("\t".join(r) + '\n')
                CysUtils.cysteine_folders(pdb_file)
                os.chdir('../')
        os.chdir('../')
        dirlist = glob.glob("PDB/*/CYS*/*/ROT*/")
        with open('dirlist', 'w') as dirf:
                for line in dirlist:
                        dirf.write(line + '\n')
        prepare_jobs = cluster.runJobs('dirlist', './prepare.sh')
        os.mkdir('Ligands/')
        lig_jobs = []
        os.chdir('Ligands')
        for lig in ligands:
                print lig
                os.mkdir(lig)
                os.chdir(lig)
                with open('smile.smi', 'w') as smile_f:
                        smile_f.write(lig_dict[lig] + '\n')
                CovUtils.build_library('smile.smi', 'frags.smi', 'lib.smi')
                if not os.path.exists('lib.smi'):
                        os.chdir('../')
                        continue
                p = Popen([os.environ["DOCKBASE"] + '/ligand/generate/fixed/build_covalent_lib.csh', 'lib.smi', '10', 'LIB'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
                lig_jobs += [line for line in p.communicate()[0].split('\n') if 'pbs' in line]
                os.chdir('../')
        os.chdir('../')
        Cluster.Cluster.wait(prepare_jobs + lig_jobs, 1000)
        dock_jobs = []
        for d in dirlist:
                lig = d.split('/')[-3]
                if len(os.listdir(ligands[lig])) == 0 or len(glob.glob(os.path.dirname(ligands[lig][:-1]) + '/LIB*')) > 0:
                        continue
                os.chdir(d)
                DOCK_Prepare.DOCK_Prepare.changeNumSave(10)
                p = Popen(['python', os.environ["SCRIPTS"] + '/DOCKovalentTask.py', 'CovLib', curr_dir + '/' + ligands[lig], 'True'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
                dock_jobs += [line for line in p.communicate()[0].split('\n') if 'pbs' in line]
                os.chdir(curr_dir)
        Cluster.Cluster.wait(dock_jobs, 1000)
        f.close()
        covlist = glob.glob("PDB/*/CYS*/*/ROT*/CovLib/")
        with open('covlist', 'w') as covf:
                for line in covlist:
                        covf.write(line + '\n')
        cov_jobs = cluster.runJobs('covlist', os.environ["SCRIPTS"] + '/Covalentizer/Final_Scripts/combine.sh')
        Cluster.Cluster.wait(cov_jobs, 1000)
        os.mkdir('Results')
        for c in covlist:
                web_folder = c + '/web_files/'
                new_folder = 'Results/' + c.split('CovLib')[0].replace('/', '_')[:-1] + '/'
                if os.path.isdir(web_folder):
                        shutil.move(web_folder, new_folder)

def print_usage(name):
        print "Usage : " + name + " <pdb_list> <ligands>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
