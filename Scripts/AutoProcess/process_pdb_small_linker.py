import glob
import subprocess
from subprocess import Popen, PIPE
import sys
import os
import shutil
sys.path.append(os.environ["COVALENTIZER"])
from Code import *

def main(name, argv):
        if len(argv) != 2:
                print_usage(name)
                return

        timeout = 12*60*60
        f = open('log.txt', 'w', 0)
        f.write('Covalentizer has started.\n')

        curr_dir = os.getcwd()
        cluster = Cluster.Cluster()
        f.write('INFO: Looking for cysteine tiol within 6A of any of the ligand atoms.\n')
        res = PYMOLUtils.env_cysteine(argv[0], argv[1], 10)
        with open('res.txt', 'w') as fres:
                for r in res:
                        fres.write("\t".join(r) + '\n')
        if len(res) == 0:
                f.write('WARNING: Did not find cysteines which are close to any of the ligand.\n')
                f.close()
                return

        f.write('INFO: Create folders for the different cysteines and rotamers.\n')
        CysUtils.cysteine_folders(argv[0])
        #dirlist = glob.glob("CYS*/*/ROT*/")
        #with open('dirlist', 'w') as dirf:
        #        for line in dirlist:
        #                dirf.write(line + '\n')
        os.mkdir('Ligands/')
        ligands = {}
        ligands_to_remove = []
        for r in res:
                if r[0] in ligands:
                        continue
                ligands[r[0]] = curr_dir + '/Ligands/' + r[0] + '/gz_files/'
                os.mkdir('Ligands/' + r[0])
                PYMOLUtils.seperate_rec_lig(argv[0], r[0], r[2])
                ChemUtils.convert('xtal-lig.pdb', 'Ligands/' + r[0] + '/' + 'smile.smi')
                mw = CovUtils.get_MW('Ligands/' + r[0] + '/' + 'smile.smi')
                if mw == None or mw < 100:
                        ligands_to_remove.append(r[0])
                os.remove('xtal-lig.pdb')
        if len(ligands_to_remove) > 0:
                f.write('INFO: Removing very small ligands (<300D) and ligand not readable by RDKit.\n')
        #Remove ligands smaller than 300D
        for lig in ligands_to_remove:
                del ligands[lig]
                shutil.rmtree('Ligands/' + lig + '/')
                lig_folders = glob.glob("CYS*/" + lig + "/")
                for l_fol in lig_folders:
                        shutil.rmtree(l_fol)
                #dirlist = [d for d in dirlist if '/' + lig + '/' in d]
                f.write('INFO: Ligand ' + lig + ' was removed.\n')
        if len(ligands) == 0:
                f.write('ERROR: All ligands with nearby cysteines are smaller than 300D. Finish the run.\n')
                f.close()
                return
        #dirlist = glob.glob("CYS*/*/ROT*/")
        #with open('dirlist', 'w') as dirf:
        #        for line in dirlist:
        #                dirf.write(line + '\n')

        f.write('INFO: Build the ligands covalent libraries and prepare them for docking.\n')
        lig_jobs = []
        ligands_to_remove = []
        os.chdir('Ligands')
        for lig in ligands:
                os.chdir(lig)
                CovUtils.build_library('smile.smi', 'frags.smi', 'lib.smi', rules = os.environ["COVALENTIZER"] + "/Code/di_amine_rules.re", linker_lib = True, linker_smiles = os.environ["COVALENTIZER"] + "/Code/di_amine_linkers.smi")
                if os.stat('frags.smi').st_size == 0 or os.stat('smile.smi').st_size == 0 or os.stat('lib.smi').st_size == 0:
                        ligands_to_remove.append(lig)
                        os.chdir('../')
                        continue
                p = Popen([os.environ["DOCKBASE"] + '/ligand/generate/build_covalent_lib_medium.csh', 'lib.smi', '100', 'LIB'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
                lig_jobs += [line for line in p.communicate()[0].split('\n') if 'pbs' in line]
                os.chdir('../')
        os.chdir('../')
        if len(ligands_to_remove) > 0:
                f.write('INFO: Removing ligands which were not fragmentable or preparable by RDKit.\n')
        #Remove ligands which were not fragmentable or preparable by RDKit
        for lig in ligands_to_remove:
                del ligands[lig]
                shutil.rmtree('Ligands/' + lig + '/')
                lig_folders = glob.glob("CYS*/" + lig + "/")
                for l_fol in lig_folders:
                        shutil.rmtree(l_fol)
                f.write('INFO: Ligand ' + lig + ' was removed.\n')
        if len(ligands) == 0:
                f.write('ERROR: All ligands with nearby cysteines are smaller than 300D or unpreparable by RDKit. Finish the run.\n')
                f.close()
                return

        dirlist = glob.glob("CYS*/*/ROT*/")
        with open('dirlist', 'w') as dirf:
                for line in dirlist:
                        dirf.write(line + '\n')

        f.write('INFO: Job IDs for ligand building: ' + str(lig_jobs) + '\n')
        f.write('INFO: Start preparing the structures for docking.\n')
        prepare_jobs = cluster.runBatchJobs('dirlist', './prepare.sh', mem='16000mb')
        f.write('INFO: Job IDs for structure prepare: ' + str(prepare_jobs) + '\n')
        waiting_done = Cluster.Cluster.wait(prepare_jobs + lig_jobs, timeout)
        if waiting_done:
                f.write('INFO: Finished preparing both proteins and molecules for docking.\n')
        else:
                f.write('ERROR: Cluster jobs have reached time limit (12h). That could be due to a cluster occupancy, problematic input, or any other unidentified problem. To rule out the first problem, try submitting the job at some other time.\n')
                f.close()
                return
        LIB_folders = glob.glob('Ligands/*/LIB.*')
        for lib in LIB_folders:
                shutil.rmtree(lib)
        f.write('INFO: Start docking processes for all covalentized fragments to all cysteine rotamers.\n')
        dock_jobs = []
        for d in dirlist:
                lig = d.split('/')[1]
                os.chdir(d)
                DOCK_Prepare.DOCK_Prepare.changeNumSave(10)
                p = Popen(['python', os.environ["COVALENTIZER"] + 'Scripts/DOCKovalentTask.py', 'CovLib', ligands[lig], 'True'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
                dock_jobs += [line for line in p.communicate()[0].split('\n') if 'pbs' in line]
                os.chdir(curr_dir)
        f.write('INFO: Job IDs for covalent docking: ' + str(dock_jobs) + '\n')
        waiting_done = Cluster.Cluster.wait(dock_jobs, timeout)
        if waiting_done:
                f.write('INFO: Finished preparing both proteins and molecules for docking.\n')
        else:   
                f.write('ERROR: Cluster jobs have reached time limit (12h). That could be due to a cluster occupancy, problematic input, or any other unidentified problem. To rule out the first problem, try submitting the job at some other time.\n')
                f.close()
                return
        f.write('INFO: Finished docking.\n')
        f.write('INFO: Start analyzing the docking results.\n')
        ##
        reslist = glob.glob("CYS*/*/ROT*/CovLib/")
        with open('reslist', 'w') as dirf:
                for line in reslist:
                        dirf.write(line + '\n')
        results_jobs = cluster.runBatchJobs('reslist', os.environ["COVALENTIZER"] + 'Scripts/Final_Scripts/combine_cluster.sh', mem='32000mb')
        f.write('INFO: Job IDs for analyzing the results: ' + str(results_jobs) + '\n')
        waiting_done = Cluster.Cluster.wait(results_jobs, timeout)
        if waiting_done:
                f.write('INFO: Finished analyzing docking results.\n')
        else:
                f.write('ERROR: Cluster jobs have reached time limit (12h). That could be due to a cluster occupancy, problematic input, or any other unidentified problem. To rule out the first problem, try submitting the job at some other time.\n')
                f.close()
                return
        ##
        #os.system(os.environ["SCRIPTS"] + '/Covalentizer/Final_Scripts/combine.sh')
        os.mkdir('Results')
        are_results = False
        for d in dirlist:
                web_folder = d + '/CovLib/web_files/'
                if os.path.isdir(web_folder):
                        with open(web_folder + 'cys_position.txt', 'r') as cys_file:
                                cys = cys_file.readline().split()[0]
                        are_results = True
                        name_list = d.split('/')[:3]
                        name_list[0] = 'CYS' + cys
                        new_folder = 'Results/' + '_'.join(name_list) + '/'
                        #new_folder = 'Results/' + d.replace('/', '_')[:-1] + '/' 
                        shutil.move(web_folder, new_folder)
        if are_results:
                f.write('INFO: Moved results into the Results folder.\n')
        else:
                f.write('INFO: No results have been found with RMSD below 1.5A\n')
        f.write('INFO: Covalentizer has finished.\n')
        f.close()

def print_usage(name):
        print "Usage : " + name + " <pdb_file> <lig_name>"

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
