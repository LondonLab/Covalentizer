#Written by Daniel Zaidman
#Code review by 

import shutil
import subprocess
import os
import PyUtils
import Cluster

class DOCKovalent:
    def __init__(self, folder_name, compound, library = False):
	self.folder = os.getcwd() + "/"
        self.name = self.folder + folder_name + "/"
        #self.compound = os.path.abspath(compound)[5:]
        self.compound = compound
        self.library = library
        PyUtils.create_folder(self.name)
        self.copyIndock()
        self.softlink()
        os.chdir(self.name)
        self.dock_command = os.environ['DOCKBASE'] + "/docking/DOCK/src/i386/dock64"
        self.DOCK()
    #Inner functions
    def copyIndock(self):
        INDOCK_old = self.folder + "INDOCK"
        INDOCK = self.name + "INDOCK"
        shutil.copy(INDOCK_old, INDOCK)
        if(not self.library):
            old = open(INDOCK_old, 'r')
            new = open(INDOCK, 'w')
            for i in range(3):
                line = old.readline()
                new.write(line)
            line = old.readline()
            new.write(line[:-21] + self.compound + "\n")
            for line in old:
                new.write(line)
            old.close()
            new.close()
        else:
            shutil.copy(INDOCK_old, INDOCK)
    def softlink(self):
        files = "dockfiles"
        PyUtils.create_softlink(self.folder + files, self.name + files)
    def DOCK(self):
        if(self.library):
            subprocess.call([os.environ['DOCKBASE'] + "/docking/setup/setup_db2.csh", self.compound])
            clu = Cluster.Cluster()
            #for line in clu.runJobs("dirlist", self.dock_command):
            for line in clu.runBatchJobs("dirlist", self.dock_command):
                print line
        else:
            subprocess.call([self.dock_command])
    def combineResults(self):
        subprocess.call([os.environ['DOCKBASE'] + "/analysis/extract_all.py", "--done"])
        subprocess.call([os.environ['DOCKBASE'] + "/analysis/getposes.py"])
