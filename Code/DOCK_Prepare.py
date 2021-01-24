#Written by Daniel Zaidman
#Code review by 
import shutil
import subprocess
import os, sys
import PyUtils

class DOCK_Prepare:
    def __init__(self, rec, lig, cov, cov_index, HG, tart=False):
        self.tart = bool(tart)
        self.rec = rec
        self.lig = lig
        self.folder = os.path.dirname(os.path.abspath(rec)) + "/"
        self.fixed_rec = self.folder + "rec.pdb"
        self.fixed_lig = self.folder + "xtal-lig.pdb"
        if os.path.isdir(self.folder + "working/"):
            shutil.rmtree(self.folder +"working/")
        if os.path.isdir(self.folder + "dockfiles/"):
            shutil.rmtree(self.folder +"dockfiles/")
        self.cov = cov
        self.cov_index = cov_index
        self.hg = HG
        PyUtils.initPythonVs()
    def blaster(self):
        self.create_fixed_names()
        if not self.tart:
            subprocess.call([os.environ['DOCKBASE'] + "/proteins/blastermaster/blastermaster.py", "--covalentResNum", self.cov_index, "--covalentResName", self.cov, "--covalentResAtoms", self.hg])
        else:
            print 'tarting'
            subprocess.call([os.environ['DOCKBASE'] + "/proteins/blastermaster/blastermaster.py", "--covalentResNum", self.cov_index, "--covalentResName", self.cov, "--covalentResAtoms", self.hg, "--chargeFile=" + self.folder + "amb.crg.oxt", "--vdwprottable=" + self.folder + "prot.table.ambcrg.ambH"])
        
    def changeIndock(self):
        INDOCK = self.folder + "INDOCK"
        old = open(INDOCK, 'r')
        new = open(INDOCK + '2', 'w')
        for i in range(17):
            line = old.readline()
            new.write(line)
        for i in range(2):
            line = old.readline()
            new.write(line[:-3] + "00.0\n")
        for i in range(24):
            line = old.readline()
            new.write(line)
        line = old.readline()
        new.write(line[:-3] + "yes\n")
        for line in old:
            new.write(line)
        old.close()
        new.close()
        os.remove(INDOCK)
        os.rename(INDOCK + '2', INDOCK)
#Inner functions
    def create_fixed_names(self):
        if(not os.path.exists(self.fixed_rec)):
            shutil.copyfile(self.rec, self.fixed_rec)
        if(not os.path.exists(self.fixed_lig)):
            shutil.copyfile(self.lig, self.fixed_lig)
    @staticmethod
    def changeNumSave(num_save):
        INDOCK = "INDOCK"
        old = open(INDOCK, 'r')
        new = open(INDOCK + '2', 'w')
        for i in range(27):
            line = old.readline()
            new.write(line)
        for i in range(2):
            line = old.readline()
            new.write(line[:-2] + str(num_save) + "\n")
        for i in range(15):
            line = old.readline()
            new.write(line)
        for line in old:
            new.write(line)
        old.close()
        new.close()
        os.remove(INDOCK)
        os.rename(INDOCK + '2', INDOCK)
