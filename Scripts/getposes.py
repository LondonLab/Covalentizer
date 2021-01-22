#!/usr/bin/env python
"""for each input ID, find the X best poses from docking
Ryan Coleman 2011-2012.
Teague Sterling 2015: Now uses ZINC REST interface
"""

import os
import sys
sys.path.append('/work/londonlab/git_dock/DOCK/analysis/')
import string
import optparse 
import extract_all
import operator
import gzip
import bsddb
import marshal
import tempfile #have to use tempfiles for cache files to prevent deadlock
import shutil

import zincapi

ZINC_FIELD_MAP = {
    'smiles': 'Substance Smiles',
    'vendor_catalogs': 'Vendor',
    'purchasability': 'Purchasability',
    'reactivity': 'Reactivity',
    'tranche': 'Tranche Name',
}


def getPosesCache(file, molnum, id, rank, cache=True):
  '''reads from the bsddb cache if possible, otherwise make the cache first'''
  tempDir, tempFile = os.path.split(file)
  bsFile = file + ".db"
  okay = "OKAY_LIGAND" #used to detect correct bsddb file state
  needToCreate = True
  #print bsFile, os.path.exists(bsFile), os.path.getsize(bsFile)
  if os.path.exists(bsFile):
    if os.path.getsize(bsFile) > 0:
      cacheDict = bsddb.hashopen(bsFile, 'r')
      #if there is no okay key, or if it isn't True
      #print cacheDict[okay], marshal.loads(cacheDict[okay]), (not cacheDict.has_key(okay)) or (not marshal.loads(cacheDict[okay]))
      if (not cacheDict.has_key(okay)) or (not marshal.loads(cacheDict[okay])):
        #have to clear the database
        cacheDict.close()
      else:
        needToCreate = False
      ##(molnum, id, rank) tuple
  if needToCreate: #construct cache
    print "creating cache for " + file
    bsFileTemp = tempfile.NamedTemporaryFile(suffix='.db', prefix=tempFile, \
        dir=tempDir, delete=False)
    bsFileTemp.close()
    cacheDict = bsddb.hashopen(bsFileTemp.name)
    cacheDict[okay] = marshal.dumps(False)
    ligandFile = gzip.open(file, 'rb')
    ligString = []
    mirTuple = ['-', '-', '-'] #molnum, id, rank tuple
    try:
      for line in ligandFile:
        if line.find("  Name:") > 0: #start of new record here
          if len(ligString) > 0:
            #print mirTuple
            prevValue = []
            if cacheDict.has_key(str(mirTuple)):
              prevValue = marshal.loads(cacheDict[str(mirTuple)])
            prevValue.append(ligString)
            cacheDict[str(mirTuple)] = marshal.dumps(prevValue)
          zincId = string.strip(string.split(line)[2])
          mirTuple = ['-', zincId, '-'] #molnum, id, rank  tuple
          ligString = [line] #reset ligstring here
        else:
          ligString.append(line)
        if line.find("Number:") > 0:
          mNum = int(string.strip(string.split(line)[2]))
          mirTuple[0] = mNum
        if line.find("Rank:") > 0:
          thisRank = int(string.strip(string.split(line)[2]))
          mirTuple[2] = thisRank
    except StopIteration:
      pass #EOF is fine
    ligandFile.close()
    if len(ligString) > 0:
      prevValue = []
      if cacheDict.has_key(str(mirTuple)):
        prevValue = marshal.loads(cacheDict[str(mirTuple)])
      prevValue.append(ligString)
      cacheDict[str(mirTuple)] = marshal.dumps(prevValue)
      #print mirTuple, "added"
    cacheDict[okay] = marshal.dumps(True) #set to true, uninterrupted run
    cacheDict.sync()
  #print [molnum, id, rank], "retrieved"
  returnItem = marshal.loads(cacheDict[str([molnum, id, rank])])
  cacheDict.close()
  if needToCreate:
    if not cache: #delete the cache to save disk space
      os.unlink(bsFileTemp.name)
    else:
      #print 'copying ', bsFileTemp.name, bsFile
      shutil.copy(bsFileTemp.name, bsFile)
      os.unlink(bsFileTemp.name)
  return returnItem

def getOnePose(file, molnum, id, rank, recList=None, realRank=None, part=None, \
        replaceScore=None, replaceRecScore=None, cache=True):
  '''in one file, get the mol2 record for the (mol#, id, rank) combo.
  if idlist is present, only save molecules if they are in that id.
  if recList is present, only save molecules if they match one receptor.
  if part is not None, only return flexible receptors containing that part.
  if replaceScore is not None, it means replace the score here 
    since it was rescored, same with replaceRecScore'''
  cachedPoses = getPosesCache(file, molnum, id, rank, cache)
  retString = []
  for cachedPose in cachedPoses: #check if this one should be returned
    thisRec = None
    for line in cachedPose:
      if realRank is not None:
        if line.find("Rank:") > 0:
          #also put the global unique rank here
          retString.append("##########           GlobalRank: " + \
              str(realRank) + "\n")
      if replaceScore is not None and line.find("Total Energy") > 0:
        line = "##########         Total Energy: " + str(replaceScore) + "\n"
      if replaceRecScore is not None and line.find("Receptor Energy") > 0:
        line = "##########      Receptor Energy: " + str(replaceRecScore) + "\n"
      if line.find("FlexRecCode") > 0:
        thisRec = string.strip(string.split(line)[2])
      retString.append(line)
    if recList is None or thisRec in recList:
      if part is None or (part in string.split(thisRec, '.')):
        return retString

def getPoses(indir, dockdir, output, idfile, topX, receptors, rankfile, part, \
        limit, zincData, dockMol2name='test.mol2.gz', idToScores=None, \
        quiet=False, cache=True, runPartPrefs=True, bbrw=True, protrmsd=True, \
        bbndist=1, bbndof=1):
  """for each id in idlist, find the best X poses, write them to output.
  if receptors is not an empty list, only get poses to those receptors in list
  if rankfile present, output a rank into the output file based on line number.
  if part is not None, only get poses with that flexible part.
  limit limits the number of ids to use
  zincData controls connecting to zinc to get smiles/vendor information
  idToScores allows passing in of the data read somewhere else, so reuse it
  """
  if zincData:
    zincData = zincapi.ZINCAPI()
    if not quiet:
      print "Attempting to connect to ZINC REST interface"
      if not quiet:
        print "Connected to ZINC database, will output Vendor information"
  if 0 == len(receptors): #means no receptors were given, so do them all.
    receptors = None #later functions expect None to mean do them all.
  idlist = extract_all.readIds(idfile, receptors, part) #read ids
  if limit > 0: #0 = keep them all
    idlist = idlist[:limit]
  cacheRanked = None
  if idToScores is None: #not passed in, have to read it
    idToScores = extract_all.getId2Scores(indir, idlist, receptors, part)
  outputFile = open(output, 'w')
  if not quiet:
    if topX > 0:
      print "getting up to ", topX, " poses for each of ", len(idlist), " ligands"
    else:
      print "getting ALL poses for each of ", len(idlist), " ligands"
  for anId in idlist:
    if not quiet:
      print "gathering poses for molecule", anId
    scoreList = idToScores[str(anId)]
    scoreList.sort(key=operator.itemgetter(extract_all.scoreCol))
    #print scoreList, type(anId), idToScores.keys()
    if 0 == topX: #0 means get them all
      topXfind = len(scoreList)
    else:
      topXfind = topX
    #for black box scoring, we need poses as well as scores. so get them.
    scorePose = {}
    if len(scoreList) > 0:
      howMany = min(topXfind, len(scoreList))
      if bbrw or runPartPrefs:
        howMany = len(scoreList) #do them all for black box reweighting
      for position in xrange(howMany):
        dirName = str(scoreList[position][extract_all.dirCol])
        mol2name = os.path.join(dockdir, dirName, dockMol2name)
        molnum = scoreList[position][extract_all.molnumCol]
        rank = scoreList[position][extract_all.rankCol]
        #scoreToCheck = float(scoreList[position][-1])
        scoreToCheck = float(scoreList[position][extract_all.scoreCol])
        newRecScore = float(scoreList[position][extract_all.recScoreCol])
        heavyAtomCount = int(scoreList[position][extract_all.heavyAtomCol])
        realRank, cacheRanked = extract_all.rankScores(scoreToCheck, indir, \
              rankfile, cacheRanked=cacheRanked)
        thisMol2 = getOnePose(mol2name, molnum, anId, rank, \
              receptors, realRank, part, scoreToCheck, newRecScore, cache) 
        scorePose[tuple(scoreList[position])] = thisMol2
    allMols = extract_all.makeAllMol2(scoreList, scorePose)
    '''clusters = allMols.breakIntoClustersByAtomCount()
    whichCluster = None
    otherClusters = []
    for cluster in clusters:
      if 0 in cluster: #find the first (best score) ligand in the cluster
        whichCluster = cluster
      else:
        otherClusters.append(cluster)
    #want to keep scorelist only for the cluster with the best scoring 
    #protonation state
    newScoreList = []
    otherClusterLists = []
    for scoreListCount, scoreListItem in enumerate(scoreList):
      if scoreListCount in whichCluster:
        newScoreList.append(scoreListItem)
      for otherCount, otherCluster in enumerate(otherClusters):
        otherClusterLists.append([])
        if scoreListCount in otherCluster:
          otherClusterLists[otherCount].append(scoreListItem)
    scoreList = newScoreList'''
    if runPartPrefs:
      partPrefs = extract_all.partPreferences(newScoreList)
    if bbrw:
      partPrefsBlackBox, ndist, dof = \
          extract_all.partPreferencesBlackBox(scoreList, scorePose, \
          protrmsd=protrmsd, ndist=bbndist, dof=bbndof)
    #infoContent = extract_all.informationContent(partPrefs.itervalues())
    if len(scoreList) > 0:
      for position in xrange(min(topXfind, len(scoreList))): #only do the topX 
        dirName = scoreList[position][extract_all.dirCol]
        mol2name = os.path.join(dockdir, dirName, dockMol2name)
        molnum = scoreList[position][extract_all.molnumCol]
        rank = scoreList[position][extract_all.rankCol]
        scoreToCheck = float(scoreList[position][extract_all.scoreCol])
        newRecScore = float(scoreList[position][extract_all.recScoreCol])
        heavyAtomCount = int(scoreList[position][extract_all.heavyAtomCol])
        realRank, cacheRanked = extract_all.rankScores(scoreToCheck, indir, \
              rankfile, cacheRanked=cacheRanked)
        thisMol2 = scorePose[tuple(scoreList[position])]
        #use zinc, rank, number & receptors not energy, this is way saner.
        for oneLine in thisMol2:
          if oneLine.find("##########        Ligand Energy:") != -1:
            if runPartPrefs:
              partKeys = partPrefs.keys()
              partKeys.sort()
              for partKey in partKeys:
                outputFile.write("##########            part" + str(partKey))
                outputFile.write("pref: %3.2f\n" % partPrefs[partKey])
            if bbrw:
              outputFile.write("##########              BBndist: " + \
                  str(ndist) + '\n')
              outputFile.write("##########                BBdof: " + \
                  str(dof) + '\n')
              partKeysBlackBox = partPrefsBlackBox.keys()
              partKeysBlackBox.sort()
              for partKey in partKeysBlackBox:
                outputFile.write("##########         BBpart" + str(partKey))
                outputFile.write("_pref: %3.2f\n" % partPrefsBlackBox[partKey])
            #always write heavy atom count
            outputFile.write('##########     heavy atom count: %4d\n' % \
                heavyAtomCount)
            outputFile.write('##########           pose count: %4d\n' % \
                len(scoreList))
          if oneLine.find("@<TRIPOS>MOLECULE") != -1:
            #write part/pref data on one line for now
            if zincData: #before the tripos molecule line, write vendor info
              zincSubstance = zincData.substances.get(anId, fields=ZINC_FIELD_MAP.keys())
              for key, field in ZINC_FIELD_MAP.items():
                value = zincSubstance[key]
                if isinstance(value, dict):
                  for vKey, vValue in value.items():
                    if not isinstance(vValue, list):
                      vValue = [vValue]
                    zincLine = "##########{}: {!s} {!s}\n".format(field.rjust(22), vKey, ' '.join(vValue))
                    outputFile.write(zincLine)
                elif isinstance(value, list):
                  for vValue in value:
                    zincLine = "##########{}: {!s}\n".format(field.rjust(22), vValue)
                    outputFile.write(zincLine)
                else:
                    zincLine = "##########{}: {!s}\n".format(field.rjust(22), value)
                    outputFile.write(zincLine)
                    
          outputFile.write(oneLine)
    else:
      if not quiet:
        print "nothing found for ligand ", anId
  outputFile.close()  
  if zincData:
    dock_zinc.close_zinc(zincDb, zincCursor)

def main(argv):
  description = "Get the best X poses of Y input codes."
  version = "%prog *version 201204* "
  usage = "%prog [options] "
  parser = optparse.OptionParser(usage=usage, description=description, \
      version=version)
  parser.add_option("-i", "--indir", type="string", action="store", \
      dest="indir", default=".", \
      help="input directory (for extract* files), (default: %default)")
  parser.add_option("-d", "--dockdir", type="string", action="store", \
      dest="dockdir", default=".", \
      help="input directory (for */test.mol2.gz files), (default: %default)")
  parser.add_option("-o", "--output", type="string", action="store", \
      dest="output", default="./poses.mol2", \
      help="output filename, (default: %default)")
  parser.add_option("-f", "--file", type="string", action="store", \
      dest="idfile", default=extract_all.uniqFileName, \
      help="file containing ids to retrieve (default: %default)")
  parser.add_option("-x", "--topx", type="long", action="store", \
      dest="topx", default=1, \
      help="how many poses of each compound, 0 = all, (default: %default)")
  parser.add_option("-l", "--toplimit", type="long", action="store", \
      dest="toplimit", default=500, \
      help="how many top compounds to save, 0 = all, (default: %default)")
  parser.add_option("-r", "--receptor", type="string", action="append", \
      dest="receptor", default=[], \
      help="which receptor(s) to get poses for, specify more with additional -r commands (default: All)")
  parser.add_option("-p", "--part", type="string", action="store", \
      dest="part", default=None, \
      help="if specified, will only use data from receptors with that part (default: All)")
  parser.add_option("--ranks", type="string", action="store", \
      dest="rankfile", default=extract_all.uniqFileName, \
      help="add ranks for each ZINC code, from file (default: %default)")
  parser.add_option("--mol2name", type="string", action="store", \
      dest="dockMol2name", default='test.mol2.gz', \
      help="DOCK mol2 output filename (default: %default)")
  parser.add_option("-z", "--zinc-data", action="store_true", default=False, \
      dest="zincData", \
      help="attempt to connect to ZINC and get SMILES/VENDOR information (default: %default)")
  parser.add_option("-q", "--quiet", action="store_true", default=False, \
      dest="quiet", \
      help="less standard out remarks desired (default: lots of comments)")
  parser.add_option("-n", "--nocache", action="store_false", default=True, \
      dest="cache", \
      help="don't make cache files for faster getposes.py next time (default: lots of cache files)")
  parser.add_option("--enablepp", action="store_true", default=False, \
      dest="runPartPrefs", \
      help="run part preference calc. (default: don't run it)")
  parser.add_option("--enablebb", action="store_true", default=False, \
      dest="bbrw", \
      help="run black box reweighting (default: don't run the black box reweighting)")
  parser.add_option("--noprotrmsd", action="store_false", default=True,
      dest='protrmsd', \
      help="for bbrw, don't use protein as well as ligand RMSD. needs 1.x.x.x.pdb files. (default: use protein and ligand RMSD)")
  parser.add_option("--bbndist", action="store", type="long", default=10, \
      dest='bbndist', \
      help="for bbrw, the number of nearby neighbors to use (ndist) (default: %default)")
  parser.add_option("--bbndof", action="store", type="long", default=1, \
      dest='bbndof', \
      help="for bbrw, the degrees of freedom/dimensionality (dof) (default: %default)")
  options, args = parser.parse_args(args=argv[1:])
  if 0 != len(args):
    parser.error("program takes no positional arguments.\n" +
                 "  Use --help for more information.")
  else:
    getPoses(options.indir, options.dockdir, options.output, options.idfile, \
        options.topx, options.receptor, rankfile=options.rankfile, \
        part=options.part, limit=options.toplimit, zincData=options.zincData, \
        dockMol2name=options.dockMol2name, quiet=options.quiet, \
        cache=options.cache, runPartPrefs=options.runPartPrefs, \
        bbrw=options.bbrw, protrmsd=options.protrmsd, \
        bbndist=options.bbndist, bbndof=options.bbndof)

if -1 != string.find(sys.argv[0], "getposes.py"):
  main(sys.argv)
