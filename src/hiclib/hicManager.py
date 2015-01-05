import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
from mirnylib.h5dict import h5dict
from os.path import join, exists
from os import listdir
import re
from mirnylib.systemutils import setExceptionHook
import shutil
from mirnylib.genome import Genome
from mirnylib.numutils import ultracorrect, adaptiveSmoothing, removeDiagonals, completeIC
from hiclib.fragmentHiC import HiCdataset
from hiclib import binnedData
import subprocess
from mirnylib.numutils_new import observedOverExpected


setExceptionHook()

def getBases(folder):
    """
    Returns all basefiles from a Hi-C folder
    """
    files = listdir(folder)
    experimentBases = [i[:-13] for i in files if i.endswith("_refined.frag")]
    return experimentBases


def checkFragment(filename):
    "Checks whether a file is a fragment-lefel h5dict"
    try:
        a = h5dict(filename, "r")
    except:
        return False
    if "chrms1" in a.keys():
        return True
    return False


def getResolution(fname):
    "Extracts resolution from a Hi-C filename"
    matches = re.findall(r"-\d{1,5}[k,M]", fname)
    if len(matches) == 0:
        raise ValueError("resolution not found")
    if len(matches) > 1:
        warnings.warn("undetermined resolution")
        print "found resolutions", matches
        print "using the last found"
    match = matches[-1]
    if match[-1] == "k":
        mult = 1000
    elif match[-1] == "M":
        mult = 1000000
    else:
        raise ValueError("Bad things happened")

    num = int(match[1:-1])
    return num * mult


def checkHeatmap(filename):
    """Checks whether a file is a heatmap
    Returns a filename as a first argument, and resolution as a second"""
    try:
        a = h5dict(filename, "r")
    except:
        return False

    keys = a.keys()
    if ("heatmap" in keys) and ("resolution" in keys):
        resolution = getResolution(filename)
        assert resolution == a["resolution"]
        return ("heatmap", a["resolution"])
    elif "heatmap" in keys:
        warnings.warn("Resolution not found in {0}".format(filename))
        resolution = getResolution(filename)
        return ("heatmap", resolution)

    if "0 0" in keys:
        resolution = getResolution(filename)
        if "0 1" in keys:
            return "byChr", resolution
        else:
            return "byChrCis", resolution
    return False




class hicExperiment(object):
    def __init__(self, fol, base, fnames, genomeName):
        """
        A class managing a Hi-C experiment.

        Parameters
        ----------

        fol : str
            Folder containing the files
        base : str
            Prefix of filenames of the experiment
        fnames : str
            Filenames with experiment files (usually start with the same prefix as base)
        genomeName : str
            Name of the genome

        """
        self.folder = fol
        self.base = base
        if not exists(join(fol, base + "_refined.frag")):
            raise ValueError("Refined data do not exists")
        self.refined = join(fol, base + "_refined.frag")
        self.heatmaps = {}
        self.byChr = {}
        self.byChrCis = {}
        self.merged = None
        self.genomeName = genomeName

        isMerged = [i.endswith("_merged.frag") for i in fnames]
        if  sum(isMerged) > 1:
            raise ValueError("Multiple merged files found")

        for fname in fnames:
            if fname.endswith("_merged.frag"):
                self.merged = join(fol, fname)
            elif fname.endswith("_refined.frag"):
                pass
            else:
                fnameFull = join(fol, fname)
                checked = checkHeatmap(fnameFull)
                if checked == False:
                    print "Filename {0} cannot be loaded", fnameFull
                    continue
                hmType, resolution = checked
                if hmType == "heatmap":
                    self.heatmaps[resolution] = fnameFull
                if hmType == "byChr":
                    self.byChr[resolution] = fnameFull
                if hmType == "byChrCis":
                    self.byChrCis[resolution] = fnameFull

    def isReplica(self):
        """
        Returns True if this is file is marked as a replica,
        (and not a combination of multiple replicas).
        Replicas are marked as -R1, -R2, -Rblabla
        """
        found = re.findall(r"-R\d", self.base)
        if len(found) == 1:
            return True
        elif len(found) == 2:
            warnings.warn("More than one match to replica code found")
        return False

    def isMaster(self):
        """
        Returns False if there is a combined experiment which includes this one and other replicas.
        Returns true if this is a combined experiment, or the only replica.
        """
        if not self.isReplica():
            return True
        found = re.findall(r"-all", self.base)
        if len(found) == 1:
            return True
        proposedBase = re.sub(r"-R\d.*-", "-all-", self.base)
        bases = getBases(self.folder)

        if proposedBase in bases:
            return False
        return True

    def getNumReads(self):
        hd = h5dict(self.refined, 'r')
        return len(hd.get_dataset("strands1"))

    def getNumCisReads(self):
        hd = h5dict(self.refined, 'r')
        mylen = len(hd.get_dataset("strands1"))
        chunks = range(0,mylen,200000000) +  mylen
        chunks = zip(chunks[:-1],chunks[1:])
        c1 = hd.get_dataset("chrms1")
        c2 = hd.get_dataset("chrms2")
        totsum = 0 
        for st,end in chunks:
            totsum += np.sum(c1[st:end] == c2[st:end])
        return totsum

    def getEnzyme(self):
        enzymes = ["HindIII", "NcoI", "BglII", "MboI"]
        for enz in enzymes:
            if enz in self.base:
                return enz
        print "Enzyme not found!"
        return None



def scanHicFolder(foldername, genomeName):
    """
    Scans folder with multiple experiments and extracts experiments
    """

    fol = foldername
    if not exists(foldername):
        raise ValueError("Folder {0} does not exist".format(foldername))
    files = listdir(foldername)
    expNames = getBases(fol)
    for i in expNames:
        fname = join(fol, i + "_refined.frag")
        if not checkFragment(fname):
            raise

    expDict = {i :[j for j in files if j.startswith(i)] for i in expNames}

    experiments = {}

    for exp in expDict:
        experiment = hicExperiment(fol, base=exp, fnames=expDict[exp], genomeName=genomeName)
        # print experiment.heatmaps
        experiments[exp + "-" + genomeName] = experiment

    return experiments


def scanHiCFolders(folderList, genomes=["hg18", "hg19", "mm9", "mm10"], baseFolder=None):
    """
    Scan folder with multiple folders with Hi-C experiments, according to the current scheme.
    (.../HiCFolder/hg19/HiCExperiment-R1-HindIII-100k.hm)
    """
    if baseFolder != None:
        folderList = [join(baseFolder, i) for i in folderList]
    globalDict = {}

    for folder in folderList:
        subfolders = os.listdir(folder)
        subfolders = [i for i in subfolders if i in genomes]
        for subfolder in subfolders:
            curExperiments = scanHicFolder(join(folder, subfolder), genomeName=subfolder)
            for experiment in curExperiments.values():
                experiment.experiment = folder
            globalDict.update(curExperiments)
            print "scanned:", folder, subfolder
    return globalDict





