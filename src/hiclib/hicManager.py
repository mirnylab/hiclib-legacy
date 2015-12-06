from __future__ import absolute_import, division, print_function, unicode_literals
import imp
import numpy as np
import os
import warnings
from mirnylib.h5dict import h5dict
from os.path import join, exists
from os import listdir
import re
from mirnylib.systemutils import setExceptionHook
from hiclib.hicShared import fileIsFragment, fileIsHeatmap
import binascii
from hiclib.fragmentHiC import HiCdataset
from mirnylib.genome import Genome
from mirnylib.numutils import completeIC



setExceptionHook()

def getBases(folder):
    """
    Returns all basefiles from a Hi-C folder
    """
    files = listdir(folder)
    experimentBases = [i[:-13] for i in files if i.endswith("_refined.frag")]
    return experimentBases


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
        self.byChrSuper = {}
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
                checked = fileIsHeatmap(fnameFull)
                if checked == False:
                    print("Filename {0} cannot be loaded", fnameFull)
                    continue
                hmType, resolution = checked
                if hmType == "heatmap":
                    self.heatmaps[resolution] = fnameFull
                if hmType == "byChr":
                    self.byChr[resolution] = fnameFull
                if hmType == "byChrCis":
                    self.byChrCis[resolution] = fnameFull
                if hmType == "byChrSuper":
                    self.byChrSuper[resolution] = fnameFull


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

    def getMinResolution(self):
        return min(list(self.heatmaps.keys()) + list(self.byChr.keys()) + list(self.byChrCis.keys()) + list(self.byChrSuper.keys()))

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
        chunks = list(range(0, mylen, 200000000)) + [mylen]
        chunks = list(zip(chunks[:-1], chunks[1:]))
        c1 = hd.get_dataset("chrms1")
        c2 = hd.get_dataset("chrms2")
        totsum = 0
        for st, end in chunks:
            totsum += np.sum(c1[st:end] == c2[st:end])
        return totsum

    def getEnzyme(self):
        enzymes = ["HindIII", "NcoI", "BglII", "MboI", "DpnII"]
        for enz in enzymes:
            if enz in self.base:
                return enz
        print("Enzyme not found!")
        return None

    def getGenomeObject(self):
        if hasattr(self, "genomeObject"):
            return self.genomeObject
        name = self.genomeName
        base = os.path.split(self.folder)[0]
        defineGenomePath = os.path.join(base, "defineGenome.py")
        assert os.path.exists(defineGenomePath)
        randomString = binascii.b2a_hex(os.urandom(15))
        genomeModule = imp.load_source(randomString, defineGenomePath)
        genomeObject = genomeModule.getGenome(name)
        self.genomeObject = genomeObject
        return genomeObject

    def getScaling(self):
        HD = HiCdataset(self.refined, self.getGenomeObject(), self.getEnzyme(), 1000, mode='r', tmpFolder="\tmp", dictToStoreIDs="h5dict")
        scal = HD.plotScaling(excludeNeighbors=2, normalize=True, mindist=2000)
        return scal

    def getByChromosomeScaling(self):
        HD = HiCdataset(self.refined, self.getGenomeObject(), self.getEnzyme(), 1000, mode='r', tmpFolder="\tmp", dictToStoreIDs="h5dict")
        scals = {}
        for chrom in range(self.getGenomeObject().chrmCount):
            for arm in [0, 1]:
                if arm == 0:
                    region = (chrom, 0, self.genomeObject.cntrMids[chrom])
                else:
                    region = (chrom, self.genomeObject.cntrMids[chrom], self.genomeObject.chrmLens[chrom])
                scal = HD.plotScaling(excludeNeighbors=2, normalize=True, mindist=2000, regions=[region])

                scals[(chrom, arm)] = scal
        return scals



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
        if not fileIsFragment(fname):
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
            for experiment in list(curExperiments.values()):
                experiment.experiment = folder
            globalDict.update(curExperiments)
            print("scanned:", folder, subfolder)
    return globalDict





