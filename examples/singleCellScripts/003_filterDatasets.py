"""
This example script does the following:

-Loads individual files from the location defined in datasets.tsv file
-Parses individual files in memory  (inMemory = True in "for onename in in_files))
   -If your system does not have enough memory, you might need to switch to hdf5 here.
-Merges files corresponding to the same experiment together, on the HDD.
-Filters datasets, builds heatmaps
-Combines multiple replicas of the same experiment together, builds heatmaps

--Datasets are defined in the datasets.tsv file
--genome is defined by genomeFolder function, and workingGenome identifyer
--output files are arranged to folders named by their workingGenome IDs

Warnings:
    Running this over NFS might cause unexpected slow-downs because NFS is
    unhappy with repeated read/write cycles to the same file

    You could do entire thing in memory, if you have RAM or your datasets are small.
    Actually, using HDF5 is then equivalent to storing compressed data in RAM,
    and might be in fact pretty fast.

    General recommendation: if you have 16+GB of RAM, and .sra (.fastq.gz) files were less than 30GB, then you should be fine with parsing things in memory.

"""
from multiprocessing import Pool
from hiclib.fragmentHiC import HiCdataset as HiCdatasetorig
from mirnylib.systemutils import fmap, setExceptionHook
import numpy as np
import os
import sys
from mirnylib.numutils import uniqueIndex
import pyximport; pyximport.install()


setExceptionHook()

def ensure(f):
    "creates directory for the file if doesn't exist, and returns filename"
    d = os.path.dirname(f)
    if os.path.isdir(d):
        return f
    else:
        try:
            os.makedirs(d)
        except:
            raise ValueError("Cannot create directory")
    return f

def genomeFolder(name):
    return os.path.join("/home/magus/HiC2011/data", name)  # Fetch genome folder by genome name

class HiCdataset(HiCdatasetorig):
    "Modification of HiCDataset to include all filters"
    def filterLessThanDistance(self):
        # This is the old function used to filter "duplicates".
        #After the final submission of the manuscript, It was replaced by a better function that does the same,
        #but at bp resolution, not 100 bp.
        M = self.N
        for i in range(5):
            for j in range(5):
                chrStrandID = 10000000 * 10000000 * (np.array(self.chrms1 * (self.strands1 + 1), dtype = np.int64) * 100 + self.chrms2 * (self.strands2 + 1))
                print(len(np.unique(chrStrandID)))
                posid = np.array((self.cuts1 + i * 100) // 500, dtype = np.int64) * 10000000 + (self.cuts2 + i * 100) // 500
                N = self.N
                self.maskFilter(uniqueIndex(posid + chrStrandID))
                print(N, "filtered to", self.N)
        self.metadata["321_quasiDuplicatesRemoved"] = M - self.N

    def filterLessThanDistanceUpgraded(self, distance=500):
        # This is a faster and more efficient way to remove reads that start/end within 500bp
        # It is written in Cython, so it operates at 1bp resolution.
        import maskRemover
        c1 = np.array(self.chrms1, dtype = np.int16)
        c2 = np.array(self.chrms2, dtype = np.int16)
        p1 = np.array(self.cuts1, dtype = np.int32)
        p2 = np.array(self.cuts2, dtype = np.int32)
        s1 = np.array(self.strands1, dtype = np.int8)
        s2 = np.array(self.strands2, dtype = np.int8)
        removeMask = maskRemover.duplicateRemoveMask(c1, c2, p1, p2, s1, s2, offset=distance, method="max")
        M = self.N
        self.maskFilter(removeMask == 0)
        self.metadata["321_quasiDuplicatesRemoved"] = M - self.N



coolResolutions = [10000000, 5000000, 2000000, 1000000, 500000, 200000, 100000, 40000, 20000, 10000, 5000, 2000, 1000]
skip = 1  # how many to skip for single replica datasets



def refineDataset(filenames):
    """
    Parameters
    ----------

    filenames[0] is a list of filenames of incoming files
    filenames[1] is a folder for outgoing file
    filenames[2] is a working genome, that is output directory
    filenames[3] is an enzyme for a given experiment


    create : bool, optional
        If True, parse each file.
        If False, assume that files were already parsed
        (e.g. if you are just playing around with filtering parameters)
    delete : bool, optional
        If True, delete parsed files after merging.
        Man, these files may be huge... if you don't have a 10TB RAID, this may be useful.
    parseInMemory : bool, optional
        Perform parsing input files in memory.

    """
    in_files = filenames[0]
    out_file = filenames[1]


    workingGenome = filenames[2]
    enzyme = 500

    def parseMappedFile(onename):
        np.random.seed()

        # create dataset in memory, parse and then save to destination
        TR = HiCdataset("bla" + str(np.random.randint(100000000000)), genome=genomeFolder(workingGenome),
                        maximumMoleculeLength=2000, enzymeName=enzyme, tmpFolder="tmp",dictToStoreIDs="dict",
                        inMemory=True, compression=None)  # remove inMemory if you don't have enough RAM
        TR.parseInputData(dictLike=onename)
        TR.save(ensure(onename + "_parsed.frag"))

    list(map(parseMappedFile, in_files))
    "Merging files alltogether, applying filters"


    TR = HiCdataset("bla" + str(np.random.randint(100000000000)),genome=genomeFolder(workingGenome),
                    enzymeName=enzyme, tmpFolder="tmp", dictToStoreIDs="dict", mode="w", inMemory=True,
                    compression=None)
    TR.merge([i + "_parsed.frag" for i in in_files])
    TR.save(ensure(out_file + "_merged.frag"))

    for delFile in [i + "_parsed.frag" for i in in_files]:
        if os.path.exists(delFile):
            os.remove(delFile)

    TR = HiCdataset("bla" + str(np.random.randint(100000000000)), enzymeName=enzyme, maximumMoleculeLength=2000,
                    genome=genomeFolder(workingGenome), tmpFolder="tmp", dictToStoreIDs="dict",
                    mode='w', inMemory=True, compression=None)

    TR.load(out_file + "_merged.frag")

    #----------------------------Set of filters applied -------------

    # This was used for the publication
    #TR.filterLessThanDistance()
    # This is what we use now. It offers better filtering.
    TR.filterLessThanDistanceUpgraded()
    if TR.N == 0:
        return

    IDs = TR.rfragAbsIdxs1 * len(TR.rFragIDs) + TR.rfragAbsIdxs2
    uind = uniqueIndex(IDs)
    print("unique read removed: ", (-uind).sum(), "out of", len(uind))
    TR.metadata["350_removedNonUniqueInteractions"] = (-uind).sum()
    TR.maskFilter(uind)

    NN = TR.N
    fs = TR.fragmentSum()
    TR.fragmentFilter(fs < 9)
    TR.metadata["360_removedMoreThan8readsPerFragment"] = NN - TR.N
    TR.writeFilteringStats()

    TR.save(out_file + "_refined.frag")
    #------------------------End set of filters applied----------
    statFolder = os.path.join("statistics", out_file)
    TR.printMetadata(saveTo=statFolder + ".stat")
    for res in coolResolutions:
        TR.saveCooler(out_file + ".{0}.cool".format(res), res)


# This code is actually parsing datasets.tsv file
try:
    filename = sys.argv[1]
except:
    filename = "datasets.tsv"

fsplit = os.path.split(filename)
if len(fsplit[0]) > 0:
    os.chdir(fsplit[0])
filename = fsplit[1]

lines = open(filename).readlines()
lines = [i for i in lines if i[0] != "#"]
lines = [i.split() for i in lines if (len(i) > 3) and (i[0] != "#")]
for line in lines:
    if len(line) != 5:
        print("bad line", line)
        raise

for i in lines:
    if len(i) == 4:
        print(i)
        # Now we assume that enzyme is fifth argument in datasets.tsv or second commandline argument
        try:
            i.append(sys.argv[2])
        except:
            print("bad line", i)
            print("Please specify enzyme as a second command line argument, or as a fifth column")
            raise ValueError("Please specify enzyme as a second command line argument, or as a fifth column")

assert False not in [len(i) == 5 for i in lines]

dataFiles = lines
# experiment is defined by experiment name, replica name, genome and enzyme
experimentNames = set((i[1], i[2], i[3], i[4]) for i in dataFiles)
byExperiment = []
combinedExperimentNames = []

for experiment in experimentNames:
    # experiment is HeLa,R1,hg19,HindIII
    workingGenome = experiment[2]
    enzyme = 500
    filenames = [i[0] for i in dataFiles if (i[1], i[2], i[3], i[4]) == experiment]
    outName = "{0}-{1}-{3}".format(*experiment)

    # Filenames, path to save, genome, enzyme
    byExperiment.append((filenames, os.path.join(workingGenome, outName), workingGenome, enzyme))
    # merged files belonging to one expriment
    combinedExperimentNames.append((experiment[0], os.path.join(workingGenome, outName), workingGenome, enzyme))


# Now running refineDataset for each experiment
mypool = Pool(16)

list(mypool.map(refineDataset, byExperiment))

