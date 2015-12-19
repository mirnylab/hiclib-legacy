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

from hiclib.fragmentHiC import HiCdataset
from mirnylib.systemutils import fmap,setExceptionHook
from mirnylib.genome import Genome 
import numpy as np 
import os
import sys
from defineGenome import getGenome


setExceptionHook()
def ensure(f):
    d = os.path.dirname(f)
    if os.path.isdir(d):
        return f
    else:
        try:
            os.makedirs(d)
        except:
            raise ValueError("Cannot create directory")
    return f


wholeGenomeResolutionsKb = [1000,500,200,100,40]
byChromosomeResolutionsKb = [40,20,10] 
HiResWithOverlapResolutionsKb = [10,5,2,1]  #add 10 here if you have more than 16GB RAM 
SuperHiResWithOverlapResolutionsKb = [] 
skip = 1 #how many to skip for single replica datasets



def refineDataset(filenames, create=True, delete=True, parseInMemory=True):
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

    statFolder = os.path.join("statistics", out_file)

    workingGenome = filenames[2]
    enzyme = filenames[3]

    if create == True:  # if we need to parse the input files (.hdf5 from mapping).
        def parse_onename(onename):
            np.random.seed()
            #Parsing individual files
            if parseInMemory == True:
                finalname = onename + "_parsed.frag"
                #if not os.path.exists(finalname):
                if True:

                    #create dataset in memory, parse and then save to destination
                    TR = HiCdataset("bla" + str(np.random.randint(100000000000)), genome=getGenome(workingGenome),
                                    maximumMoleculeLength=500,enzymeName = enzyme,tmpFolder = "tmp",
                                    inMemory=True)  # remove inMemory if you don't have enough RAM

                    TR.parseInputData(dictLike=onename)
                    folder = os.path.split(onename)[0]
                    print(onename)
                    TR.save(ensure(finalname))
                    folder, fname = os.path.split(onename)
                    statSubFolder = os.path.join(statFolder, folder)
                   

                    TR.printMetadata(saveTo=ensure(os.path.join(statSubFolder, fname + ".stat")))
                else:
                    print("skipping parsed: ", onename)
            else:
                #Create dataset at destination, parse on HDD, then no need to save.
                TR = HiCdataset(ensure(onename + "_parsed.frag"),
                                genome=getGenome(workingGenome),enzymeName = enzyme,tmpFolder = "tmp",
                                maximumMoleculeLength=500, mode='w')
                TR.parseInputData(dictLike=onename, enzymeToFillRsites=enzyme)
                TR.printMetadata(saveTo=ensure(os.path.join(statFolder, onename + ".stat")))
        list(map(parse_onename, in_files))
        "Merging files alltogether, applying filters"
        TR = HiCdataset(ensure(out_file + "_merged.frag"),
                        genome=getGenome(workingGenome),enzymeName = enzyme,tmpFolder = "tmp",dictToStoreIDs="h5dict",
                        mode="w")
        TR.merge([i + "_parsed.frag" for i in in_files])
            #Merge in all parsed files from one experiment

        if delete == True:  # cleaning up parsed files
            for delFile in [i + "_parsed.frag" for i in in_files]:
                os.remove(delFile)

        "Now opening new dataset for refined data, and performing all the filtering "
        TR = HiCdataset(out_file + "_refined.frag",enzymeName = enzyme,
                        genome=getGenome(workingGenome),tmpFolder = "tmp",dictToStoreIDs="h5dict",
                        mode='w')
        TR.load(out_file + "_merged.frag")

        #----------------------------Set of filters applied -------------
        TR.filterDuplicates()
        #TR.save(out_file+".dat")
        TR.filterLarge(10000,10)
        TR.filterExtreme(cutH=0.001, cutL=0)
        TR.writeFilteringStats()
        TR.printMetadata(saveTo=statFolder + ".stat")

        #------------------------End set of filters applied----------

    else:
        #If merging & filters has already been done, just load files
        TR = HiCdataset(out_file + "_working.frag",enzymeName = enzyme,
                        mode='w', genome=getGenome(workingGenome))
        TR.load(out_file + "_refined.frag")
        TR.printMetadata(saveTo=statFolder + ".stat")

    print("----->Building Raw heatmap at different resolutions")
    TR.printStats()
    for res in wholeGenomeResolutionsKb:    
        TR.saveHeatmap(out_file + "-{0}k.hm".format(res), res*1000)
    for res in byChromosomeResolutionsKb: 
        TR.saveByChromosomeHeatmap(out_file + "-{0}k.byChr".format(res), res*1000)
    for res in HiResWithOverlapResolutionsKb[:-skip]:
        TR.saveHiResHeatmapWithOverlaps(out_file + "-{0}k_HighRes.byChr".format(res), res*1000)        
    for res in SuperHiResWithOverlapResolutionsKb[:-skip]:
        TR.saveSuperHighResMapWithOverlaps(out_file + "-{0}k_HighRes.byChr".format(res), res*1000)            


#This code is actually parsing datasets.tsv file

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
        #Now we assume that enzyme is fifth argument in datasets.tsv or second commandline argument
        try:
            i.append(sys.argv[2])
        except:
            print("bad line", i)
            print("Please specify enzyme as a second command line argument, or as a fifth column")
            raise ValueError("Please specify enzyme as a second command line argument, or as a fifth column")

assert False not in [len(i) == 5 for i in lines]

dataFiles = lines
#experiment is defined by experiment name, replica name, genome and enzyme 
experimentNames = set((i[1], i[2], i[3], i[4]) for i in dataFiles)
byExperiment = []
combinedExperimentNames = []

for experiment in experimentNames:
    #experiment is HeLa,R1,hg19,HindIII
    workingGenome = experiment[2]
    enzyme = experiment[3]
    filenames = [i[0] for i in dataFiles if (i[1], i[2], i[3], i[4]) == experiment]
    outName = "{0}-{1}-{3}".format(*experiment)

    #Filenames, path to save, genome, enzyme
    byExperiment.append((filenames, os.path.join(workingGenome, outName), workingGenome, enzyme))
    # merged files belonging to one expriment
    combinedExperimentNames.append((experiment[0], os.path.join(workingGenome, outName), workingGenome, enzyme))


#Now running refineDataset for each experiment
for i in byExperiment:
    refineDataset(i, create=True, delete=True)
    pass


#Now merging different experiments alltogether
#note that the first column is not here, as it is a replica 
experiments = set([(i[0], i[2], i[3]) for i in combinedExperimentNames])
print(experiments)

for experiment in experiments:
    workingGenome = experiment[1]
    myExperimentNames = [i[1] + "_refined.frag" for i in combinedExperimentNames if (i[0], i[2], i[3]) == (experiment[0], experiment[1],experiment[2])]    
    assert len(myExperimentNames) > 0
    if len(myExperimentNames) > 0:
        #If we have more than one experiment (replica) for the same data, we can combine. 
        TR = HiCdataset(os.path.join(workingGenome, "%s-all-%s_refined.frag" %
                                     (experiment[0],experiment[2])), genome=getGenome(workingGenome),
                                     enzymeName = experiment[2],tmpFolder = "tmp",dictToStoreIDs="h5dict")
        statSaveName = os.path.join("statistics", workingGenome, "%s-all-%s_refined.stat" % (experiment[0], experiment[2]))

        TR.merge(myExperimentNames)
        TR.printMetadata(saveTo=statSaveName)
        for res in wholeGenomeResolutionsKb:    
            TR.saveHeatmap(os.path.join(workingGenome, "%s-all-%s-{0}k.hm" % (experiment[0], experiment[2])).format(res), res*1000)
        for res in byChromosomeResolutionsKb: 
            TR.saveByChromosomeHeatmap(os.path.join(workingGenome, "%s-all-%s-{0}k.byChr" % (experiment[0], experiment[2])).format(res), res*1000)
        for res in HiResWithOverlapResolutionsKb:
            TR.saveHiResHeatmapWithOverlaps(os.path.join(workingGenome, "%s-all-%s-{0}k_HighRes.byChr" % (experiment[0], experiment[2])).format(res), res*1000)
        for res in SuperHiResWithOverlapResolutionsKb:
            TR.saveSuperHighResMapWithOverlaps(os.path.join(workingGenome, "%s-all-%s-{0}k_SuperHighRes.byChr" % (experiment[0], experiment[2])).format(res), res*1000)

