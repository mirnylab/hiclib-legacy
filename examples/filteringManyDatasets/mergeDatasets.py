"""
This example script does the following: 

-Loads individual files from the location defined by source(ID) 
-Parses individual files in memory  (inMemory = True in "for onename in in_files))
   -If your system does not have enough memory, you might need to switch to hdf5 here.
-Merges files corresponding to the same experiment together, on the HDD. 
-Filters datasets, builds heatmaps 
-Combines multiple replicas of the same experiment together, builds heatmaps

--Locations of individual files are defined by source(ID) function.
--Datasets are defined in the datasets.tsv file
--genome is defined by genomeFolder function, and workingGenome constant.
--Resulting files are placed to the workingDirectory folder  

Warnings: 
    Running this over NFS might cause unexpected slow-downs because NFS is
    unhappy with repeated read/write cycles to the same file 
    
    You could do entire thing in memory, if you have RAM or your datasets are small. 
    Actually, using HDF5 is then equivalent to storing compressed data in RAM, 
    and might be in fact pretty fast. 
 
"""

from hiclib.fragmentHiC import HiCdataset   
import os 

def genomeFolder(name):
    return os.path.join("/home/magus/HiC2011/data",name)   #Fetch genome folder by genome name

workingGenome = "hg19"
workingDirectory = "%s" % workingGenome

def source(ID):
    return "/net/tamm/home/magus/gigareadHiC/data/%s.fastq.hdf5" % ID


if os.path.exists(workingDirectory) == False:
    print "working directory does not exist"
    exit()    

 
def refineDatasets(filenames,create = True, delete = True ,parseInMemory = True):
    """
    Parameters
    ----------
    
    filenames[0] is a list of filenames of incoming files 
    filenames[1] is a folder for outgoing file
    
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
        
    if create == True:  #if we need to parse the input files (.hdf5 from mapping).         
        for onename in in_files:
            #Parsing individual files
            if not os.path.exists(source(onename)): raise StandardError("path not found: %s" % onename)
            if parseInMemory == True: 
                TR = HiCdataset("bla",genome = genomeFolder(workingGenome),
                                maximumMoleculeLength=500,override = True,
                                inMemory = True)   #remove inMemory if you don't have enough RAM
                #create dataset in memory, parse and then save to destination 
                
                TR.parseInputData(dictLike = source(onename))
                TR.save(onename + "_parsed.frag")
            else: 
                #Create dataset at destination, parse on HDD, then no need to save. 
                TR = HiCdataset(onename + "_parsed.frag",
                                genome = genomeFolder(workingGenome),
                                maximumMoleculeLength=500,override = True)
                TR.parseInputData(dictLike = source(onename))
                
        
        "Merging files alltogether, applying filters"
        TR = HiCdataset(out_file+"_merged.frag",
                        genome = genomeFolder(workingGenome),
                        override = True)
        
        TR.merge([i+"_parsed.frag" for i in in_files])   #Merge in all parsed files from one experiment
        
        if delete == True:   #cleaning up parsed files
            for delFile in [i+"_parsed.frag" for i in in_files]:
                os.remove(delFile)                 
        TR.flush()
        
        "Now opening new dataset for refined data, and performing all the filtering "         
        TR = HiCdataset(out_file+"_refined.frag",
                        genome = genomeFolder(workingGenome),
                        override = True,autoFlush = False) #because we do many operations, we disable autoFlush here                  
        TR.load(out_file+"_merged.frag")
        TR.filterRsiteStart(offset = 5)                
        TR.filterDuplicates()
        #TR.save(out_file+".dat")
        TR.filterLarge()    
        TR.filterExtreme(cutH = 0.005, cutL = 0)
        TR.flush()        
         
    else: 
        #If merging & filters has already been done, just load files                 
        TR = HiCdataset(out_file+"_working.frag",override = True,genome = genomeFolder(workingGenome))
        TR.load(out_file +"_refined.frag")        

    print "----->Building Raw heatmap at two resolutions"
    TR.printStats()    
    TR.saveHeatmap(out_file + "-200k.hm",200000)
    TR.saveHeatmap(out_file + "-500k.hm",500000)
    TR.saveHeatmap(out_file + "-1M.hm",1000000)
    TR.saveHeatmap(out_file + "-2M.hm",2000000)
    
    
    print "----->Building RB heatmap"
    TR = HiCdataset(out_file + "_breaks.frag",genome = genomeFolder(workingGenome), override = True)
    TR.load(out_file + "_refined.frag")    
    TR.maskFilter((TR.dists1 > TR.maximumMoleculeLength) + (TR.dists2 > TR.maximumMoleculeLength) * TR.DS)
    TR.printStats()    
    TR.saveHeatmap(out_file + "-200k-breaks.hm",200000)
    TR.saveHeatmap(out_file + "-500k-breaks.hm",500000)
    TR.saveHeatmap(out_file + "-1M-breaks.hm",1000000)
    TR.saveHeatmap(out_file + "-2M-breaks.hm",2000000)
    

#This code is actually parsing datasets.tsv file
dataFiles = open("datasets.tsv").readlines() 
dataFiles = [i.split() for i in dataFiles if (len(i) > 3) and (i[0] != "#")]
assert False not in [len(i) == 3 for i in dataFiles]

experimentNames  = set((i[1],i[2]) for i in dataFiles)
byExperiment = []
newExperimentNames = []
for experiment in experimentNames:
    filenames = [i[0] for i in dataFiles if (i[1],i[2]) == experiment] 
    outName = str(experiment[0]) + "-" + str(experiment[1])
    byExperiment.append((filenames,os.path.join(workingDirectory,outName),"HindIII"))
    newExperimentNames.append((experiment[0],os.path.join(workingDirectory,outName)))
    

for i in byExperiment: refineDatasets(i, create = True)    
experiments = set(i[0] for i in newExperimentNames)

for experiment in experiments:
    myExperimentNames = [i[1] + "_refined.frag" for i in newExperimentNames if i[0] == experiment]
    assert len(myExperimentNames) > 0 
    if len(myExperimentNames) > 1:
        TR = HiCdataset(os.path.join(workingDirectory,"%s-all_refined.frag" % experiment),genome = genomeFolder(workingGenome))        
        TR.merge(myExperimentNames)
        TR.saveHeatmap(os.path.join(workingDirectory,"%s-all-100k.hm" % experiment),100000)
        TR.saveHeatmap(os.path.join(workingDirectory,"%s-all-200k.hm" % experiment),200000)
        TR.saveHeatmap(os.path.join(workingDirectory,"%s-all-500k.hm" % experiment),500000)
        TR.saveHeatmap(os.path.join(workingDirectory,"%s-all-1M.hm" % experiment),1000000)



#map(refine_paper,
#        [((source("SRR027961"),
#       source("SRR027960")),   os.path.join(workingDirectory, "GM-NcoI-%s" % workingGenome ),"NcoI"),  
#      ((source("SRR027956"),
#        source("SRR027957"),
#        source("SRR027958"),
#        source("SRR027959")),  os.path.join(workingDirectory, "GM-HindIII-%s" % workingGenome ),"HindIII")])

