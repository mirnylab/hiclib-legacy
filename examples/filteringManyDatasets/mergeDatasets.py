from hiclib.fragmentHiC import HiCdataset   
import os 

 


def genomeFolder(name):
    return os.path.join("/home/magus/HiC2011/data",name)   #Fetch genome folder by genome name

workingGenome = "hg19"
workingDirectory = "%s" % workingGenome

def source(ID):
    return "data/%s.fastq.hdf5" % ID


if os.path.exists(workingDirectory) == False:
    print "working directory does not exist"
    exit()    




 
def refineDatasets(filename,create = True):
    """filename[0] is a list of filenames of incoming files 
    filename[1] is a folder for outgoing file"""
    if create == True:        
        for onename in filename[0]:
            #Parsing individual files
            if not os.path.exists(onename): raise StandardError("path not found: %s" % onename)
            TR = HiCdataset(onename+"_parsed.frag",genome = genomeFolder(workingGenome),maximumMoleculeLength=500,override = True)             
            TR.parseInputData(dictLike = onename)              
        
        #Merging files alltogether, applying filters
        TR = HiCdataset(filename[1]+"_merged.frag",genome = genomeFolder(workingGenome),override = True)
        TR.merge([i+"_parsed.frag" for i in filename[0]])
        TR.flush()
                 
        TR = HiCdataset(filename[1]+"_refined.frag",genome = genomeFolder(workingGenome),override = True,autoFlush = False) 
        #because we do many operations, we disable autoFlush here 
        TR.load(filename[1]+"_merged.frag")
        TR.filterRsiteStart(offset = 5)                
        TR.filterDuplicates()
        #TR.save(filename[1]+".dat")
        TR.filterLarge()    
        TR.filterExtreme(cutH = 0.005, cutL = 0)
        TR.flush() 
    else: 
        #If merging & filters has already been done, just load files
         
        TR = HiCdataset(filename[1]+"_working.frag",override = True,genome = genomeFolder(workingGenome))
        TR.load(filename[1] +"_refined.frag")
        TR.rebuildFragments()

    print "----->Building Raw heatmap at two resolutions"
    TR.printStats()
    TR.saveHeatmap(filename[1] + "-100k.hm",100000)
    TR.saveHeatmap(filename[1] + "-200k.hm",200000)
    TR.saveHeatmap(filename[1] + "-500k.hm",500000)
    TR.saveHeatmap(filename[1] + "-1M.hm",1000000)
    
    
    print "----->Building RB heatmap"
    TR = HiCdataset(filename[1] + "_breaks.frag",genome = genomeFolder(workingGenome), override = True)
    TR.load(filename[1] + "_refined.frag")    
    TR.maskFilter((TR.dists1 > TR.maximumMoleculeLength) + (TR.dists2 > TR.maximumMoleculeLength) * TR.DS)
    TR.printStats()
    TR.saveHeatmap(filename[1] + "-100k-breaks.hm",100000)
    TR.saveHeatmap(filename[1] + "-200k-breaks.hm",200000)
    TR.saveHeatmap(filename[1] + "-500k-breaks.hm",500000)
    TR.saveHeatmap(filename[1] + "-1M-breaks.hm",1000000)



dataFiles = open("datasets.tsv").readlines() 
dataFiles = [i.split() for i in dataFiles if len(i) > 3]
assert False not in [len(i) == 3 for i in dataFiles]

experimentNames  = set((i[1],i[2]) for i in dataFiles)
byExperiment = []
newExperimentNames = []
for experiment in experimentNames:
    filenames = [source(i[0]) for i in dataFiles if (i[1],i[2]) == experiment] 
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

