from hiclib.fragmentHiC import HiCdataset 
from mirnylib.h5dict import h5dict  
import os,sys 
from mirnylib.plotting import mat_img
import numpy

if os.path.exists("test-1M.hm"): os.remove("test-1M.hm")
workingGenome = "hg18"
genomeFolder = "../../../data/hg18"
if not os.path.exists(genomeFolder):
    
    genomeFolder = sys.argv[1]        
    if not os.path.exists(genomeFolder):
        raise StandardError("Please provide hg18 Genome folder in the code or as a first argument") 

def source(ID):
    return os.path.join("%s-%s.hdf5" % (ID, workingGenome))   #determines path to the parsed file by ID



 
def refine_paper(filename,create = True):
    """filename[0] is a list of filenames of incoming files 
    filename[1] is a folder for outgoing file"""
    if create == True:        
        for onename in filename[0]:
            #Parsing individual files
            if not os.path.exists(onename): raise StandardError("path not found: %s" % onename)
            TR = HiCdataset(onename+"_parsed.frag",genome = genomeFolder,maximumMoleculeLength=500,override = True)             
            TR.parseInputData(dictLike = onename, enzymeToFillRsites = "HindIII")
            assert len(TR.DS) == 523790
            assert len(TR.ufragments) == 509379
              
        
        #Merging files alltogether, applying filters
        TR = HiCdataset(filename[1]+"_merged.frag",genome = genomeFolder,override = True)
        TR.merge([i+"_parsed.frag" for i in filename[0]])
        TR.flush()
                 
        TR = HiCdataset(filename[1]+"_refined.frag",genome = genomeFolder,override = True,autoFlush = False) 
        #because we do many operations, we disable autoFlush here 
        TR.load(filename[1]+"_merged.frag")
        TR.filterRsiteStart(offset = 5)
        assert len(TR.DS) == 509248                
        TR.filterDuplicates()
        assert len(TR.DS) == 508892
        #TR.save(filename[1]+".dat")
        TR.filterLarge()
        assert len(TR.DS) == 506081    
        TR.filterExtreme(cutH = 0.005, cutL = 0)
        assert len(TR.DS) == 490313
        TR.flush() 
    else: 
        #If merging & filters has already been done, just load files
         
        TR = HiCdataset(filename[1]+"_working.frag",override = True,genome = genomeFolder)
        TR.load(filename[1] +"_refined.frag")
        TR.rebuildFragments()

    print "----->Building Raw heatmap at two resolutions"
    TR.printStats()
    TR.saveHeatmap(filename[1] + "-1M.hm",1000000)
    from mirnylib.h5dict import h5dict
    a = h5dict(filename[1] + "-1M.hm")
    assert  a["heatmap"][::10,::10].sum()  == 12726
    
    print "---->Testing updateGenome method"
    from mirnylib.genome import Genome
    removeChromIDs = numpy.array([0,1,1,1,1]+[0]*17 + [1]+[0])    
    #print ((removeChromIDs[TR.chrms1] == 1) + (removeChromIDs[TR.chrms2] == 1) ).sum()    
    t =  ((removeChromIDs[TR.chrms1] == 1) * (removeChromIDs[TR.chrms2] == 1) ).sum() + ((removeChromIDs[TR.chrms1] == 1) * (TR.chrms2 == -1)).sum()
    newGenome = Genome(genomePath = genomeFolder, readChrms = ["2","3","4","5","X"])    
    TR.updateGenome(newGenome,removeSSreads = "trans")
    assert  len(TR.DS) == t
    
    a = h5dict(filename[1] + "-1M.hm")["heatmap"]
    

    print "----->Building RB heatmap"
    TR = HiCdataset(filename[1] + "_breaks.frag",genome = genomeFolder, override = True)
    TR.load(filename[1] + "_refined.frag")    
    TR.maskFilter((TR.dists1 > TR.maximumMoleculeLength) + (TR.dists2 > TR.maximumMoleculeLength) * TR.DS)
    
    print len(TR.DS)
    assert len(TR.DS) == 16848


map(refine_paper,
        [
      [(source("test"),
), "test","HindIII"]])

os.remove("test_breaks.frag")
os.remove("test-hg18.hdf5_parsed.frag")
os.remove("test_merged.frag")
os.remove("test_refined.frag")
print "Test finished successfully!"

