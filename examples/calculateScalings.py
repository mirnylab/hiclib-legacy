from hiclib.fragmentHiC import HiCdataset,HiCStatistics
import hiclib.binnedData
import mirnylab.numutils 
import hiclib.domainSearch 
import numpy
import matplotlib.pyplot as plt   
import os
import mirnylab.plotting 
import scipy.stats

dataDir = "../.."   #base folder with project data



genomeFolder = "../../data/hg19"
dataFile = "../../ErezPaperData/hg19/GM-HindIII-hg19_refined.frag"
dataFile2 = "hg19/GM-NcoI-hg19_refined.frag"

if os.path.exists("working") == False: os.mkdir("working")

def plotScalingsByChromosome():
    """ An example script that plots scalings for different chromosomes with different color.
    This shows the use of "plotScaling" with user-defined regions.  
    It also shows how easy it is to get genome data from the Genome class. 
    This method is absolutely transparent and can be used with any genome/any data
    """
    TR = HiCStatistics("working/working.frag",genome = genomeFolder, maximumMoleculeLength=500,override = True )    
    TR.load(dataFile)
    """Filter cis reads from the same strand to reduce # reads 
    (it would be taken care of anyways, it's just faster this way)"""
    TR.maskFilter((TR.strands1 == TR.strands2) * (TR.chrms1 == TR.chrms2))
    
    TR.plotScaling( color = "black", linewidth = 2,excludeNeighbors=3,enzyme = "HindIII",mindist = 20000)
     
    for chromosome in xrange(TR.genome.chrmCount):
        regions = [(chromosome,0,TR.genome.cntrMids[chromosome]),
                   (chromosome,TR.genome.cntrMids[chromosome], TR.genome.chrmLens[chromosome])]
        TR.plotScaling(regions = regions,excludeNeighbors=3,enzyme = "HindIII",
                       color = (plt.cm.get_cmap("jet")(0.3 + (0.5 * chromosome ) / TR.genome.chrmCount)),mindist = 20000)
         
    
    plt.xlabel("Genomic distance (bp)")
    plt.ylabel("Normalized Pc")        
    plt.show()
        
plotScalingsByChromosome()
exit()        

def compareWeightsNoWeights():
    TR = HiCStatistics("working/working.frag",genome = genomeFolder, maximumMoleculeLength=500)
    TR.load(dataFile)
    TR.maskFilter((TR.strands1 == TR.strands2) * (TR.chrms1 == TR.chrms2))    
    TR.plotScaling( excludeNeighbors=3,enzyme = "HindIII", mindist = 100,label = "No weights, Hind")
    TR.plotScaling( useWeights= True, excludeNeighbors=3,enzyme = "HindIII", mindist = 100,label = "weights, Hind")

    TR = HiCStatistics("working/working.frag",genome = genomeFolder, maximumMoleculeLength=500)
    TR.load(dataFile2)
    TR.maskFilter((TR.strands1 == TR.strands2) * (TR.chrms1 == TR.chrms2))    
    TR.plotScaling( excludeNeighbors=3,enzyme = "NcoI", mindist = 100 ,label = "No weights, NcoI")
    TR.plotScaling( useWeights= True, excludeNeighbors=3,enzyme = "NcoI", mindist = 100 ,label = "weights, NcoI")
    mirnylab.plotting.niceShow()
    
#compareWeightsNoWeights()
           
def plotScalingsByDomains(build = True):
    """An example script that shows calculation of domains, and then plotting inter-domain and intra-domain scaling. 
    This script uses advanced features, including heatmap creation and low-level analysis of domains.
    It assumes that data has been fragment-corrected already 
    """    
    resolution = 2000000 #resolution to find domains 
    
 
    TR = HiCStatistics("working/working.frag",genomeFolder = genomeFolder, maximumMoleculeLength=500)
    TR.load(dataFile)
    TR.saveHeatmap("working/workingHeatmap", resolution = resolution)   #Simply load corrected dataset and save it as a heatmap
    
    HM = hiclib.binnedData.binnedDataAnalysis(resolution = resolution, genomeFolder = genomeFolder)  #load heatmap 
    #Following standart protocol of obtaining domains with advanced statistical model 
    HM.simpleLoad("working/workingHeatmap", "GM")    
    HM.removePoorRegions(cutoff = 5)
    HM.removeStandalone(5)
    HM.loadGC()
    HM.TrunkTrans(high = 0.0002)       
    HM.fakeCisOnce()   
    HM.removeZeros()
    data = HM.dataDict["GM"]
    #plotting.mat_img(data)              #check that we're feeding what we want to domain finder 
    domains = hiclib.domainSearch.exponentialDomains(data)   #most advanced domain finder
    
    cr = scipy.stats.spearmanr(domains,HM.trackDict["GC"])   #be sure that positive domains are always open
    print "correlation with GC is %lf" % cr[0]
    if cr[0] < 0: domains = -domains                         #inverting domains if spearmanr(GC,dom) < 0    
             
    domainIndex = numpy.zeros(len(domains),int)   #Building 2 groups of bins: first third and last third
    p5,p35,p65,p95 = numpy.percentile(domains,[5,35,65,95])
    domainIndex[(p5 < domains)* (domains < p35)] = -1 
    domainIndex[(p65 < domains) * (domains < p95)] = 1   
    HM.trackDict["domains"] = domainIndex
    HM.restoreZeros(value = 0)      #now we have track "domainIndex" that stores indexes of megabases with proper domains
    domainIndexRescaled = HM.trackDict["domains"]
    chromosomeStarts = numpy.array(HM.genome.chromosomeStarts)    #chromosome starts
    
    fragments= TR.ufragments   #Now filtering fragments 
    pos = fragments % TR.fragIDmult
    chr = fragments / TR.fragIDmult   
    mask1 = domainIndexRescaled[(pos/resolution) + chromosomeStarts[chr]] == 1
    mask2 = domainIndexRescaled[(pos/resolution) + chromosomeStarts[chr]] == -1    #fragment masks for open/closed regions 

    TR.plotScaling(fragids1 = mask2, fragids2 = mask2, useWeights = False,
                    excludeNeighbors = 2, enzyme = "HindIII",  withinArms = True, 
                    mindist = 20000 , appendReadCount= False, label = "Closed-Closed", linewidth = 2)
    
    
    TR.plotScaling(fragids1 = mask1, fragids2 = mask2, useWeights = False,
                    excludeNeighbors = 2, enzyme = "HindIII",  withinArms = True, 
                    mindist = 20000 , appendReadCount= False,  label = "Open-Closed", linewidth = 2)  

    TR.plotScaling(fragids1 = mask1, fragids2 = mask1, useWeights = False,
                    excludeNeighbors = 2, enzyme = "HindIII",  withinArms = True, 
                    mindist = 20000 , appendReadCount= False, label = "Open-Open", linewidth = 2)  
    
    plt.xlabel("genomic separation")
    plt.ylabel("Pc")    
    
    mirnylab.plotting.niceShow()
    raise 
    
plotScalingsByDomains()  
    
    
    
    
    
    
    
    
     
     
    
    
    
    
    
     
    
     
    


TR.calculateWeights()
TR.plotScaling(useWeights = False, excludeNeighbors = 3, label = "NcoI - no weights", enzyme = "NcoI",mindist = 500)
TR.plotScaling(useWeights = True, excludeNeighbors = 3, label = "NcoI - weights", enzyme = "NcoI",mindist = 500)

TR = HiCStatistics("working2",genomeFolder = genomeFolder(workingGenome), maximumMoleculeLength=500,override = True)
TR.load("hg19/GM-NcoI-hg19_refined")
TR.maskFilter((TR.strands1 == TR.strands2) * (TR.chrms1 == TR.chrms2))
TR.calculateWeights()
TR.maskFilter(TR.chrms1 == TR.chrms2)
TR.plotScaling(useWeights = False, excludeNeighbors = 3, label = "HindIII - no weights", enzyme = "HindIII",mindist = 500)
TR.plotScaling(useWeights = True, excludeNeighbors = 3, label = "HindIII - weights", enzyme = "HindIII",mindist = 500)
plotting.niceShow()



TR.plotScaling(useWeights = True, excludeNeighbors = 0, label = "None", enzyme = "NcoI",mindist = 100)
TR.plotScaling(useWeights = True, excludeNeighbors = 1, label = "1", enzyme = "NcoI",mindist = 100)

TR.plotScaling(useWeights = True, excludeNeighbors = 4, label = "4", enzyme = "NcoI",mindist = 100)
TR.plotScaling(useWeights = True, excludeNeighbors = 6, label = "6", enzyme = "NcoI",mindist = 100)
plt.show()

exit()
plt.show()
f1 = TR.fragids1
f2 = TR.fragids2

#f1mod, f2mod = TR.genome.getPairsLessThanDistance(f1,f2,cutoffDistance = 3,enzymeName = "HindIII")
raise  
dists = TR.genome.getFragmentDistance(f1, f2,"HindIII")
mask = dists < 1000
plt.scatter((f1-f2)[mask],dists[mask])
plt.show()
