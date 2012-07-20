"""
This is an example file that shows different domain search algorithms 
"""

import mirnylib.systemutils 
mirnylib.systemutils.setExceptionHook()
import mirnylib.plotting
plotting = mirnylib.plotting 
import mirnylib.numutils
numutils = mirnylib.numutils  
import numpy 
na = numpy.array 

import matplotlib.pyplot as plt 

import hiclib 
from hiclib.binnedData import binnedData, binnedDataAnalysis
from hiclib.fragmentHiC import HiCStatistics 
from hiclib.domainSearch import gradientDescentDomains
from mirnylib.plotting import mat_img

import scipy.stats  
#os.chdir("/old/home/magus/HiC2011")

hmfile = "hg19/tcc-HindIII-hg19-1M.hm"
hmfileBreaks = "hg19/tcc-HindIII-hg19-1M-breaks.hm"

hmfile200 = "hg19/tcc-HindIII-hg19-200k.hm"
hmfileBreaks200 = "hg19/tcc-HindIII-hg19-200k-breaks.hm"

fragfile = "hg19/tcc-HindIII-hg19_refined.frag"

genome = "../data/hg19"

def showHeatmap():
    HM = binnedData(1000000,genome)
    HM.simpleLoad(hmfile,"tcc")
    HM.fakeCisOnce()
    mat_img(numpy.log(HM.dataDict["tcc"] + 1))
    
def correctedScalingPlot():
    
    plt.figure(figsize = (4,4))
    Tanay = binnedDataAnalysis(200000,genome)
    Tanay.simpleLoad(hmfileBreaks200,"GM-all")    
    Tanay.removePoorRegions()
    Tanay.plotScaling("GM-all", label = "Raw data",color = "#A7A241")
    Tanay.iterativeCorrectWithSS(M = 10)
    Tanay.plotScaling("GM-all",label = "Corrected",color = "#344370")
    ax = plt.gca() 
    mirnylib.plotting.removeAxes(shift = 0 )
    fs = 6
    plt.xlabel("Genomic distance (MB)",fontsize = 6 )
    plt.ylabel("Contact probability",fontsize = 6 )
    for xlabel_i in ax.get_xticklabels(): xlabel_i.set_fontsize(fs)
    for xlabel_i in ax.get_yticklabels(): xlabel_i.set_fontsize(fs)
    legend = plt.legend(loc=0,prop={"size":6})
    legend.draw_frame(False)
    plt.xscale("log")
    plt.yscale("log")
    plt.show() 

def averagePC2TransMap():
    "plot average trans map reconstructed from PC2 only"    
    Tanay = binnedDataAnalysis(1000000,genome)
    #Tanay.loadNcoI()
    Tanay.simpleLoad(hmfile,"GM-all")    
    Tanay.removePoorRegions() 
    Tanay.removeZeros()    
    Tanay.trunkTrans()
    Tanay.fakeCis()
    
    Tanay.doEig()
    Tanay.restoreZeros(value = 0 )    
    PC2 = Tanay.EigDict["GM-all"][1]
    proj = PC2[:,None] * PC2[None,:]
    
    mm = PC2 != 0 
    mask = mm[:,None] * mm[None,:]
    proj[mask] += (2 * abs(proj.min()) )      
    Tanay.dataDict["BLA"] = proj
    Tanay.trunkTrans(0.001)     
    Tanay.averageTransMap("BLA",chromMax = 22)
    plt.show()


def saddlePlot():
    "plot of values ordered by Eig1GW"
    
    #plt.figure(figsize = (1.5,1.5))
    plt.figure(figsize = (3,3))
    Tanay = binnedData(1000000,genome)    
    Tanay.simpleLoad(hmfile,"GM-all")
    Tanay.loadGC()
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()
    Tanay.removeZeros()
    Tanay.fakeCis()        
    Tanay.iterativeCorrectWithoutSS(M = 25)
    Tanay.doEig()
    PC = Tanay.EigDict["GM-all"][0]
    if PC[0]>0: PC = -PC    
    
    def reorder(data,array = PC):
        inds = numpy.argsort(array)
        ndata = data[inds,:]
        return ndata[:,inds]
    toplot = (mirnylib.numutils.coarsegrain(reorder(Tanay.dataDict["GM-all"]),60))
    toplot /= toplot.mean()
    toplot = numpy.log(toplot)
    sh = toplot.shape
    toplot = toplot.reshape((-1))
    ag = numpy.argmax(toplot)
    toplot[ag] = 0 
    toplot[ag] = numpy.max(toplot)
    toplot.shape = sh
    toplot[0,-1] = toplot[0,-2]
    toplot[-1,0] = toplot[-2,0]  
    
    plt.imshow(toplot,vmin = toplot.min(),vmax = toplot.max(),interpolation="nearest")
    cbar = plt.colorbar(orientation = "vertical")
    #labels = ["10","100","1000","10000"]        
    #cbar.ax.set_xticklabels(labels)    
    cbar.ax.set_xlabel("Log(relative contact probability)", fontsize = 6)
    for xlabel_i in cbar.ax.get_xticklabels(): xlabel_i.set_fontsize(6)    
    cbar.set_ticks([-0.5,0,0.5,1])
    plotting.removeBorder()
    plotting.niceShow()
    

def doArmPlot():
    
    Tanay = binnedDataAnalysis(1000000,genome)    
    Tanay.simpleLoad(hmfile,"GM-all")
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()
    Tanay.fakeCis()     
    Tanay.iterativeCorrectWithoutSS(M = 20 )
    plt.figure(figsize = (3.6,3.6))
    Tanay.averageTransMap()

    plotting.removeBorder()
    cb = plt.colorbar(orientation = "vertical")
    cb.set_ticks([-0.05,0.05,0.15])
    for xlabel_i in cb.ax.get_xticklabels(): xlabel_i.set_fontsize(6)    
    plotting.niceShow()    
         
def plotScalingsByDomains(build = True):
    """An example script that shows calculation of domains, and then plotting inter-domain and intra-domain scaling. 
    This script uses advanced features, including heatmap creation and low-level analysis of domains.
    It assumes that data has been fragment-corrected already 
    """    
    resolution = 350000 #resolution to find domains 
    
 
    TR = HiCStatistics("working/working.frag",genome, maximumMoleculeLength=500)
    TR.load(fragfile)
    TR.saveHeatmap("working/workingHeatmap", resolution = resolution)   #Simply load corrected dataset and save it as a heatmap
    TR.maskFilter((TR.chrms1 == TR.chrms2) * (TR.strands1 == TR.strands2))   #making it easier to calculate scaling
    
    HM = hiclib.binnedData.binnedDataAnalysis(resolution , genome)  #load heatmap 
    #Following standart protocol of obtaining domains with advanced statistical model 
    HM.simpleLoad("working/workingHeatmap", "GM")    
    HM.removePoorRegions(cutoff = 3)
    HM.removeStandalone(4)
    HM.loadGC()
    HM.trunkTrans(high = 0.0002)       
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
    chromosomeStarts = numpy.array(HM.genome.chrmStartsBinCont)    #chromosome starts
    
    fragments= TR.ufragments   #Now filtering fragments 
    pos = fragments % TR.fragIDmult
    chrom = fragments / TR.fragIDmult   
    mask1 = domainIndexRescaled[(pos/resolution) + chromosomeStarts[chrom]] == 1
    mask2 = domainIndexRescaled[(pos/resolution) + chromosomeStarts[chrom]] == -1    #fragment masks for open/closed regions 
    
     

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
    
    mirnylib.plotting.niceShow()

    
    
plotScalingsByDomains()

