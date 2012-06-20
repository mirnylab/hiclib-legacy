import os,sys
from mirnylib.systemutils import setExceptionHook
sys.path.append(os.path.split(os.getcwd())[0])
from hiclib.binnedData import binnedData, binnedDataAnalysis,\
    experimentalBinnedData

import mirnylib.systemutils
mirnylib.systemutils
from hiclib.fragmentHiC import HiCdataset
from mirnylib.numutils import EIG, coarsegrain, project
import numpy
import mirnylib.plotting               
import scipy.stats
cr = scipy.stats.spearmanr
import cPickle
from mirnylib.plotting import mat_img,removeAxes,removeBorder, niceShow

import matplotlib.pyplot as plt 
import matplotlib

genomeVersion = "hg18"
myGenome = "../../data/%s"% genomeVersion

GM1M = "../../ErezPaperData/%s/GM-HindIII-%s-1M.hm"  % (genomeVersion, genomeVersion)
GM1MNcoI = "../../ErezPaperData/%s/GM-NcoI-%s-1M.hm"  % (genomeVersion, genomeVersion)
GM200k = "../../ErezPaperData/%s/GM-HindIII-%s-200k.hm"  % (genomeVersion, genomeVersion)
GM200kBreaks = "../../ErezPaperData/%s/GM-HindIII-%s-200k-breaks.hm"  % (genomeVersion, genomeVersion)
GM1MBreaks = "../../ErezPaperData/%s/GM-HindIII-%s-1M-breaks.hm"  % (genomeVersion, genomeVersion)
GM200kNcoI = "../../ErezPaperData/%s/GM-NcoI-%s-200k.hm"  % (genomeVersion, genomeVersion)
GMFrag = "../../ErezPaperData/%s/GM-NcoI-%s_refined.frag"  % (genomeVersion, genomeVersion)

tcc1M = "../../tcc/%s/tcc-HindIII-%s-1M.hm" % (genomeVersion, genomeVersion)
tcc200k = "../../tcc/%s/tcc-HindIII-%s-200k.hm" % (genomeVersion, genomeVersion)
workingFile1 = "../../ErezPaperData/working1"
workingFile2 = "../../ErezPaperData/working2"
workingFile3 = "../../ErezPaperData/working3"

def correctedScalingPlot():
    "Paper figure to compare scaling before/after correction"    
    plt.figure(figsize = (4,4))
    Tanay = binnedDataAnalysis(200000,genome = myGenome)
    Tanay.simpleLoad(GM200kBreaks,"GM-all")    
    Tanay.removePoorRegions()
    Tanay.removeDiagonal()
    Tanay.plotScaling("GM-all", label = "Raw data",color = "#A7A241")
    Tanay.iterativeCorrectWithSS()
    Tanay.plotScaling("GM-all",label = "Corrected",color = "#344370")
    ax = plt.gca() 
    mirnylib.plotting.removeAxes()
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

#correctedScalingPlot()


def doArmPlot():
    "Perform an average trans map for one dataset - paper figure"    
    
    Tanay = binnedDataAnalysis(resolution = 1000000, genome = myGenome)    
    Tanay.simpleLoad(GM1M,"GM-all")
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()
    Tanay.fakeCis()     
    Tanay.iterativeCorrectWithoutSS()
    plt.figure(figsize = (1.6,1.6))
    Tanay.averageTransMap("GM-all")
    mirnylib.plotting.removeBorder()
    cb = plt.colorbar(orientation = "vertical")
    cb.set_ticks([-0.05,0.05,0.15])
    for xlabel_i in cb.ax.get_xticklabels(): xlabel_i.set_fontsize(6)    
    mirnylib.plotting.niceShow()    
    


def averagePC2TransMap():
    """plot average trans map reconstructed from PC2 only 
    paper supplemental figure"""    
    Tanay = binnedDataAnalysis(1000000,genome = myGenome)
    #Tanay.loadNcoI()
    Tanay.simpleLoad(GM1M,"GM-all")    
    Tanay.removePoorRegions() 
    Tanay.removeZeros()    
    Tanay.truncTrans()
    Tanay.fakeCis()    
    Tanay.doEig()
    
    Tanay.restoreZeros(value = 0 )    
    PC2 = Tanay.EigDict["GM-all"][1]
    proj = PC2[:,None] * PC2[None,:]
    
    mm = PC2 != 0 
    mask = mm[:,None] * mm[None,:]
    proj[mask] += (2 * abs(proj.min()) )      
    Tanay.dataDict["BLA"] = proj
    Tanay.truncTrans(0.001)     
    Tanay.averageTransMap("BLA",chromMax = 22)
    plt.show()





def compareInterarmMaps():
    "plots witn 8 inetrarm maps - paper supplement figure"
    Tanay = binnedDataAnalysis(1000000,myGenome)

    Tanay.simpleLoad(GM1M,"GM-all")
    Tanay.simpleLoad(GM1MNcoI,"GM-NcoI")    
    Tanay.removeDiagonal()
    Tanay.removePoorRegions(cutoff  = 2)
    #Tanay.removeStandalone(3)
    fs = 10
    vmin = None
    vmax = None

    plt.subplot(421)
    plt.title("GM, HindIII, raw", fontsize = fs)
    Tanay.averageTransMap("GM-all",vmin = vmin, vmax = vmax)
    plt.colorbar()
    plt.subplot(422)
    plt.title("GM, NcoI, raw", fontsize = fs)
    Tanay.averageTransMap("GM-NcoI",vmin = vmin, vmax = vmax)
    plt.colorbar()
    
    Tanay.iterativeCorrectWithSS()
    vmin = None
    vmax = None
        
    plt.subplot(425)
    
    plt.title("GM, HindIII, with SS reads", fontsize = fs)
    Tanay.averageTransMap("GM-all",vmin = vmin, vmax = vmax)
    plt.colorbar()
    plt.subplot(426)
    
    plt.title("GM, NcoI, with SS reads", fontsize = fs)
    Tanay.averageTransMap("GM-NcoI",vmin = vmin, vmax = vmax)
    plt.colorbar()
    
    Tanay.iterativeCorrectWithoutSS()
    vmin = None
    vmax = None
    plt.subplot(423)    
    plt.title("GM, HindIII, no SS reads", fontsize = fs)
    Tanay.averageTransMap("GM-all",vmin = vmin, vmax=vmax)
    plt.colorbar()
    
    plt.subplot(424)    
    plt.title("GM, NcoI, no ss reads", fontsize = fs)
    Tanay.averageTransMap("GM-NcoI",vmin = vmin, vmax=vmax)
    plt.colorbar()
    Tanay.fakeCis()
    
    vmin = None
    vmax = None 
    plt.subplot(427)
    
    plt.title("GM, HindIII, trans only", fontsize = fs)
    Tanay.averageTransMap("GM-all",vmin = vmin, vmax=vmax)
    plt.colorbar()
    plt.subplot(428)
    
    plt.title("GM, NcoI, trans only", fontsize = fs)
    Tanay.averageTransMap("GM-NcoI",vmin = vmin, vmax=vmax)
    plt.colorbar()

    plt.show()    

def plotDiagonalCorrelation():
    "Correlation of diagonal bins - paper figure"
    S = 50
    x = numpy.arange(2,S)
    Tanay = binnedData(200000,myGenome)
    Tanay.simpleLoad(GM200k,"GM-HindIII")
    Tanay.simpleLoad(GM200kNcoI,"GM-NcoI")
    Tanay.simpleLoad(tcc200k,"TCC")
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()    
    Tanay.removeZeros()
    pairs = [("GM-HindIII","GM-NcoI"),("GM-HindIII","TCC"),("GM-NcoI","TCC")]
    cors = [[] for _ in pairs]
    for i in x:
        for j,pair in enumerate(pairs): 
            cors[j].append(cr(
                       numpy.diagonal(Tanay.dataDict[pair[0]],i),
                       numpy.diagonal(Tanay.dataDict[pair[1]],i)
                       )[0])
    
    Tanay.iterativeCorrectWithoutSS(M = 1)
    cors2 = [[] for _ in pairs]
    for i in x:
        for j,pair in enumerate(pairs):
            cors2[j].append(cr(
                       numpy.diagonal(Tanay.dataDict[pair[0]],i),
                       numpy.diagonal(Tanay.dataDict[pair[1]],i)
                       )[0])    
    Tanay.iterativeCorrectWithoutSS(M = 20) 
    cors3 = [[] for _ in pairs]
    for i in x:
        for j,pair in enumerate(pairs):
            cors3[j].append(cr(
                       numpy.diagonal(Tanay.dataDict[pair[0]],i),
                       numpy.diagonal(Tanay.dataDict[pair[1]],i)
                       )[0])


    matplotlib.rcParams['font.sans-serif']='Arial'
    

    #plt.figure(figsize = (2.3,1.8))
    print cors
    print cors2
    print cors3
    plt.figure(figsize = (10,3))
    ax = plt.gca()
    for j,pair in enumerate(pairs):
        plt.subplot(1,len(pairs),j)        
        fs = 8
        for xlabel_i in ax.get_xticklabels(): xlabel_i.set_fontsize(fs)
        for xlabel_i in ax.get_yticklabels(): xlabel_i.set_fontsize(fs)
        plt.title("%s vs %s" % pair)
        plt.plot(x/5.,cors3[j],color = "#E5A826",label = "Iterative")
        plt.plot(x/5.,cors2[j],color = "#28459A",label = "Single")
        plt.plot(x/5.,cors[j],color = "#E55726",label = "Raw")
        plt.xlabel("Genomic Separation, MB",fontsize = 8)
        plt.ylabel("Spearman correlation",fontsize = 8)
        plt.legend()
                        
        legend = plt.legend(prop={"size":6},loc = 9,handlelength=2)
        legend.draw_frame(False)
        plt.ylim((0,1))
        removeAxes(shift = 0)
    
    plt.show()     
     
#plotDiagonalCorrelation()

 

def plotCrossValidation():
    "main figure subplot with corss-validation"
    matplotlib.rcParams['font.sans-serif']='Arial'
    plt.figure(figsize = (1,1))
    FG = HiCdataset(workingFile1,myGenome)
    FG.load(GMFrag)     
    
    Tanay = binnedData(1000000)    
    Tanay.simpleLoad("GM-all-10p","GM-1")   #need to create these datasets using fragment-level analysis
    Tanay.simpleLoad("GM-all-90p","GM-9")
    Tanay.removePoorRegions()
    Tanay.iterativeCorrectWithSS()
    Tanay.removeZeros()
    b1,b2 = (Tanay.biasDict["GM-1"],Tanay.biasDict["GM-9"])
    cPickle.dump((b1,b2),open("CrossValidatioN",'wb'))
    ax = plt.gca() 
    b1,b2 = cPickle.load(open("CrossValidatioN",'rb'))
    print cr(b1,b2)
    plt.scatter(b1,b2,s=.7,color = "k",linewidth = 0)
    plt.xlabel(r"10% reads",fontsize = 8)
    plt.ylabel(r"90% reads",fontsize = 8)
    plt.xlim((0,1.5))
    plt.ylim((0,1.5))
    plt.xticks([0,0.5,1,1.5])
    plt.yticks([0,0.5,1,1.5])
    removeAxes(shift = 0 )
    fs = 6
    for xlabel_i in ax.get_xticklabels(): xlabel_i.set_fontsize(fs)
    for xlabel_i in ax.get_yticklabels(): xlabel_i.set_fontsize(fs)    
    plt.show() 

    
#plotCrossValidation()

def saddlePlot():
    "plot of values ordered by Eig1GW"
    
    #plt.figure(figsize = (1.5,1.5))
    plt.figure(figsize = (3,3))
    Tanay = binnedData(1000000)    
    Tanay.simpleLoad("../data/GM-all-hg18-1M","GM-all")
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()    
    Tanay.removeZeros()
    Tanay.fakeCis()        
    Tanay.iterativeCorrectWithoutSS()
    Tanay.doEig()
    PC = Tanay.EIG["GM-all"][:,0]
    if PC[0]>0: PC = -PC
    
    def reorder(data,array = PC):
        inds = numpy.argsort(array)
        ndata = data[inds,:]
        return ndata[:,inds]
    toplot = (coarsegrain(reorder(Tanay.dataDict["GM-all"]),60))
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
    removeBorder()
    mirnylib.plotting.niceShow()
      
#saddlePlot() 


def doCartoonPlot():
    "simple eigenvector decomposition cartoon"
    "educational purposes only"
    
    x  = 1.*numpy.arange(30) 
    a = numpy.sin(0.5 * x)
    b = 0.07* (15 - x)
    
    #mat_img(a[:,None] * a[None,])
    #mat_img(b[:,None] * b[None,:])
    noise = 0.2*numpy.random.randn(len(x),len(x))
    noise2 = 0.5* ( noise + noise.T)
    signal = a[:,None] * a[None,] + b[:,None] * b[None,:] + noise2 
    PC = EIG(signal)
    s1 = project(signal,PC[:,0])
    s2 = project(signal,PC[:,1])
    
    datas = [signal,s2,s1,s1+s2]
    datas = [i - 0.5*(i.min() + i.max()) for i in datas]
    mins = min([i.min() for i in datas])
    maxs = max([i.max() for i in datas])
    mylen = 5 
    begs = [0, 0, mylen,mylen ] 
    begs2 = [1,mylen+2, 1, mylen + 2]
    
    for i in xrange(4):
        ax = plt.subplot2grid((2*mylen+2,2*mylen), (begs2[i], begs[i]), rowspan=mylen, colspan = mylen )
        plt.imshow(datas[i],vmin = mins,vmax = maxs,interpolation = "nearest")
        removeBorder(ax = ax)
        
    ax = plt.subplot2grid((2*mylen+2,2*mylen),(0,mylen),colspan =mylen)
    removeBorder(ax = ax)
    plt.plot(PC[:,0],'k') 
    ax = plt.subplot2grid((2*mylen+2,2*mylen),(mylen+1,0),colspan =mylen)
    removeBorder(ax = ax)
    plt.plot(PC[:,1],'k') 
    
    plt.show()
    plt.imshow(signal - s1 - s2)
    mat_img(s1)
    mat_img(s2)
    mat_img(signal - s1 - s2)  


def compareWithGenomicFeatures():
    mirnylib.systemutils.setExceptionHook()
    Tanay = experimentalBinnedData(1000000,myGenome)
    Tanay.simpleLoad(GM1M,"GM-all")
    Tanay.simpleLoad(tcc1M,"tcc")
    
    datasets = {"CTCF": "wgEncodeBroadChipSeqSignalGm12878Ctcf.wig" ,"H3K27me3":    "wgEncodeBroadChipSeqSignalGm12878H3k27me3.wig",
                "H3K4me1" : "wgEncodeBroadChipSeqSignalGm12878H3k4me1.wig", "H3K4me3" : "wgEncodeBroadChipSeqSignalGm12878H3k4me3.wig",
                "H4K20me1" :   "wgEncodeBroadChipSeqSignalGm12878H4k20me1.wig", "H3K27ac" : "wgEncodeBroadChipSeqSignalGm12878H3k27ac.wig",
                "H3K36me3" : "wgEncodeBroadChipSeqSignalGm12878H3k36me3.wig", "H3K4me2" : "wgEncodeBroadChipSeqSignalGm12878H3k4me2.wig",
                "H3K9ac" :   "wgEncodeBroadChipSeqSignalGm12878H3k9ac.wig"
                }
    controlFile = "../../histoneMarks/hg18/wgEncodeBroadChipSeqSignalGm12878Control.wig"

    for key in datasets.keys(): 
        datasets[key] = os.path.join("../../histoneMarks/hg18",datasets[key])
        datasets[key] = (datasets[key],controlFile)
    datasets["DNAse"] =   ("../../histoneMarks/hg18/wgEncodeUwDnaseSeqRawSignalRep1Gm06990.bigWig",None)
    
     
        
    keys = datasets.keys()  
    

    Tanay.loadErezEigenvector1MB("../../ErezPaperData/eigenvectors")
    for key in keys: 
        Tanay.loadWigFile(datasets[key][0],key,control = datasets[key][1])
    Tanay.removeChromosome(22)
    Tanay.removeDiagonal()
    Tanay.removePoorRegions()
    Tanay.truncTrans()
    Tanay.fakeCis()    
    Tanay.removeZeros() 
    Tanay.doEig()
        
    E1 = Tanay.EigDict["GM-all"][0]
    E2 = Tanay.EigDict["tcc"][0]
    
    if scipy.stats.spearmanr(Tanay.trackDict["GC"],E1)[0] < 0: E1 *= -1
    if scipy.stats.spearmanr(Tanay.trackDict["GC"],E2)[0] < 0: E2 *= -1
    
    print "Eig-GC", scipy.stats.spearmanr(Tanay.trackDict["GC"],E1)
    print "Eigtcc-GC", scipy.stats.spearmanr(Tanay.trackDict["GC"],E2)
    
    Tanay.trackDict["Eig1GW_Lieberman"] = E1
    Tanay.trackDict["Eig1GW_TCC"] = E2
    Tanay.trackDict["EigLieberman2009"] = Tanay.trackDict["Erez"]
    keys = keys + ["GC","EigLieberman2009","Eig1GW_Lieberman","Eig1GW_TCC"]
    print "Key\t","key-GC\t","key-EigLieberman2009\t","key-EigGW_Lieberman\t","key-EigGW_TCC\t"
    for key in keys: 
        c1 = cr(Tanay.trackDict["GC"],Tanay.trackDict[key])
        c2 = cr(E1,Tanay.trackDict[key])
        c3 = cr(E2,Tanay.trackDict[key])
        c4 = cr(Tanay.trackDict["Erez"],Tanay.trackDict[key])
        print key, "\t%lf\t" % c1[0] ,"%lf\t" % c4[0] ,"%lf\t" % c2[0] ,"%lf" % c3[0]
    print
    #for key in keys:
    #    print key,"partial correlation is:",mirnylib.numutils.partialCorrelation(
    #                                                Tanay.trackDict[key], E1, Tanay.trackDict["GC"])
    raise 

#compareWithGenomicFeatures()

def plotTanayGenomicFeature():
    """Shows how genomic feature is spawned by Eig1, not Tanay domains
    paper supplementary  figure"""
    Tanay = experimentalBinnedData(1000000,myGenome)
    Tanay.simpleLoad(GM1M,"GM-all")        
    Tanay.loadTanayDomains()
    
    Tanay.loadWigFile("../../histoneMarks/hg18/wgEncodeUwDnaseSeqRawSignalRep1Gm06990.bigWig", label = "feature") 
    #Tanay.loadWigFile("../../histoneMarks/hg18/wgEncodeBroadChipSeqSignalGm12878H3k9ac.wig", label = "feature",
    #control = "../../histoneMarks/hg18/wgEncodeBroadChipSeqSignalGm12878Control.wig")
    #Tanay.loadWigFile("../../histoneMarks/hg18/wgEncodeBroadChipSeqSignalGm12878H3k4me3.wig", label = "feature",
    #control = "../../histoneMarks/hg18/wgEncodeBroadChipSeqSignalGm12878Control.wig")
        
    Tanay.removeDiagonal()
    Tanay.removePoorRegions()
    Tanay.removeZeros()
    Tanay.fakeCis()
    
    Tanay.doEig()
    E1 = Tanay.EigDict["GM-all"][0] 
    E2 = Tanay.EigDict["GM-all"][1]
    GC = Tanay.trackDict["GC"]    
    
    if scipy.stats.spearmanr(E1,GC)[0] < 0: E1 = -E1
    if scipy.stats.spearmanr(E2,GC)[0] < 0: E2 = -E2         
    
    TD = Tanay.trackDict["TanayDomains"]
    print scipy.stats.spearmanr(Tanay.trackDict["feature"],E1)
    
    plt.scatter(Tanay.trackDict["feature"],E1,c = TD,s=4,linewidth = 0)
    cm = plt.cm.get_cmap("jet")
        
    print "Our 2r is",(numpy.corrcoef(Tanay.trackDict["feature"],E1)[0,1])**2
    tset = set(TD)
    tfeature = numpy.zeros_like(TD,dtype = float)
    feature = Tanay.trackDict["feature"]
    for i in tset:
        #print i  
        #print numpy.mean(feature[TD==i])
        tfeature[TD == i] = numpy.mean(feature[TD==i])
        #print tfeature
    print "Tanay r2 is", (numpy.corrcoef(tfeature,E1)[0,1])**2  
    
    plt.legend([matplotlib.lines.Line2D([0],[0],color = cm(i),marker = "o",markersize = 8,linewidth = 0) for i in [0.333,0.666,0.999]],
               ["Active","Centromere proximal","Centromere distal"],loc = 2)
    
    plt.ylabel("Eig1GW")
    #plt.xlabel("UWashington DNAse")
    #plt.xlabel("H3K9ac ChIP-seq")
    plt.xlabel("H3K4me3 ChIP-seq")
    plt.title("Color represents domain from (Yaffe 2011)")
    niceShow(subplotAdjust = (0.13,0.11,0.97,0.92))
    
        