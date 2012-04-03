import os,sys
sys.path.append(os.path.split(os.getcwd())[0])
from hiclib.binnedData import binnedData, binnedDataAnalysis
from hiclib.fragmentHiC import HiCdataset
from mirnylab.numutils import EIG, coarsegrain, project
import numpy
from mirnylab.genome import Genome
import mirnylab.plotting              
import scipy.stats
cr = scipy.stats.spearmanr
import cPickle
from mirnylab.plotting import mat_img,removeAxes,removeBorder

import matplotlib.pyplot as plt 
import matplotlib
#def exampleOfTanayUsage():        
#    a = TanayCorrection()
#    a.loadTanayFragments()
#    a.loadTanayCorrections()
#    a.loadObserved()
#    expected,number = a.calculateTanayExpected(21,22)
#    number
#    mat_img(expected)   #Tanay expected data
#    mat_img(iterativeCorrectWithoutSS(expected))  #Factorized Tanay expected data
#    
#exampleOfTanayUsage()

myGenome = Genome("../../data/hg19",  readChrms = ["#","X"])

GM1M = "../../ErezPaperData/hg19/GM-HindIII-hg19-1M.hm"
GM1MNcoI = "../../ErezPaperData/hg19/GM-NcoI-hg19-1M.hm"
GM200k = "../../ErezPaperData/hg19/GM-HindIII-hg19-200k.hm"
GM200kBreaks = "../../ErezPaperData/hg19/GM-HindIII-hg19-200k-breaks.hm"
GM200kNcoI = "../../ErezPaperData/hg19/GM-NcoI-hg19-200k.hm"
GMFrag = "../../ErezPaperData/hg19/GM-NcoI-hg19_refined.frag"
workingFile1 = "../../ErezPaperData/working1"
workingFile2 = "../../ErezPaperData/working2"
workingFile3 = "../../ErezPaperData/working3"


def correctedScalingPlot():
    
    plt.figure(figsize = (4,4))
    Tanay = binnedDataAnalysis(200000,genome = myGenome)
    Tanay.simpleLoad(GM200kBreaks,"GM-all")    
    Tanay.removePoorRegions()
    Tanay.removeDiagonal()
    Tanay.plotScaling("GM-all", label = "Raw data",color = "#A7A241")
    Tanay.iterativeCorrectWithSS()
    Tanay.plotScaling("GM-all",label = "Corrected",color = "#344370")
    ax = plt.gca() 
    mirnylab.plotting.removeAxes()
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
    
    Tanay = binnedDataAnalysis(resolution = 1000000, genome = myGenome)    
    Tanay.simpleLoad(GM1M,"GM-all")
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()
    Tanay.fakeCis()     
    Tanay.iterativeCorrectWithoutSS()
    plt.figure(figsize = (3.6,3.6))
    Tanay.averageTransMap()

    mirnylab.plotting.removeBorder()
    cb = plt.colorbar(orientation = "vertical")
    cb.set_ticks([-0.05,0.05,0.15])
    for xlabel_i in cb.ax.get_xticklabels(): xlabel_i.set_fontsize(6)    
    mirnylab.plotting.niceShow()    
#doArmPlot()     

    

def checkCorrelationBetweenInterchromosomal(): 
    resolution = 1000000                    
    Tanay = binnedDataAnalysis(resolution,genome = myGenome)     
    Tanay.loadGC()
    Tanay.simpleLoad(GM1M,"GM-all")
    Tanay.simpleLoad(GM1MNcoI,"GM-NcoI")
    Tanay.removePoorRegions()
    Tanay.removeDiagonal()
    Tanay.removeZeros()    
    t =  Tanay.interchromosomalValues()
    t2 = Tanay.interchromosomalValues("GM-NcoI")
    print "raw: ",numpy.corrcoef(t,t2) 
    
    #Tanay.iterativeCorrectWithSS(M = 55)
    Tanay.correct() 
    
    t =  Tanay.interchromosomalValues()
    t2 = Tanay.interchromosomalValues("GM-NcoI")
    print "single correct, ",numpy.corrcoef(t,t2) 
    
    Tanay.iterativeCorrectWithoutSS()
    
    t =  Tanay.interchromosomalValues()
    t2 = Tanay.interchromosomalValues("GM-NcoI")
    print "iterative correct, ",numpy.corrcoef(t,t2) 
    

    Tanay.fakeCis()
    Tanay.iterativeCorrectWithoutSS()
    t =  Tanay.interchromosomalValues()
    t2 = Tanay.interchromosomalValues("GM-NcoI")
    print "trans only iterative correct, ",numpy.corrcoef(t,t2)


#checkCorrelationBetweenInterchromosomal()

def averagePC2TransMap():
    "plot average trans map reconstructed from PC2 only"    
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


#averagePC2TransMap()




def compareInterarmMaps():
    "plots supplementary figure witn 8 inetrarm maps"
    Tanay = binnedDataAnalysis(1000000,myGenome)
    #Tanay.loadK562()
    #Tanay.loadRachel1MB()
    #Tanay.loadRachel1MB("Rachel-1M-breaks", "Rachel-breaks")
    Tanay.simpleLoad(GM1M,"GM-all")
    Tanay.simpleLoad(GM1MNcoI,"GM-NcoI")    
    Tanay.removeDiagonal()
    Tanay.removePoorRegions(cutoff = 1)
    Tanay.removeStandalone(3)
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

#compareInterarmMaps()


def plotDiagonalCorrelation():
    "main paper plot with correlations of diagonals"
    S = 50
    x = numpy.arange(2,S)
    Tanay = binnedData(200000,myGenome)
    Tanay.simpleLoad(GM200k,"GM-all")
    Tanay.simpleLoad(GM200kNcoI,"GM-NcoI")
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()    
    Tanay.removeZeros()
    cors = []
    for i in x:
        cors.append(cr(
                       numpy.diagonal(Tanay.dataDict["GM-all"],i),
                       numpy.diagonal(Tanay.dataDict["GM-NcoI"],i)
                       )[0])
    
    Tanay.iterativeCorrectWithSS(M = 1)
    cors2 = []
    for i in x:
        cors2.append(cr(
                       numpy.diagonal(Tanay.dataDict["GM-all"],i),
                       numpy.diagonal(Tanay.dataDict["GM-NcoI"],i)
                       )[0])
    del Tanay
    Tanay = binnedData(200000,myGenome)
    Tanay.simpleLoad(GM200k,"GM-all")
    Tanay.simpleLoad(GM200kNcoI,"GM-NcoI")
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()
    Tanay.removeZeros()
    Tanay.iterativeCorrectWithSS() 
    cors3 = []
    for i in x:
        cors3.append(cr(
                       numpy.diagonal(Tanay.dataDict["GM-all"],i),
                       numpy.diagonal(Tanay.dataDict["GM-NcoI"],i)
                       )[0])


    matplotlib.rcParams['font.sans-serif']='Arial'
    

    #plt.figure(figsize = (2.3,1.8))
    plt.figure(figsize = (4,4))
    ax = plt.gca()
    fs = 8
    for xlabel_i in ax.get_xticklabels(): xlabel_i.set_fontsize(fs)
    for xlabel_i in ax.get_yticklabels(): xlabel_i.set_fontsize(fs)
    plt.plot(x/5.,cors3,color = "#E5A826",label = "Iterative")
    plt.plot(x/5.,cors2,color = "#28459A",label = "Single")
    plt.plot(x/5.,cors,color = "#E55726",label = "Raw")
    plt.xlabel("Genomic Separation, MB",fontsize = 8)
    plt.ylabel("Spearman correlation",fontsize = 8)
    plt.legend()
                    
    legend = plt.legend(prop={"size":6},loc = 9,handlelength=2)
    legend.draw_frame(False)
    plt.ylim((0,0.6))
    removeAxes(shift = 0)
    
    plt.show()     
     
#plotDiagonalCorrelation()




def plotCrossValidation():
    "main figure subplot with corss-validation"
    matplotlib.rcParams['font.sans-serif']='Arial'
    plt.figure(figsize = (1,1))
    FG = HiCdataset(workingFile1,myGenome)
    FG.load(GMFrag)
    mask = numpy.random.random()< 0.1 
    
    
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
    Tanay.loadGC()
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
    mirnylab.plotting.niceShow()
      
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

 
#doCartoonPlot()
