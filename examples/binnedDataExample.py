import os,sys
sys.path.append(os.path.split(os.getcwd())[0])
from binnedData import binnedData, binnedDataAnalysis
from mirnylab.numutils import *
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
#    mat_img(ultracorrect(expected))  #Factorized Tanay expected data
#    
#exampleOfTanayUsage()

myGenome = Genome("../../data/hg19",  readChrms = ["#","X"])

GM1M = "../../ErezPaperData/hg19/GM-HindIII-hg19-1M"
GM1MNcoI = "../../ErezPaperData/hg19/GM-NcoI-hg19-1M"
GM200k = "../../ErezPaperData/hg19/GM-HindIII-hg19-200k"
GM200kBreaks = "../../ErezPaperData/hg19/GM-HindIII-hg19-200k-breaks"
GM200kNcoI = "../../ErezPaperData/hg19/GM-NcoI-hg19-200k"

def correctedScalingPlot():
    
    plt.figure(figsize = (4,4))
    Tanay = binnedDataAnalysis(200000,genome = myGenome)
    Tanay.simpleLoad(GM200kBreaks,"GM-all")    
    Tanay.removePoorRegions()
    Tanay.plotScaling("GM-all", label = "Raw data",color = "#A7A241")
    Tanay.ultracorrectAll()
    Tanay.plotScaling("GM-all",label = "Corrected",color = "#344370")
    ax = plt.gca() 
    mirnylab.plotting.removeAxes(shift = 0 )
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
    Tanay.ultracorrect()
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
    Tanay.removeZeros()    
    t =  Tanay.interchromosomalValues()
    t2 = Tanay.interchromosomalValues("GM-NcoI")
    print "raw: ",numpy.corrcoef(t,t2) 
    
    
    
    #Tanay.ultracorrectAll(M = 55)
    Tanay.correct() 
    
    t =  Tanay.interchromosomalValues()
    t2 = Tanay.interchromosomalValues("GM-NcoI")
    print "single correct, ",numpy.corrcoef(t,t2) 
    
    Tanay.ultracorrect()
    
    t =  Tanay.interchromosomalValues()
    t2 = Tanay.interchromosomalValues("GM-NcoI")
    print "iterative correct, ",numpy.corrcoef(t,t2) 
    
    Tanay.ultracorrectAll()
    t =  Tanay.interchromosomalValues()
    t2 = Tanay.interchromosomalValues("GM-NcoI")
    print "iterative + sm correct, ",numpy.corrcoef(t,t2) 
    
    
    Tanay.fakeCis()
    Tanay.ultracorrect()
    t =  Tanay.interchromosomalValues()
    t2 = Tanay.interchromosomalValues("GM-NcoI")
    print "trans only iterative correct, ",numpy.corrcoef(t,t2)


#checkCorrelationBetweenInterchromosomal()


def averagePC2TransMap():
    "plot average trans map reconstructed from PC2 only"    
    Tanay = binnedData(resolution = 1000000, genome = myGenome)
    #Tanay.loadNcoI()
    Tanay.simpleLoad(GM1M,"GM-all")    
    Tanay.removePoorRegions() 
    s = Tanay.removeZeros()
    arg = numpy.nonzero(s)[0]
    Tanay.TrunkTrans()
    Tanay.fakeCis()
    Tanay.ultracorrect()
    Tanay.doPCA()
    PC2 = Tanay.doEig()["GM-all"][:,1]
    c = numpy.zeros(arg.max() +1)
    c[arg] = PC2

    proj = c[:,None] * c[None,:]
    mm = c != 0 
    mask = mm[:,None] * mm[None,:]
    proj[mask] += (2 * abs(proj.min()) ) 
    Tanay = binnedData() 
    Tanay.dataDict["BLA"] = proj
    Tanay.TrunkTrans(0.001)     
    Tanay.averageTransMap("BLA",chromMax = 22)
    plt.show()
#averagePC2TransMap()


def doTanayDomainsTanayDataPlot():
    matplotlib.rcParams['font.sans-serif']='Arial'
    Tanay = binnedData()
    Tanay.loadTanaySmoothed()
    Tanay.dataDict.pop("TanayExpected")
    Tanay.dataDict["GM-all"] = Tanay.dataDict.pop("TanayObserved")
    "load GM-all data here to check Tanay domains vs our eigenvector"
    #Tanay.simpleLoad("GM-all-1M","GM-all")
    Tanay.loadGC()
    Tanay.loadTanayDomains()
    Tanay.removePoorRegions()
    Tanay.removeZeros()
    Tanay.fancyFakeCis()
    Tanay.ultracorrect() 
    Tanay.doEig()
    PC = Tanay.EIG["GM-all"]
    PC1 = PC[:,0]
    PC2 = PC[:,1]
    if PC1[0] <0: 
        PC1 = -PC1
    if PC2[0] >0:
        PC2 = -PC2

    colors = ["#28459A","#E5A826","#E55726"]
    labels = ["Active","Centromere proximal","Centromere distal"]
    
    #plt.figure(figsize = (1.8,1.8))
    plt.figure(figsize = (4,4))
    fs1 = 6
    fs2 = 8 
    
    fs1 = 10
    fs2 = 12
    def onePlot(tan,PC1,PC2,overrideColor = None):
        ax = plt.gca()
        for xlabel_i in ax.get_xticklabels(): xlabel_i.set_fontsize(fs1)  
        for xlabel_i in ax.get_yticklabels(): xlabel_i.set_fontsize(fs1)
        #plt.xticks([0,0.05,0.1])
        #plt.yticks([-0.1,-0.05,0,0.05])
        
        P1 = []
        P2 = []
        C1 = []
    
        for i in xrange(3):
            mask = tan == i
            print mask.sum()
            P1+= list(PC1[mask])
            P2+= list(PC2[mask])
            for j in xrange(mask.sum()):
                C1.append(colors[i])
            plt.scatter([],[],c = colors[i],label = labels[i], linewidth = 0,s=3)
        inds = numpy.arange(len(P1))
        numpy.random.shuffle(inds)
        P1 = [P1[i] for i in inds]
        P2 = [P2[i] for i in inds]
        C = [C1[i] for i in inds]
        if overrideColor == True: C = "grey"
        plt.scatter(P1,P2,c = C,linewidth = 0,s=2)
        if overrideColor == None: 
            legend = plt.legend(loc=0,prop={"size":fs1},handlelength = 4)
            legend.draw_frame(False)
        removeAxes(shift = 0)

    tan = Tanay.trackDict["TanayDomains"]
    plt.subplot(211)
    plt.ylim((-0.06,0.06))
    plt.xlim((-0.05,0.2))
    onePlot(tan,PC1,PC2)
    plt.xlabel("Eig1",fontsize = fs1)
    plt.ylabel("Eig2",fontsize = fs1)
    plt.subplot(212)
    plt.ylim((-0.06,0.06))
    plt.xlim((-0.05,0.2))
    onePlot(tan,PC1,PC2,True)
    plt.xlabel("Eig1",fontsize = fs1)
    plt.ylabel("Eig2",fontsize = fs1)
    
    #plt.xlabel("PC1",fontsize = 6)
    #plt.ylabel("PC2",fontsize = 6)
    plt.gcf().subplots_adjust(left=0.07, bottom=0.12, top=0.98, right=0.98)
    plt.show()

    raise 
    
#doTanayDomainsTanayDataPlot()


def compareInterarmMaps():
    "plots supplementary figure witn 8 inetrarm maps"
    Tanay = binnedData()
    #Tanay.loadK562()
    #Tanay.loadRachel1MB()
    #Tanay.loadRachel1MB("Rachel-1M-breaks", "Rachel-breaks")
    Tanay.simpleLoad("../data/GM-all-hg18-1M","GM-all")
    Tanay.simpleLoad("../data/GM-NcoI-hg18-1M","GM-NcoI")
    Tanay.removeStandalone(0)
    Tanay.removePoorRegions()
    fs = 10
    vmin = None
    vmax = None
    plt.subplot(421)
    plt.title("GM, HindIII, raw", fontsize = fs)
    Tanay.averageTransMap("GM-all",vmin = vmin, vmax = vmax)
    plt.subplot(422)
    plt.title("GM, NcoI, raw", fontsize = fs)
    Tanay.averageTransMap("GM-NcoI",vmin = vmin, vmax = vmax)
    Tanay.ultracorrectAll()
    vmin = None
    vmax = None
    plt.subplot(425)
    plt.title("GM, HindIII, with SS reads", fontsize = fs)
    m1 = Tanay.averageTransMap("GM-all",vmin = vmin, vmax = vmax)
    plt.subplot(426)
    plt.title("GM, NcoI, with SS reads", fontsize = fs)
    m2 = Tanay.averageTransMap("GM-NcoI",vmin = vmin, vmax = vmax)
    
    Tanay.ultracorrect()
    vmin = None
    vmax = None
    
    plt.subplot(423)
    plt.title("GM, HindIII, no SS reads", fontsize = fs)
    Tanay.averageTransMap("GM-all",vmin = vmin, vmax=vmax)
    plt.subplot(424)
    plt.title("GM, NcoI, no ss reads", fontsize = fs)
    Tanay.averageTransMap("GM-NcoI",vmin = vmin, vmax=vmax)
    
    Tanay.fancyFakeCis()
    plt.subplot(427)
    plt.title("GM, HindIII, trans only", fontsize = fs)
    Tanay.averageTransMap("GM-all",vmin = vmin, vmax=vmax)
    plt.subplot(428)
    plt.title("GM, NcoI, trans only", fontsize = fs)
    Tanay.averageTransMap("GM-NcoI",vmin = vmin, vmax=vmax)
    
    
    


    plt.show()

#compareInterarmMaps()


def plotDiagonalCorrelation():
    "main paper plot with correlations of diagonals"
    S = 50
    x = numpy.arange(2,S)
    Tanay = binnedData(200000)
    Tanay.simpleLoad("../data/GM-all-hg18-200k","GM-all")
    Tanay.simpleLoad("../data/GM-NcoI-hg18-200k","GM-NcoI")
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()    
    Tanay.removeZeros()
    cors = []
    for i in x:
        cors.append(cr(
                       numpy.diagonal(Tanay.dataDict["GM-all"],i),
                       numpy.diagonal(Tanay.dataDict["GM-NcoI"],i)
                       )[0])
    
    Tanay.ultracorrectAll(M = 1)
    cors2 = []
    for i in x:
        cors2.append(cr(
                       numpy.diagonal(Tanay.dataDict["GM-all"],i),
                       numpy.diagonal(Tanay.dataDict["GM-NcoI"],i)
                       )[0])
    del Tanay
    Tanay = binnedData(200000)
    Tanay.simpleLoad("../data/GM-all-hg18-200k","GM-all")
    Tanay.simpleLoad("../data/GM-NcoI-hg18-200k","GM-NcoI")
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()
    Tanay.removeZeros()
    Tanay.ultracorrectAll() 
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
     
plotDiagonalCorrelation()


def plotCrossValidation():
    "main figure subplot with corss-validation"
    matplotlib.rcParams['font.sans-serif']='Arial'
    plt.figure(figsize = (1,1))
    Tanay = binnedData(1000000)    
    Tanay.simpleLoad("GM-all-10p","GM-1")   #need to create these datasets using fragment-level analysis
    Tanay.simpleLoad("GM-all-90p","GM-9")
    Tanay.removePoorRegions()
    Tanay.ultracorrectAll()
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
    Tanay.fancyFakeCis()        
    Tanay.ultracorrect()
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
    ax = plt.gca()
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
