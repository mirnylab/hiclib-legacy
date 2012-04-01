import base
base #Eclipse warning damper 

import plotting 
import dnaUtils 
import numutils  
import numpy 
na = numpy.array 
from numpy import r_ 

import matplotlib.pyplot as plt 

import cPickle
from math import exp


from binnedData  import binnedData
from expression import strangeGradientDescent
#os.chdir("/old/home/magus/HiC2011")


def oneScaling(datas,centromeres, mode = "within" ,label = "label",color = "r"):

     
    bins = numutils.logbins(2,max([len(i) for i in datas]),1.17)
    observed = [0 for i in bins[:-1]]
    expected = [0 for i in bins[:-1]]
    for x,centromere in zip(datas,centromeres):
        x = x.copy()     
        s = numpy.sum(x,axis = 0) > 0 
        mask = s[:,None] * s[None,:]
        if mode == "within": 
            x[:centromere,:centromere] = 0
            mask[:centromere,:centromere] = 0
            x[centromere:,centromere:]= 0 
            mask[centromere:,centromere:] = 0
        if mode == "between":
            x[:centromere,centromere:] = 0
            mask[:centromere,centromere:] = 0 
            x[centromere:,:centromere] = 0 
            mask[centromere:,:centromere] = 0  
        for i in xrange(len(bins) - 1):
            low = bins[i]
            high = bins[i+1]
            obs = 0
            exp = 0
            high2 = min(high,len(x))
            for k in xrange(low,high2):
                obs += numpy.sum(numpy.diag(x,k))
                exp += numpy.sum(numpy.diag(mask,k))                    
            observed[i] += (obs)
            expected[i] += (exp)
        observed = na(observed,float)
        expected = na(expected,float)
        print observed 
        print expected 
        values = observed/(expected + 0.00001)
        bins = na(bins,float)
        bins2 = 0.5 * (bins[:-1] + bins[1:])

        

    plt.plot(bins2,values,color, linewidth = 2)
        

#
#TR = binnedData(200000,"HG18")
#chroms = numpy.array(TR.chromosomeStarts)
#numpy.savetxt("ChromosmeStartPositions-hg18-200k",chroms)
#TR = binnedData(1000000,"HG18")
#chroms = numpy.array(TR.chromosomeStarts)
#numpy.savetxt("ChromosmeStartPositions-hg18-1M",chroms)
#exitProgram()

#TR = binnedData(200000,"HG18")
#TR.simpleLoad("GM-all-hg18-200k-raw.npz","GM-all")
#numpy.savetxt("GM-raw-HindIII-200k-hg18-SS",TR.singlesDict["GM-all"])
#
#TR = binnedData(200000,"HG18")
#TR.simpleLoad("GM-NcoI-hg18-200k-raw.npz","GM-all")
#numpy.savetxt("GM-raw-NcoI-200k-hg18-SS",TR.singlesDict["GM-all"])
#
#TR = binnedData(1000000,"HG18")
#TR.simpleLoad("GM-all-hg18-1M-raw.npz","GM-all")
#numpy.savetxt("GM-raw-HindIII-1M-hg18-SS",TR.singlesDict["GM-all"])
#
#TR = binnedData(1000000,"HG18")
#TR.simpleLoad("GM-NcoI-hg18-1M-raw.npz","GM-all")
#numpy.savetxt("GM-raw-NcoI-1M-hg18-SS",TR.singlesDict["GM-all"])
#
#
#
TR = binnedData(1000000,"HG18")
#TR.simpleLoad("GM-all-hg18-200k.npz","GM-all")
#
#

TR.simpleLoad("GM-NcoI-hg18-1M.npz","GM-all")
TR.removePoorRegions()
#numpy.savetxt("GM-original-HindIII-1M-hg18-SS",TR.singlesDict["GM-all"])
mapping = TR.removeZeros()
TR.fancyFakeCis()
TR.doEig()
EIG = TR.EIG["GM-all"]
E1 = EIG[:,0]
E2 = EIG[:,1]
E1new = numpy.zeros(len(mapping),float) * numpy.NAN
E2new = numpy.zeros(len(mapping),float) * numpy.NAN
E1new[mapping] = E1
E2new[mapping] = E2
numpy.savetxt("GM-NcoI-Eig1-1M",E1new)
numpy.savetxt("GM-NcoI-Eig2-1M",E2new)
exit()

numpy.savetxt("GM-HindIII-1M-Eig1")
numpy.savetxt("GM-corrected-HindIII-1M-hg18-SS",TR.singlesDict["GM-all"])
chroms = numpy.array(TR.chromosomeStarts)
numpy.savetxt("GM-HindIII-chromosmeStartPositions-hg18",chroms)

TR = binnedData(200000,"HG18")
TR.simpleLoad("GM-NcoI-hg18-200k.npz","GM-all")
TR.removePoorRegions()
numpy.savetxt("GM-original-NcoI-200k-hg18-SS",TR.singlesDict["GM-all"])
TR.ultracorrectAll()
numpy.savetxt("GM-corrected-NcoI-200k-hg18-SS",TR.singlesDict["GM-all"])
chroms = numpy.array(TR.chromosomeStarts)
numpy.savetxt("GM-NcoI-chromosmeStartPositions-hg18",chroms)


TR = binnedData(1000000,"HG18")
TR.simpleLoad("GM-NcoI-hg18-1M.npz","GM-all")
TR.removePoorRegions()
numpy.savetxt("GM-original-NcoI-1M-hg18-SS",TR.singlesDict["GM-all"])
TR.ultracorrectAll()
numpy.savetxt("GM-corrected-NcoI-1M-hg18-SS",TR.singlesDict["GM-all"])
chroms = numpy.array(TR.chromosomeStarts)
numpy.savetxt("GM-NcoI-chromosmeStartPositions-hg18",chroms)


exit()

def drawToyModel(M = 14):
    "draws cartoon of iterative correction of size M"    
    d = numpy.random.random((M,M)) * 1
    for i in xrange(M):
        for j in xrange(i,M):
            d[i,j] = d[j,i]
    a = d      
    a += numpy.diag(numpy.ones(M)) * 1.5
    a += numpy.diag(numpy.ones(M))[::-1] * 1.5    
    bias = numpy.random.random(M) * 3 + 0.35   #strength of the biases is here

    ax = plt.subplot(2,3,1)
    plotting.removeBorder(ax = ax)
    plt.imshow(a.copy(),interpolation = "nearest")
    plt.title("Initial signal")
    a *= bias[:,None]
    a *= bias[None,:] 

    ax = plt.subplot(2,3,2)
    plt.title("With mult biases")
    plotting.removeBorder(ax = ax)
    plt.imshow(a.copy(),interpolation = "nearest")
    
    ax = plt.subplot(2,3,3)
    plt.title("Biases extracted by\n iterative correction")
    plotting.removeBorder(ax = ax)
    pa = a/numutils.ultracorrect(a)    
    plt.imshow(pa,interpolation = "nearest")    
    
    ax = plt.subplot(2,3,4)
    plt.title("Iterative corrected")
    plotting.removeBorder(ax = ax)
    plt.imshow(numutils.ultracorrect(a),interpolation = "nearest")
    
    ax = plt.subplot(2,3,5)
    plt.title("Corrected by sum")
    plotting.removeBorder(ax = ax)
    plt.imshow(numutils.correct(a),interpolation = "nearest")
    
    ax = plt.subplot(2,3,6)
    plt.title("Biases")
    plotting.removeBorder(ax = ax)
    plt.plot(r_[bias,bias[-1]],drawstyle = "steps-post")
    
    plt.show()

def saveDataForDomains():
    TR = binnedData(genomeType="HG19",resolution = 200000)
    TR.simpleLoad("GM-all-hg19-200k.npz","GM-all")
    TR.simpleLoad("GM-NcoI-hg19-200k.npz","GM-NcoI")
    TR.removePoorRegions()
    TR.removeStandalone(5)
    mapping = TR.removeZeros()
    TR.fancyFakeCisOld()
    numpy.savez_compressed("200k_analyze",**{"mapping":mapping,"data":TR.dataDict["GM-all"]})
    numpy.savez_compressed("200k_analyze_NcoI",**{"mapping":mapping,"data":TR.dataDict["GM-NcoI"]})

#saveDataForDomains()
    
#numpy.savez_compressed("map1",mymap)

#mymap = numpy.load("map1.npz")["arr_0"]

  

def doDomains(filename):
    data = dict(numpy.load(filename))
    mymap = data["data"]

    mymap = numpy.array(mymap,dtype = "int64",order = "C")         
    GD = strangeGradientDescent(mymap)
    GD.doSearch()
    return GD.vec
 
#dom = doDomains("200k_analyze_NcoI.npz")
#cPickle.dump(dom,open("200kDomNcoI",'wb'))

#dom = doDomains("200k_analyze.npz")
#cPickle.dump(dom,open("200kDom",'wb'))


def analyzeEigenvector():
    "does domain analysis for Magda probes"
    positions = """17    5973934     6027745         wscd1
17    13972846    14111994    8.04    cox10
17    21279699    21320404    7.25    kcnj12
17    29718642    29865236    8.49    rab11fip4
17    37844393    37884915    8.07    erbb2
17    45963516    46016323    8.12    sp2
17    53818356    53864748    7.85    pctp
17    62010914    62055278    8.19    scn4a
17    70097161    70142561    8.08    sox9
17    78010442    78074412    7.92    ccdc40"""
    positions = [(int(i.split()[1]) + int(i.split()[2]))/2 for i in positions.split("\n")]
    mapping = numpy.load("200k_analyze.npz")["mapping"]
    vector = cPickle.load(open("200kDom",'rb'))
    #vector2 = cPickle.load(open("200kDom"))
    vector = vector[len(vector)/2:]
    print len(mapping), len(vector)
    TR = binnedData(200000,"HG19")
    TR.loadGC()    
    allGenome = numpy.zeros(len(mapping),float) - 3    
    allGenome[mapping] = vector
    #------------------
    allGenome2 = TR.trackDict["GC"] 
    #-----------------    
    HG19 = dnaUtils.Genome("HG19")
    HG19.createMapping(200000)
    ch17 = allGenome[HG19.chromosomeStarts[16]:HG19.chromosomeEnds[16]]
    ch17gc = allGenome2[HG19.chromosomeStarts[16]:HG19.chromosomeEnds[16]]    
    locations = na(positions)/200000
    values = ch17[locations]
    svector = numpy.sort(ch17)
    inds = numpy.searchsorted(svector,values) / float(len(ch17))
    print inds 
    p1 = numpy.percentile(ch17,5)
    p3 = numpy.percentile(ch17,35)
    p7 = numpy.percentile(ch17,65)
    p9 = numpy.percentile(ch17,95)

    p1g = numpy.percentile(ch17gc,5)
    p3g = numpy.percentile(ch17gc,35)
    p7g = numpy.percentile(ch17gc,65)
    p9g = numpy.percentile(ch17gc,95)
    
    mask = numpy.zeros(len(ch17),bool)
    mask[(ch17>p1) * (ch17 < p3) * (ch17gc > p1g) * (ch17gc < p3g)] = 1
    mask[(ch17>p7) * (ch17 < p9) * (ch17gc > p7g) * (ch17gc < p9g)] = 1
    print mask 
    
    r = numutils.rank(ch17)
    x = numpy.arange(len(ch17))*.2+.1
    plt.plot(x,r,color = 'k',label = "rank of eigenvector")
    plt.plot(x[mask],r[mask],'.',markersize = 10,color =  "#88E526",label = "permitted regions")
    plt.scatter([i/1000000 for i in positions],[len(ch17)/2 for i in positions],label = "old probes",color = "#28459A")
    positions2 = [4.6,14,19,22.1,27,32.5,37,39.3,45.9,53,62,68.5,75]
    plt.scatter(positions2,[1.2*len(ch17)/2 for i in positions2],color ="#E55726" ,label = "new probes")
    plt.xlabel("position")
    plt.ylabel("rank")
    plotting.niceShow()
    
    r[(ch17>p3) * (ch17 < p7)] = len(ch17)/2
    print ch17<p1
    r[ch17<p1] = min(r) 
    r[ch17>p9] = max(r) 
    plt.plot(numpy.arange(len(ch17))*.2+.1,r,label = 'eigenvector') 
    plt.scatter([i/1000000 for i in positions],[len(ch17)/2 for i in positions],label = "old probes")
    plt.show()

    plt.hist(allGenome,bins = 30)
    plt.show()
          
#analyzeEigenvector()
