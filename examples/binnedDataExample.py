"""
This is a collection of utils to reconstruct almost any plot in the paper.
It requires pre-processed datasets for every plot.
Most of the filenames are defined within the text,
but you can redefine them in the beginning.

Please contact imakaev@mit.edu with all questions.
"""
import os
import sys
from mirnylib.systemutils import setExceptionHook
sys.path.append(os.path.split(os.getcwd())[0])
from hiclib.binnedData import binnedData, binnedDataAnalysis,\
    experimentalBinnedData

import mirnylib.systemutils
mirnylib.systemutils
from hiclib.fragmentHiC import HiCdataset
from mirnylib.numutils import EIG, coarsegrain, project, arrayInArray
import numpy
import mirnylib.plotting
import scipy.stats
import scipy.ndimage
from hiclib import fragmentHiC
cr = scipy.stats.spearmanr
import cPickle
from mirnylib.plotting import mat_img, removeAxes, removeBorder, niceShow

import matplotlib.pyplot as plt
import matplotlib

genomeVersion = "hg18"
myGenome = "../../data/%s" % genomeVersion

GM1M = "../../ErezPaperData/%s/GM-HindIII-%s-1M.hm" % (
    genomeVersion, genomeVersion)
GM1MNcoI = "../../ErezPaperData/%s/GM-NcoI-%s-1M.hm" % (
    genomeVersion, genomeVersion)
GM200k = "../../ErezPaperData/%s/GM-HindIII-%s-200k.hm" % (
    genomeVersion, genomeVersion)
GM200kBreaks = "../../ErezPaperData/%s/GM-HindIII-%s-200k-breaks.hm" % (
    genomeVersion, genomeVersion)
GM1MBreaks = "../../ErezPaperData/%s/GM-HindIII-%s-1M-breaks.hm" % (
    genomeVersion, genomeVersion)
GM200kNcoI = "../../ErezPaperData/%s/GM-NcoI-%s-200k.hm" % (
    genomeVersion, genomeVersion)
GMFrag = "../../ErezPaperData/%s/GM-NcoI-%s_refined.frag" % (
    genomeVersion, genomeVersion)

tcc1M = "../../tcc/%s/tcc-HindIII-%s-1M.hm" % (genomeVersion, genomeVersion)
tcc200k = "../../tcc/%s/tcc-HindIII-%s-200k.hm" % (
    genomeVersion, genomeVersion)
workingFile1 = "../../ErezPaperData/working1"
workingFile2 = "../../ErezPaperData/working2"
workingFile3 = "../../ErezPaperData/working3"


def correctedScalingPlot():
    "Paper figure to compare scaling before/after correction"
    plt.figure(figsize=(4, 4))
    Tanay = binnedDataAnalysis(200000, genome=myGenome)
    Tanay.simpleLoad(GM200kBreaks, "GM-all")
    Tanay.removePoorRegions()
    Tanay.removeDiagonal()
    Tanay.plotScaling("GM-all", label="Raw data", color="#A7A241")
    Tanay.iterativeCorrectWithSS()
    Tanay.plotScaling("GM-all", label="Corrected", color="#344370")
    ax = plt.gca()
    mirnylib.plotting.removeAxes()
    fs = 6
    plt.xlabel("Genomic distance (MB)", fontsize=6)
    plt.ylabel("Contact probability", fontsize=6)
    for xlabel_i in ax.get_xticklabels():
        xlabel_i.set_fontsize(fs)
    for xlabel_i in ax.get_yticklabels():
        xlabel_i.set_fontsize(fs)
    legend = plt.legend(loc=0, prop={"size": 6})
    legend.draw_frame(False)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()


def doArmPlot(filename=GM1M, genome=myGenome, mouse=False, **kwargs):
    "Plot an single interarm map - paper figure"
    Tanay = binnedDataAnalysis(1000000, genome)
    Tanay.simpleLoad(filename, "GM-all")
    if mouse == True:
        Tanay.fakeTranslocations([(0, 0, None, 12, 52000000, None),
                                  (4, 45000000, None, 12, 0, 30000000),
                                  (9, 0, 50000000, 12, 0, 35000000)])
        Tanay.removeChromosome(19)
    else:
        Tanay.removeChromosome(22)
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()
    Tanay.truncTrans()
    Tanay.fakeCis()
    #mat_img(Tanay.dataDict["GM-all"])
    #plt.figure(figsize = (3.6,3.6))
    Tanay.averageTransMap("GM-all", **kwargs)

    #plotting.removeBorder()
    cb = plt.colorbar(orientation="vertical")
    #cb.set_ticks([-0.05,0.05,0.15])
    for xlabel_i in cb.ax.get_xticklabels():
        xlabel_i.set_fontsize(6)

doArmPlot()


def doReconstructedArmPlot(filename=GM1M, genome=myGenome,
                           usePCs=[0, 1], mouse=False,
                           **kwargs):
    "Plot an PC2-PC3 interarm map - supp paper figure"
    Tanay = binnedDataAnalysis(1000000, genome)
    Tanay.simpleLoad(filename, "GM-all")
    if mouse == True:
        Tanay.fakeTranslocations([(0, 0, None, 12, 52000000, None),
                                  (4, 45000000, None, 12, 0, 30000000),
                                  (9, 0, 50000000, 12, 0, 35000000)])
        Tanay.removeChromosome(19)
    else:
        Tanay.removeChromosome(22)
    Tanay.removeDiagonal()
    Tanay.removePoorRegions()
    Tanay.removeZeros()
    Tanay.truncTrans()
    Tanay.fakeCis()
    Tanay.doEig(numPCs=max(usePCs) + 1)
    print Tanay.eigEigenvalueDict

    Tanay.restoreZeros(value=0)
    PCs = Tanay.EigDict["GM-all"][usePCs]
    eigenvalues = Tanay.eigEigenvalueDict["GM-all"][usePCs]

    proj = reduce(lambda x, y: x + y,
                  [PCs[i][:, None] * PCs[i][None, :] * \
                   eigenvalues[i] for i in xrange(len(PCs))])
    mask = PCs[0] != 0
    mask = mask[:, None] * mask[None, :]
    data = Tanay.dataDict["GM-all"]
    datamean = numpy.mean(data[mask])
    proj[mask] += datamean
    Tanay.dataDict["BLA"] = proj
    Tanay.averageTransMap("BLA", **kwargs)
    cb = plt.colorbar(orientation="vertical")
    for xlabel_i in cb.ax.get_xticklabels():
        xlabel_i.set_fontsize(6)


def differentArmPlots():
    "Different inetrarm maps - supp figure"
    plt.figure(figsize=(6, 6))
    plt.subplot(221)
    doArmPlot("../../ErezPaperData/hg18/GM-HindIII-hg18-1M.hm",
              "../../data/hg18")
    plt.title("Human, Hi-C, Lieberman 2009")
    plt.subplot(222)
    doArmPlot("../../tcc/hg18/tcc-HindIII-hg18-1M.hm", "../../data/hg18")
    plt.title("Human, TCC, Kalhor 2011")
    plt.subplot(223)
    doArmPlot("../../mouse/data/combined/mouse_all-1M.hm",
              "../../data/mm9", mouse=True)
    plt.title("Mouse, Hi-C, McCord 2012")
    plt.show()


def differentArmPlotsWithReconstructedHeatmaps():
    "Many reconstructed heatmaps - supp paper figure"
    plt.figure(figsize=(8, 8))
    plt.subplot(531)
    doArmPlot("../../ErezPaperData/hg18/GM-HindIII-hg18-1M.hm",
              "../../data/hg18", vmin=-0.05, vmax=0.15)
    plt.title("Human HiC 2009")
    plt.subplot(532)
    doArmPlot("../../tcc/hg18/tcc-HindIII-hg18-1M.hm",
              "../../data/hg18", vmin=-0.15, vmax=0.3)
    plt.title("Human TCC 2011")
    plt.subplot(533)
    doArmPlot("../../mouse/data/combined/mouse_all-1M.hm",
              "../../data/mm9", mouse=True, vmin=-0.15, vmax=0.15)
    plt.title("Mouse 2012")

    for num, chromSet in enumerate([[0], [1], [2], [1, 2]]):
        plt.subplot(5, 3, 3 * num + 4)
        doReconstructedArmPlot("../../ErezPaperData/hg18/"\
        "GM-HindIII-hg18-1M.hm",
                               "../../data/hg18", chromSet,
                               vmin=-0.05, vmax=0.15)
        plt.title("From E" + "".join([str(i + 1) +
                                      ", " for i in chromSet])[:-2])
        plt.subplot(5, 3, 3 * num + 5)
        doReconstructedArmPlot("../../tcc/hg18/tcc-HindIII-hg18-1M.hm",
                               "../../data/hg18", chromSet,
                               vmin=-0.15, vmax=0.3)
        plt.title("From E" + "".join([str(i + 1) +
                                      ", " for i in chromSet])[:-2])
        plt.subplot(5, 3, 3 * num + 6)
        doReconstructedArmPlot("../../mouse/data/combined/mouse_all-1M.hm",
                               "../../data/mm9", chromSet, mouse=True,
                                vmin=-0.15, vmax=0.15)
        plt.title("From E" + "".join([str(i + 1) +
                                      ", " for i in chromSet])[:-2])
    plt.show()

#differentArmPlotsWithReconstructedHeatmaps()


def compareInterarmMaps():
    "plots witn 8 inetrarm maps - paper supplement figure"
    Tanay = binnedDataAnalysis(1000000, myGenome)

    Tanay.simpleLoad(GM1M, "GM-all")
    Tanay.simpleLoad(GM1MNcoI, "GM-NcoI")
    Tanay.removeDiagonal()
    Tanay.removePoorRegions(cutoff=2)
    #Tanay.removeStandalone(3)
    fs = 10
    vmin = None
    vmax = None

    plt.subplot(421)
    plt.title("GM, HindIII, raw", fontsize=fs)
    Tanay.averageTransMap("GM-all", vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.subplot(422)
    plt.title("GM, NcoI, raw", fontsize=fs)
    Tanay.averageTransMap("GM-NcoI", vmin=vmin, vmax=vmax)
    plt.colorbar()

    Tanay.iterativeCorrectWithSS()
    vmin = None
    vmax = None

    plt.subplot(425)

    plt.title("GM, HindIII, with SS reads", fontsize=fs)
    Tanay.averageTransMap("GM-all", vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.subplot(426)

    plt.title("GM, NcoI, with SS reads", fontsize=fs)
    Tanay.averageTransMap("GM-NcoI", vmin=vmin, vmax=vmax)
    plt.colorbar()

    Tanay.iterativeCorrectWithoutSS()
    vmin = None
    vmax = None
    plt.subplot(423)
    plt.title("GM, HindIII, no SS reads", fontsize=fs)
    Tanay.averageTransMap("GM-all", vmin=vmin, vmax=vmax)
    plt.colorbar()

    plt.subplot(424)
    plt.title("GM, NcoI, no ss reads", fontsize=fs)
    Tanay.averageTransMap("GM-NcoI", vmin=vmin, vmax=vmax)
    plt.colorbar()
    Tanay.fakeCis()

    vmin = None
    vmax = None
    plt.subplot(427)

    plt.title("GM, HindIII, trans only", fontsize=fs)
    Tanay.averageTransMap("GM-all", vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.subplot(428)

    plt.title("GM, NcoI, trans only", fontsize=fs)
    Tanay.averageTransMap("GM-NcoI", vmin=vmin, vmax=vmax)
    plt.colorbar()

    plt.show()


def compareCorrelationOfEigenvectors():
    """Plot correlation figure with eigenvector correlation between datasets
    paper figure """
    Tanay = binnedDataAnalysis(1000000, "../../data/hg18")
    Tanay.simpleLoad("../../ErezPaperData/hg18/GM-HindIII-hg18-1M.hm", "Erez")
    Tanay.simpleLoad("../../ErezPaperData/hg18/GM-NcoI-hg18-1M.hm", "NcoI")
    Tanay.simpleLoad("../../tcc/hg18/tcc-HindIII-hg18-1M.hm", "TCC")
    Tanay.removeDiagonal()
    Tanay.removePoorRegions()
    Tanay.removeZeros()
    Tanay.truncTrans()
    Tanay.fakeCis()
    M = 10
    Tanay.doEig(numPCs=M)

    E1 = Tanay.EigDict["Erez"]
    E2 = Tanay.EigDict["NcoI"]
    E3 = Tanay.EigDict["TCC"]

    data = numpy.zeros((M, M))
    data2 = numpy.zeros((M, M))
    data3 = numpy.zeros((M, M))
    for i in xrange(M):
        for j in xrange(M):
            data[i][j] = abs(numpy.corrcoef(E2[i], E1[j])[0, 1])
            data2[i][j] = abs(numpy.corrcoef(E3[i], E1[j])[0, 1])
            data3[i][j] = abs(numpy.corrcoef(E3[i], E2[j])[0, 1])
    plt.figure(figsize=(7.5, 2.5))
    plt.gcf().subplots_adjust(0.2, 0.2, 0.85, 0.85)
    plt.subplot(131)
    plt.xlabel("HiC 2009, HindIII")
    plt.ylabel("HiC 2009, NcoI")
    #plt.title("Abs. correlation between eigenvectors")
    plt.imshow(data, interpolation="nearest", vmin=0, vmax=1)
    plt.colorbar()
    plt.subplot(132)
    plt.xlabel("HiC 2009, HindIII")
    plt.ylabel("TCC 2011, HindIII")
    #plt.title("Abs. correlation between eigenvectors")
    plt.imshow(data2, interpolation="nearest", vmin=0, vmax=1)
    plt.colorbar()
    plt.subplot(133)
    plt.xlabel("HiC 2009, NcoI")
    plt.ylabel("TCC 2011, HindIII")
    #plt.title("Abs. correlation between eigenvectors")
    plt.imshow(data3, interpolation="nearest", vmin=0, vmax=1)
    plt.colorbar()
    plt.show()
    raise


def plotDiagonalCorrelation():
    "Correlation of diagonal bins - paper figure"
    S = 50
    x = numpy.arange(2, S)
    Tanay = binnedData(200000, myGenome)
    Tanay.simpleLoad(GM200k, "GM-HindIII")
    Tanay.simpleLoad(GM200kNcoI, "GM-NcoI")
    Tanay.simpleLoad(tcc200k, "TCC")
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()
    Tanay.removeZeros()
    pairs = [("GM-HindIII", "GM-NcoI"), ("GM-HindIII", "TCC"), (
        "GM-NcoI", "TCC")]
    cors = [[] for _ in pairs]
    for i in x:
        for j, pair in enumerate(pairs):
            cors[j].append(cr(
                           numpy.diagonal(Tanay.dataDict[pair[0]], i),
                           numpy.diagonal(Tanay.dataDict[pair[1]], i)
                           )[0])

    Tanay.iterativeCorrectWithoutSS(M=1)
    cors2 = [[] for _ in pairs]
    for i in x:
        for j, pair in enumerate(pairs):
            cors2[j].append(cr(
                            numpy.diagonal(Tanay.dataDict[pair[0]], i),
                            numpy.diagonal(Tanay.dataDict[pair[1]], i)
                            )[0])
    Tanay.iterativeCorrectWithoutSS(M=20)
    cors3 = [[] for _ in pairs]
    for i in x:
        for j, pair in enumerate(pairs):
            cors3[j].append(cr(
                            numpy.diagonal(Tanay.dataDict[pair[0]], i),
                            numpy.diagonal(Tanay.dataDict[pair[1]], i)
                            )[0])

    matplotlib.rcParams['font.sans-serif'] = 'Arial'

    #plt.figure(figsize = (2.3,1.8))
    print cors
    print cors2
    print cors3
    plt.figure(figsize=(10, 3))
    ax = plt.gca()
    for j, pair in enumerate(pairs):
        plt.subplot(1, len(pairs), j)
        fs = 8
        for xlabel_i in ax.get_xticklabels():
            xlabel_i.set_fontsize(fs)
        for xlabel_i in ax.get_yticklabels():
            xlabel_i.set_fontsize(fs)
        plt.title("%s vs %s" % pair)
        plt.plot(x / 5., cors3[j], color="#E5A826", label="Iterative")
        plt.plot(x / 5., cors2[j], color="#28459A", label="Single")
        plt.plot(x / 5., cors[j], color="#E55726", label="Raw")
        plt.xlabel("Genomic Separation, MB", fontsize=8)
        plt.ylabel("Spearman correlation", fontsize=8)
        plt.legend()

        legend = plt.legend(prop={"size": 6}, loc=9, handlelength=2)
        legend.draw_frame(False)
        plt.ylim((0, 1))
        removeAxes(shift=0)

    plt.show()

#plotDiagonalCorrelation()


def plotCrossValidation():
    "main figure subplot with corss-validation"
    matplotlib.rcParams['font.sans-serif'] = 'Arial'
    plt.figure(figsize=(1, 1))
    FG = HiCdataset(workingFile1, myGenome)
    FG.load(GMFrag)

    Tanay = binnedData(1000000)
    Tanay.simpleLoad("GM-all-10p", "GM-1")
        #need to create these datasets using fragment-level analysis
    Tanay.simpleLoad("GM-all-90p", "GM-9")
    Tanay.removePoorRegions()
    Tanay.iterativeCorrectWithSS()
    Tanay.removeZeros()
    b1, b2 = (Tanay.biasDict["GM-1"], Tanay.biasDict["GM-9"])
    cPickle.dump((b1, b2), open("CrossValidatioN", 'wb'))
    ax = plt.gca()
    b1, b2 = cPickle.load(open("CrossValidatioN", 'rb'))
    print cr(b1, b2)
    plt.scatter(b1, b2, s=.7, color="k", linewidth=0)
    plt.xlabel(r"10% reads", fontsize=8)
    plt.ylabel(r"90% reads", fontsize=8)
    plt.xlim((0, 1.5))
    plt.ylim((0, 1.5))
    plt.xticks([0, 0.5, 1, 1.5])
    plt.yticks([0, 0.5, 1, 1.5])
    removeAxes(shift=0)
    fs = 6
    for xlabel_i in ax.get_xticklabels():
        xlabel_i.set_fontsize(fs)
    for xlabel_i in ax.get_yticklabels():
        xlabel_i.set_fontsize(fs)
    plt.show()


#plotCrossValidation()

def saddlePlot():
    "plot of values ordered by Eig1GW"

    #plt.figure(figsize = (1.5,1.5))
    plt.figure(figsize=(3, 3))
    Tanay = binnedData(1000000)
    Tanay.simpleLoad("../data/GM-all-hg18-1M", "GM-all")
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()
    Tanay.removeZeros()
    Tanay.fakeCis()
    Tanay.iterativeCorrectWithoutSS()
    Tanay.doEig()
    PC = Tanay.EIG["GM-all"][:, 0]
    if PC[0] > 0:
        PC = -PC

    def reorder(data, array=PC):
        inds = numpy.argsort(array)
        ndata = data[inds, :]
        return ndata[:, inds]
    toplot = (coarsegrain(reorder(Tanay.dataDict["GM-all"]), 60))
    toplot /= toplot.mean()
    toplot = numpy.log(toplot)
    sh = toplot.shape
    toplot = toplot.reshape((-1))
    ag = numpy.argmax(toplot)
    toplot[ag] = 0
    toplot[ag] = numpy.max(toplot)
    toplot.shape = sh
    toplot[0, -1] = toplot[0, -2]
    toplot[-1, 0] = toplot[-2, 0]
    plt.imshow(toplot, vmin=toplot.min(), vmax=toplot.max(),
               interpolation="nearest")
    cbar = plt.colorbar(orientation="vertical")
    #labels = ["10","100","1000","10000"]
    #cbar.ax.set_xticklabels(labels)
    cbar.ax.set_xlabel("Log(relative contact probability)", fontsize=6)
    for xlabel_i in cbar.ax.get_xticklabels():
        xlabel_i.set_fontsize(6)
    cbar.set_ticks([-0.5, 0, 0.5, 1])
    removeBorder()
    mirnylib.plotting.niceShow()

#saddlePlot()


def doCartoonPlot():
    "simple eigenvector decomposition cartoon"
    "educational purposes only"

    x = 1. * numpy.arange(30)
    a = numpy.sin(0.5 * x)
    b = 0.07 * (15 - x)

    #mat_img(a[:,None] * a[None,])
    #mat_img(b[:,None] * b[None,:])
    noise = 0.2 * numpy.random.randn(len(x), len(x))
    noise2 = 0.5 * (noise + noise.T)
    signal = a[:, None] * a[None, ] + b[:, None] * b[None, :] + noise2
    PC = EIG(signal)
    s1 = project(signal, PC[:, 0])
    s2 = project(signal, PC[:, 1])

    datas = [signal, s2, s1, s1 + s2]
    datas = [i - 0.5 * (i.min() + i.max()) for i in datas]
    mins = min([i.min() for i in datas])
    maxs = max([i.max() for i in datas])
    mylen = 5
    begs = [0, 0, mylen, mylen]
    begs2 = [1, mylen + 2, 1, mylen + 2]

    for i in xrange(4):
        ax = plt.subplot2grid((2 * mylen + 2, 2 * mylen), (begs2[
            i], begs[i]), rowspan=mylen, colspan=mylen)
        plt.imshow(datas[i], vmin=mins, vmax=maxs,
                   interpolation="nearest")
        removeBorder(ax=ax)

    ax = plt.subplot2grid(
        (2 * mylen + 2, 2 * mylen), (0, mylen), colspan=mylen)
    removeBorder(ax=ax)
    plt.plot(PC[:, 0], 'k')
    ax = plt.subplot2grid(
        (2 * mylen + 2, 2 * mylen), (mylen + 1, 0), colspan=mylen)
    removeBorder(ax=ax)
    plt.plot(PC[:, 1], 'k')

    plt.show()
    plt.imshow(signal - s1 - s2)
    mat_img(s1)
    mat_img(s2)
    mat_img(signal - s1 - s2)


def compareWithGenomicFeatures():
    mirnylib.systemutils.setExceptionHook()
    Tanay = experimentalBinnedData(1000000, myGenome)
    Tanay.simpleLoad(GM1M, "GM-all")
    Tanay.simpleLoad(tcc1M, "tcc")

    datasets = {"CTCF": "wgEncodeBroadChipSeqSignalGm12878Ctcf.wig",
                "H3K27me3": "wgEncodeBroadChipSeqSignalGm12878H3k27me3.wig",
                "H3K4me1": "wgEncodeBroadChipSeqSignalGm12878H3k4me1.wig",
                "H3K4me3": "wgEncodeBroadChipSeqSignalGm12878H3k4me3.wig",
                "H4K20me1": "wgEncodeBroadChipSeqSignalGm12878H4k20me1.wig",
                "H3K27ac": "wgEncodeBroadChipSeqSignalGm12878H3k27ac.wig",
                "H3K36me3": "wgEncodeBroadChipSeqSignalGm12878H3k36me3.wig",
                "H3K4me2": "wgEncodeBroadChipSeqSignalGm12878H3k4me2.wig",
                "H3K9ac": "wgEncodeBroadChipSeqSignalGm12878H3k9ac.wig"
                }
    controlFile = "../../histoneMarks/hg18/wgEncodeBroadChipSeqSignal"\
    "Gm12878Control.wig"

    for key in datasets.keys():
        datasets[key] = os.path.join("../../histoneMarks/hg18", datasets[key])
        datasets[key] = (datasets[key], controlFile)
    datasets["DNAse"] = ("../../histoneMarks/hg18/wgEncodeUwDnaseSeqRawSignal"\
    "Rep1Gm06990.bigWig", None)
    keys = datasets.keys()
    Tanay.loadErezEigenvector1MB("../../ErezPaperData/eigenvectors")
    for key in keys:
        Tanay.loadWigFile(datasets[key][0], key, control=datasets[key][1])
    Tanay.removeChromosome(22)
    Tanay.removeDiagonal()
    Tanay.removePoorRegions()
    Tanay.truncTrans()
    Tanay.fakeCis()
    Tanay.removeZeros()
    Tanay.doEig()

    E1 = Tanay.EigDict["GM-all"][0]
    E2 = Tanay.EigDict["tcc"][0]

    if scipy.stats.spearmanr(Tanay.trackDict["GC"], E1)[0] < 0:
        E1 *= -1
    if scipy.stats.spearmanr(Tanay.trackDict["GC"], E2)[0] < 0:
        E2 *= -1

    print "Eig-GC", scipy.stats.spearmanr(Tanay.trackDict["GC"], E1)
    print "Eigtcc-GC", scipy.stats.spearmanr(Tanay.trackDict["GC"], E2)

    Tanay.trackDict["Eig1GW_Lieberman"] = E1
    Tanay.trackDict["Eig1GW_TCC"] = E2
    Tanay.trackDict["EigLieberman2009"] = Tanay.trackDict["Erez"]
    keys = keys + ["GC", "EigLieberman2009", "Eig1GW_Lieberman", "Eig1GW_TCC"]
    print "Key\t", "key-GC\t", "key-EigLieberman2009\t",
    print "key-EigGW_Lieberman\t", "key-EigGW_TCC\t"
    for key in keys:
        c1 = cr(Tanay.trackDict["GC"], Tanay.trackDict[key])
        c2 = cr(E1, Tanay.trackDict[key])
        c3 = cr(E2, Tanay.trackDict[key])
        c4 = cr(Tanay.trackDict["Erez"], Tanay.trackDict[key])
        print key, "\t%lf\t" % c1[0], "%lf\t" % c4[
            0], "%lf\t" % c2[0], "%lf" % c3[0]
    print

compareWithGenomicFeatures()


def plotTanayGenomicFeature():
    """Shows how genomic feature is spawned by Eig1, not Tanay domains
    paper supplementary  figure"""
    Tanay = experimentalBinnedData(1000000, myGenome)
    Tanay.simpleLoad(GM1M, "GM-all")
    Tanay.loadTanayDomains()

    Tanay.loadWigFile("../../histoneMarks/hg18/wgEncodeUwDnaseSeqRawSignal"\
    "Rep1Gm06990.bigWig", label="feature")
    #Tanay.loadWigFile("../../histoneMarks/hg18/wgEncodeBroadChipSeqSignal"\
    #"Gm12878H3k9ac.wig", label = "feature",
    #control = "../../histoneMarks/hg18/wgEncodeBroadChipSeqSignal"\
    #"Gm12878Control.wig")
    #Tanay.loadWigFile("../../histoneMarks/hg18/wgEncodeBroadChipSeqSignal"\
    #"Gm12878H3k4me3.wig", label = "feature",
    #control = "../../histoneMarks/hg18/wgEncodeBroadChipSeqSignal"\
    #"Gm12878Control.wig")

    Tanay.removeDiagonal()
    Tanay.removePoorRegions()
    Tanay.removeZeros()
    Tanay.fakeCis()

    Tanay.doEig()
    E1 = Tanay.EigDict["GM-all"][0]
    E2 = Tanay.EigDict["GM-all"][1]
    GC = Tanay.trackDict["GC"]

    if scipy.stats.spearmanr(E1, GC)[0] < 0:
        E1 = -E1
    if scipy.stats.spearmanr(E2, GC)[0] < 0:
        E2 = -E2

    TD = Tanay.trackDict["TanayDomains"]
    print scipy.stats.spearmanr(Tanay.trackDict["feature"], E1)

    plt.scatter(Tanay.trackDict["feature"], E1, c=TD, s=4, linewidth=0)
    cm = plt.cm.get_cmap("jet")

    print "Our 2r is", (
        numpy.corrcoef(Tanay.trackDict["feature"], E1)[0, 1]) ** 2
    tset = set(TD)
    tfeature = numpy.zeros_like(TD, dtype=float)
    feature = Tanay.trackDict["feature"]
    for i in tset:
        #print i
        #print numpy.mean(feature[TD==i])
        tfeature[TD == i] = numpy.mean(feature[TD == i])
        #print tfeature
    print "Tanay r2 is", (numpy.corrcoef(tfeature, E1)[0, 1]) ** 2

    plt.legend([matplotlib.lines.Line2D([0], [0], color=cm(i), marker="o",
                                        markersize=8, linewidth=0) \
                for i in [0.333, 0.666, 0.999]],
               ["Active", "Centromere proximal", "Centromere distal"], loc=2)

    plt.ylabel("Eig1GW")
    #plt.xlabel("UWashington DNAse")
    #plt.xlabel("H3K9ac ChIP-seq")
    plt.xlabel("H3K4me3 ChIP-seq")
    plt.title("Color represents domain from (Yaffe 2011)")
    niceShow(subplotAdjust=(0.13, 0.11, 0.97, 0.92))


def plotCorrelationAtDifferentBinning():
    """Plots figure with correlation at different binning.
    Note the caching and creating of binned heatmaps flags below.
    Suppplementary paper figure
    """

    sizes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    setExceptionHook()

    cache = False
    create = False

    if create == True:
        if cache == True:
            #-------------------standard version code-----------------
            FR = fragmentHiC.HiCdataset("bla", "../../data/hg18",
                                        override=False, inMemory=True)
            FR.load("../../ErezPaperData/hg18/GM-HindIII-hg18_refined.frag")

            FR3 = fragmentHiC.HiCdataset("bla", "../../data/hg18",
                                         override=False, inMemory=True)
            FR3.load("../../ErezPaperData/hg18/GM-HindIII-hg18_refined.frag")

            FR2 = fragmentHiC.HiCdataset("bla", "../../data/hg18",
                                         override=False, inMemory=True)
            FR2.load("../../ErezPaperData/hg18/GM-NcoI-hg18_refined.frag")

            #----------------------cross-check code----------------
#            FR = fragmentHiC.HiCdataset("bla", "../data/hg18",
#                                        override=False, inMemory=True)
#            FR.load("../ErezPaperData/hg18/GM-NcoI-hg18_refined.frag")
#
#            FR3 = fragmentHiC.HiCdataset("bla", "../data/hg18",
#                                         override=False, inMemory=True)
#            FR3.load("../ErezPaperData/hg18/GM-NcoI-hg18_refined.frag")
#
#            FR2 = fragmentHiC.HiCdataset("bla", "../data/hg18",
#                                         override=False, inMemory=True)
#            FR2.load("../ErezPaperData/hg18/GM-HindIII-hg18_refined.frag")
            #-------end corss-check code ---------------------------------
            #--------Filter only trans DS reads-----------------
            FR.maskFilter(FR.DS * (FR.chrms1 != FR.chrms2))
            FR2.maskFilter(FR2.DS * (FR2.chrms1 != FR2.chrms2))
            FR3.maskFilter(FR3.DS * (FR3.chrms1 != FR3.chrms2))

            #Now create two halfs of one dataset and down-sample second dataset
            #----------------------standard version code--------
            fraction = 0.5 * len(FR.DS) / float(len(FR2.DS))

            rarray = numpy.random.random(len(FR.DS))
            mask1 = rarray < 0.5
            mask3 = rarray >= 0.5
            mask2 = numpy.random.random(len(FR2.DS)) < fraction

            #-------------------- cross-check code---------
            #fraction = 0.5 * len(FR2.DS) / float(len(FR.DS))

            #rarray = numpy.random.random(len(FR.DS))
            #mask1 =  rarray  < fraction
            #mask3 = (rarray > fraction) * (rarray < fraction * 2)
            #mask2 =  numpy.random.random(len(FR2.DS)) > 0.5
            #-----------------------------------------

            FR.maskFilter(mask1)
            FR2.maskFilter(mask2)
            FR3.maskFilter(mask3)

            FR.save("../../tcc/working/cache1")
            FR2.save("../../tcc/working/cache2")
            FR3.save("../../tcc/working/cache3")
        else:
            FR = fragmentHiC.HiCdataset("bla", "../../data/hg18",
                                        override=False, inMemory=True)
            FR.load("../../tcc/working/cache1")

            FR3 = fragmentHiC.HiCdataset("bla", "../../data/hg18",
                                         override=False, inMemory=True)
            FR3.load("../../tcc/working/cache3")

            FR2 = fragmentHiC.HiCdataset("bla", "../../data/hg18",
                                         override=False, inMemory=True)
            FR2.load("../../tcc/working/cache2")

        for size in sizes:
            FR.saveHeatmap("../../tcc/working/HindIII_%d.hm" %
                           size, size * 1000000)
            FR2.saveHeatmap("../../tcc/working/NcoI_%d.hm" %
                            size, size * 1000000)
            FR3.saveHeatmap("../../tcc/working/control_%d.hm" %
                            size, size * 1000000)

    p1 = []
    p2 = []
    p3 = []
    p4 = []
    evs = []
    for size in sizes:

        BD = binnedDataAnalysis(size * 1000000, "../../data/hg18")
        BD.simpleLoad("../../tcc/working/HindIII_%d.hm" % size, "HindIII")
        BD.simpleLoad("../../tcc/working/NcoI_%d.hm" % size, "NcoI")
        BD.simpleLoad("../../tcc/working/control_%d.hm" % size, "control")
        BD.removeDiagonal()
        BD.removePoorRegions(cutoff=2)
        BD.removeCis()

        data1 = BD.dataDict["HindIII"]
        data2 = BD.dataDict["NcoI"]
        data3 = BD.dataDict["control"]

        mask = (numpy.sum(
            data1, axis=0) > 0) * (numpy.sum(data2, axis=0) > 0)
        validMask = mask[:, None] * mask[None, :]
        transmask = BD.chromosomeIndex[:, None] != BD.chromosomeIndex[None, :]
        cormask = transmask * validMask

        c1 = scipy.stats.spearmanr(data1[cormask], data2[cormask])[0]
        c4 = scipy.stats.spearmanr(data1[cormask], data3[cormask])[0]

        if size == 1:
            evs.append(BD.interchromosomalValues("HindIII"))
            evs.append(BD.interchromosomalValues("NcoI"))
            evs.append(BD.interchromosomalValues("control"))
        p4.append(c4)
        p1.append(c1)

        print "size\t%d\traw:" % size, c1,
        BD.removeZeros()
        BD.fakeCis()  # does iterative correction as well
        BD.restoreZeros(value=0)

        data1 = BD.dataDict["HindIII"]
        data2 = BD.dataDict["NcoI"]
        data3 = BD.dataDict["control"]
        c2 = scipy.stats.spearmanr(data1[cormask], data2[cormask])[0]
        c3 = scipy.stats.spearmanr(data1[cormask], data3[cormask])[0]

        if size == 1:
            evs.append(BD.interchromosomalValues("HindIII"))
            evs.append(BD.interchromosomalValues("NcoI"))
            evs.append(BD.interchromosomalValues("control"))
            print evs

        p3.append(c3)
        p2.append(c2)

        print "\tcorrected:", c2, "\tcontrol", c3

    plt.plot(sizes, p1, label="Raw data, between enzymes")
    plt.plot(sizes, p2, label="Iteratively corrected, between")
    plt.plot(sizes, p3, label="IC, within")
    plt.xlabel("Bin size, MB")
    plt.xticks(range(1, 11))
    plt.ylabel("Spearman correlation coefficient")
    plt.legend()
    niceShow()

    setExceptionHook()
    0 / 0


def plotSixHeatmaps():
    "Plots 6 heatmaps to a correlation supplementary figure"
    plt.figure(figsize=(3, 3))
    for size in [10]:
        BD = binnedDataAnalysis(size * 1000000, "../../data/hg18")
        BD.simpleLoad("../../tcc/working/HindIII_%d.hm" % size, "HindIII")
        BD.simpleLoad("../../tcc/working/NcoI_%d.hm" % size, "NcoI")
        BD.simpleLoad("../../tcc/working/control_%d.hm" % size, "control")
        BD.removeDiagonal()
        BD.removePoorRegions(cutoff=5)
        BD.removeZeros()

        beg1 = BD.chromosomeStarts[0]
        beg2 = BD.chromosomeStarts[1]
        end1 = BD.chromosomeEnds[0] - 1
        end2 = BD.chromosomeEnds[1] - 1

        data1 = BD.dataDict["HindIII"][beg1:end1, beg2:end2]
        data2 = BD.dataDict["NcoI"][beg1:end1, beg2:end2]
        data3 = BD.dataDict["control"][beg1:end1, beg2:end2]

        def minmax(*args):
            mi = min([i.min() for i in args])
            ma = max([i.max() for i in args])
            return mi, ma
        vmin, vmax = minmax(data1, data2)

        plt.subplot(321)
        plt.imshow(data1, interpolation="nearest", vmin=vmin, vmax=vmax)
        plt.title("HindIII, raw")
        plt.colorbar()
        plt.subplot(322)
        plt.imshow(data2, interpolation="nearest", vmin=vmin, vmax=vmax)
        plt.colorbar()
        plt.title("NcoI, raw")

        BD.iterativeCorrectWithoutSS()
        data1 = BD.dataDict["HindIII"][beg1:end1, beg2:end2]
        data2 = BD.dataDict["NcoI"][beg1:end1, beg2:end2]
        data3 = BD.dataDict["control"][beg1:end1, beg2:end2]
        vmin, vmax = minmax(data1, data2, data3)

        plt.subplot(323)
        plt.imshow(data1, interpolation="nearest", vmin=vmin, vmax=vmax)
        plt.colorbar()
        plt.title("HindIII, IC",)
        plt.subplot(324)
        plt.imshow(data2, interpolation="nearest", vmin=vmin, vmax=vmax)
        plt.colorbar()
        plt.title("NcoI, IC")

        plt.subplot(325)
        plt.imshow(data1, interpolation="nearest", vmin=vmin, vmax=vmax)
        plt.colorbar()
        plt.title("HindIII, 50%")
        plt.subplot(326)
        plt.imshow(data3, interpolation="nearest", vmin=vmin, vmax=vmax)
        plt.colorbar()
        plt.title("HindIII, 50%")
        plt.show()


def heatmapFromHotFragments(dataset="../../mouse/data/combined/\
                    mouse1_merged.frag",
                    workingFile="../../tcc/working/workingMouse.frag",
                    cacheFile="../../tcc/working/workingMouseFiltered.frag",
                    genomeFolder="../../data/mm9",
                    label=""):
    mirnylib.systemutils.setExceptionHook()

    if not os.path.exists(cacheFile):
        FH = HiCdataset(workingFile, genomeFolder)
        FH.load(dataset)
        FH.filterRsiteStart(offset=5)
        FH.filterDuplicates()
        #TR.save(filename[1]+".dat")
        FH.filterLarge()
        FH.maskFilter(FH.DS)
        FH.save(cacheFile)

    FH = HiCdataset(workingFile, genomeFolder)
    FH.load(cacheFile)
    fs = FH.fragmentSum()
    p996, p998 = numpy.percentile(fs, [99.6, 99.8])
    frag96 = FH.ufragments[fs > p996]
    frag98 = FH.ufragments[(fs > p998)]

    FH.maskFilter(arrayInArray(
        FH.fragids1, frag96) + arrayInArray(FH.fragids2, frag96))
    hm98100 = FH.buildAllHeatmap(5000000)
    FH.maskFilter(arrayInArray(
        FH.fragids1, frag98) + arrayInArray(FH.fragids2, frag98))
    hm99100 = FH.buildAllHeatmap(5000000)

    plt.subplot(121)

    plt.imshow(numpy.log(hm99100 + 1), interpolation="nearest")
    plt.colorbar()
    plt.title("log # counts for top .2 % of fragments, " + label)

    plt.subplot(122)
    plt.imshow(numpy.log(hm98100 - hm99100 + 1), interpolation="nearest")
    plt.title("log # counts for second top .2%  (.996 - .998), " + label)
    plt.colorbar()
    plt.show()


def plotCisToTransHotFragments(dataset="../../mouse/data/combined/\
                   mouse1_merged.frag",
                   workingFile="../../tcc/working/workingMouse.frag",
                   cacheFile="../../tcc/working/workingMouseFiltered.frag",
                   genomeFolder="../../data/mm9", label=None):
    mirnylib.systemutils.setExceptionHook()
    if not os.path.exists(cacheFile):
        print "caching parsed data"
        FH = HiCdataset(workingFile, genomeFolder)
        FH.load(dataset)
        FH.filterRsiteStart(offset=5)
        FH.filterDuplicates()
        #TR.save(filename[1]+".dat")
        FH.filterLarge()
        FH.maskFilter(FH.DS)
        FH.save(cacheFile)

    FH = HiCdataset(workingFile, genomeFolder)
    FH.load(cacheFile)

    fs = FH.fragmentSum()

    FH.saveFragments()
    FH.maskFilter(FH.chrms1 == FH.chrms2)

    FH.originalFragments()

    fsCis = FH.fragmentSum()
    args = numpy.argsort(fs)

    fsSort = 1. * fs[args]
    fsCisSort = 1. * fsCis[args]

    cisToTrans = fsCisSort / fsSort

    p1, p2, p3 = numpy.percentile(fsSort, [99, 99.5, 99.9])

    bins = mirnylib.numutils.logbins(1, fsSort.max(), 1.08)
    counts = numpy.histogram(fsSort, bins)
    values = numpy.histogram(fsSort, bins, weights=cisToTrans)

    plt.plot(0.5 * (values[1][:-1] + values[1][1:]), values[0] /
             counts[0], '.', label=label)

    for linep in p1, p2, p3:
        plt.vlines(linep, 0, 1)

    plt.xlabel("Counts per fragment")
    plt.ylabel("Cis-to-trans ratio")
    plt.title("Vertical lines are at 99%,99.5% and 99.9% reads per fragment")

    niceShow()


def doAllPlotsForHotFragments():
    plotCisToTransHotFragments(label="Mouse, 2012")  # mouse by default
    heatmapFromHotFragments(label="(Mouse, 2012)")

    plotCisToTransHotFragments(dataset="../../tcc/hg18/\
               tcc-HindIII-hg18_merged.frag",
               workingFile="../../tcc/working/workingTcc.frag",
               cacheFile="../../tcc/working/workingTccFiltered.frag",
               genomeFolder="../../data/hg18",
               label="Kalhor, 2011")

    heatmapFromHotFragments(dataset="../../tcc/hg18/\
                tcc-HindIII-hg18_merged.frag",
                workingFile="../../tcc/working/workingTcc.frag",
                cacheFile="../../tcc/working/workingTccFiltered.frag",
                genomeFolder="../../data/hg18",
                label="(Kalhor, 2011)")

    plotCisToTransHotFragments(dataset="../../ErezPaperData/hg18\
                /GM-HindIII-hg18_merged.frag",
                workingFile="../../tcc/working/workingErez.frag",
                cacheFile="../../tcc/working/workingErezFiltered.frag",
                genomeFolder="../../data/hg18",
                label="Lieberman, 2009")

    heatmapFromHotFragments(dataset="../../ErezPaperData/hg18/\
                GM-HindIII-hg18_merged.frag",
                workingFile="../../tcc/working/workingErez.frag",
                cacheFile="../../tcc/working/workingErezFiltered.frag",
                genomeFolder="../../data/hg18",
                label="(Lieberman, 2009)")


def compareMouseWithGenomicFeatures():
    mirnylib.systemutils.setExceptionHook()
    Tanay = experimentalBinnedData(1000000, "../../data/mm9")
    Tanay.simpleLoad("../../mouse/data/combined/mouse_all-1M.hm", "mouse")
    Tanay.fakeTranslocations([(0, 0, None, 12, 52000000, None),
                              (4, 45000000, None, 12, 0, 30000000),
                              (9, 0, 50000000, 12, 0, 35000000)])

    datasets = {"H3k73me3": ("wgEncodeCaltechHistC2c12Ab2621F"\
    "Cntrl50bPcr1xSigRep1.bigWig",
    "wgEncodeCaltechHistC2c12InputFCntrl50bPcr1xSigRep1.bigWig"),
                "H3k79me2": ("wgEncodeCaltechHistC2c12Ab3594FCntrl"\
                "50bPcr1xSigRep1.bigWig",
                "wgEncodeCaltechHistC2c12InputFCntrl50bPcr1xSigRep1.bigWig"),
                "H3ac": ("wgEncodeCaltechHistC2c12H3ac06599FCntrl"\
                "32bPcr2xSigRep1.bigWig",
                "wgEncodeCaltechHistC2c12InputFCntrl32bPcr2xSigRep1.bigWig"),
                "H3k04me3": ("wgEncodeCaltechHistC2c12H3k04me3FCntrl"\
                "50bPcr1xSigRep1.bigWig",
                "wgEncodeCaltechHistC2c12InputFCntrl50bPcr1xSigRep1.bigWig"),
                "H3k27me3": ("wgEncodeCaltechHistC2c12H3k27me3FCntrl"\
                "32bPcr2xSigRep1.bigWig",
                "wgEncodeCaltechHistC2c12InputFCntrl32bPcr2xSigRep1.bigWig"),
                "H3k36me3": ("wgEncodeCaltechHistC2c12H3k36me3FCntrl"\
                "50bPcr1xSigRep1.bigWig",
                "wgEncodeCaltechHistC2c12InputFCntrl50bPcr1xSigRep1.bigWig"),
                "UW_DNAse": ("wgEncodeUwDnaseBcellcd43nC57bl6MAdult"\
                "8wksSigRep1.bigWig", None)
                }

    for key in datasets.keys():
        mark = os.path.join("../../histoneMarks/mm9", datasets[key][0])
        if datasets[key][1] is not None:
            control = os.path.join("../../histoneMarks/mm9", datasets[key][1])
        else:
            control = None
        datasets[key] = (mark, control)

    keys = datasets.keys()

    for key in keys:
        Tanay.loadWigFile(datasets[key][0], key, control=datasets[key][1])

#    H3 = Tanay.trackDict["H3k73me3"]
#    con = Tanay.trackDict["Control"]
#
#    plt.scatter(H3,con,s=1)
#    plt.show()

    keys = keys + ["GC"]

    Tanay.removeDiagonal()
    Tanay.removePoorRegions(cutoff=2)
    Tanay.truncTrans()
    Tanay.removeZeros()
    Tanay.fakeCis()
    Tanay.doEig()
    Tanay.restoreZeros()
    E1 = Tanay.EigDict["mouse"][0]
    #setExceptionHook()
    #0/0
    cPickle.dump(E1, open("mouseEV.pkl", 'w'))

    """
    #Scatter plot for some data
    GC = Tanay.trackDict["GC"]
    plt.scatter(GC,E1,s=4,label = "All genome")
    beg = Tanay.genome.chrmStartsBinCont[12]
    end = Tanay.genome.chrmEndsBinCont[12]
    plt.scatter(GC[beg:end],E1[beg:end],label = "Chromosome 13",color = "g")
    plt.xlabel("GC content")
    plt.ylabel("Eig1GW")
    niceShow()
    """

    #from hiclib.domainSearch import exponentialDomains
    #E1 = exponentialDomains(Tanay.dataDict["GM-all"])
    cr = scipy.stats.spearmanr

    GCCorr = cr(Tanay.trackDict["GC"], E1)
    if GCCorr[0] < 0:
        E1 = -E1
    print "Eig-GC", cr(Tanay.trackDict["GC"], E1)

    print "key\tkey-GC\t,key-feature"
    for key in keys:
        print "%s %.3lf %.3lf" % (key, cr(Tanay.trackDict["GC"],
                                          Tanay.trackDict[key])[0],
                                  cr(E1, Tanay.trackDict[key])[0])

compareMouseWithGenomicFeatures()


def calculateTanayCorrelation():
    "Calculates correlation between datasets, smoothed in a Tanay way"
    BD = binnedData(1000000, "../../data/hg18")
    BD.simpleLoad("../../ErezPaperData/hg18/GM-HindIII-hg18-1M.hm", "HindIII")
    BD.simpleLoad("../../ErezPaperData/hg18/GM-NcoI-hg18-1M.hm", "NcoI")

    def tanaySmooth(matrix):
        matrix = numpy.array(matrix, dtype="double")
        a = numpy.arange(-9, 10)
        mat = 1 / (1. + numpy.abs(a[:, None]) + numpy.abs(a[None, :]))
        return scipy.ndimage.filters.convolve(input=matrix,
                                              weights=mat,
                                              mode="constant")

    def propagateSmooth(data):
        mask1 = numpy.sum(data, axis=0) > 0
        mask = mask1[:, None] * mask1[None, :]
        ret = numpy.zeros_like(data, dtype=float)
        for i in xrange(BD.genome.chrmCount):
            for j in xrange(BD.genome.chrmCount):
                beg1 = BD.chromosomeStarts[i]
                beg2 = BD.chromosomeStarts[j]
                end1 = BD.chromosomeEnds[i]
                end2 = BD.chromosomeEnds[j]
                mymask = mask[beg1:end1, beg2:end2]
                d = data[beg1:end1, beg2:end2]
                toret = tanaySmooth(d) / tanaySmooth(mymask)
                toret[mymask == 0] = 0
                ret[beg1:end1, beg2:end2] = toret
        return ret

    BD.removePoorRegions(cutoff=2)

    BD.removeCis()

    BD.iterativeCorrectWithoutSS()
    data1 = BD.dataDict["HindIII"]
    data2 = BD.dataDict["NcoI"]

    mask = (numpy.sum(data1, axis=0) > 0) * (numpy.sum(data2, axis=0) > 0)
    validMask = mask[:, None] * mask[None, :]
    transmask = BD.chromosomeIndex[:, None] != BD.chromosomeIndex[None, :]
    cormask = transmask * validMask

    d1 = propagateSmooth(data1)
    d2 = propagateSmooth(data2)
    print scipy.stats.spearmanr(d1[cormask], d2[cormask])
