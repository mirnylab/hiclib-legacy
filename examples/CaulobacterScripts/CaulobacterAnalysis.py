"""
This is a working code which contains many functions used to analyze Caulobacter paper.
It is not well-organized, but has some approaches to Hi-C analysis.

"""


import os
from mirnylib.systemutils import setExceptionHook
import hiclib.fragmentHiC
import mirnylib.numutils
from hiclib.binnedData import binnedData
import cPickle
from mirnylib.genome import Genome
from mirnylib.numutils import ultracorrect, adaptiveSmoothing, trunc, \
    fillDiagonal
from mirnylib.h5dict import h5dict
from mirnylib.plotting import mat_img, removeAxes, removeBorder
from scipy.stats.stats import spearmanr
from scipy.stats import pearsonr
from scipy.ndimage.filters import gaussian_filter
numutils = mirnylib.numutils
import numpy as np
import matplotlib.pyplot as plt
from hiclib.fragmentHiC import HiCStatistics
from matplotlib.patches import Polygon
import matplotlib

def scaling(heatmap):
    """
    returns P(s) of a 10kb-binned heatmap
    """
    N = len(heatmap)
    inds = range(1, int(0.9 * N))
    values = [np.mean(np.diagonal(heatmap, i)) for i in inds]
    return [10 * i for i in inds], values



datasets = [  # "dnaA_0hr", "dnaA_4hr"
             'Wildtype_0min_BglII_rep1',
              #'Wildtype_0min_BglII_rep2',
             #'0min_cell_cycle',
              #"Wildtype_0min_NcoI",
            #  "empty_plasmid_van_locus",
            # "rsaA_relocated_to_van_locus",
            # "rsaA_relocated_to_xyl_locus",
             #"WT_0min_better_sync_seq_5",

            #, '30min_cell_cycle', '45min_cell_cycle', '60min_cell_cycle', '75min_cell_cycle', '105min_cell_cycle'
            #, "105min_cell_cycle"
             'Deletion_HU',
             'Novobiocin_new',  # "Novobiocin_treated",
                "Rifampicin",
              #"simulated_main_model", "simulated_gyrase",
             # "simulated_rif",
               #, '105min_cell_cycle'
             "Deletion_SMC",
             #  "Rifampicin",
             ]
#uncomment to include/exclude datasets



names = {i:i for i in datasets}
#by default each dataset is names by filename
#but some names are changed through the dictionary below.
names.update({'Wildtype_0min_BglII_rep1' : "Wild type rep1",
              #'Wildtype_0min_BglII_rep2' : "Wild type rep2",
               "Wildtype_0min_NcoI" : "Wild type, NcoI",
              "simulated_realDomains" : "Simulated wild type",
             '0min_cell_cycle' : "0 min cell cycle",
             # "Wildtype_0min_NcoI",
              "empty_plasmid_van_locus" : "Empty plasmid at van locus",
             "rsaA_relocated_to_van_locus" : "rsaA reolcated to van locus",
             "rsaA_relocated_to_xyl_locus" : "rsaA relocated to xyl locus",
             #"WT_0min_better_sync_seq_5",

             '30min_cell_cycle' : "30 min cell cycle",
             '45min_cell_cycle' : "45 min cell cycle",
              '60min_cell_cycle' : "60 min cell cycle",
               '75min_cell_cycle' : "75 min cell cycle",
                '105min_cell_cycle' : "105 min cell cycle",
             'Deletion_HU' : "hup1hup2 deletion",
              'simulated_HU' : "Simulated hup1hup2 deletion",
             'Novobiocin_new' : "Gyrase inhibition",  # "Novobiocin_treated",
              #"simulated_main_model" ,  #"simulated_gyrase",

              "simulated_gyrase" : "Simulated gyrase inhibition",
              "Rifampicin" : "Rifampicin treated",
               "simulated_rif" : "Simulated rifampicin treated",
             "Deletion_SMC" : "SMC deletion",
             })


# extra functions to get quick access to some filenames

def getGenome():
    return Genome("../data/caul", chrmFileTemplate="%s.fa", readChrms=[])


def frag(name):
    return os.path.join("data", name + "_refined.frag")


def hm(name):
    return os.path.join("data", name + "-10k_overlap.hm")


def enzyme(dataset):
    if "NcoI" in dataset:
        return "NcoI"
    return "BglII"


# example of fragment-based P(s) calculations
def plotScalingsForDifferentDatasets():
    setExceptionHook()  # for debug purposes only

    for  dataset in datasets:

        #Creating one inMemory dataset to get scaling
        FR = hiclib.fragmentHiC.HiCStatistics("bla", getGenome(), inMemory=True)

        #Load data. Do not build fragments - maskFilter will do it anyways
        FR.load(frag(dataset), buildFragments=False)

        #Keep read paris from the same strand - see (Imakaev 2012)
        FR.maskFilter((FR.chrms1 == FR.chrms2) * (FR.strands1 == FR.strands2))

        #perform fragment-based IC
        FR.iterativeCorrectionFromMax()

        #These are coordiates of chromosomal arms
        regions = [(0, 200000, 1750000), (0, 2300000, 3800000)]

        #main P(s) code. Do not consider direct neighbors, but kees second
        #neighbors. Caulobacter is not badly affected by inefficient restriction
        #so we set excludeNeighbors to 1. In Eucaryotes it would be 2-4.
        #useWeights uses weights written by fragment-based IC
        scaling = FR.plotScaling(excludeNeighbors=1, enzyme=enzyme(dataset),
                                 normalize=True,
                                 useWeights=True,
                                 regions=regions,
                                 appendReadCount=True, mindist=5000,
                                 label=dataset, linewidth=2)
        #saving P(s)
        cPickle.dump(scaling, open(os.path.join("scalings", os.path.split(dataset)[-1]), 'w'))
    #Plotting the result
    plt.xlim((5000, 2000000))
    plt.legend()
    plt.show()

#plotScalingsForDifferentDatasets()


def saveAllDatasets():
    """
    An example which saves the heatmap in different colormaps.
    It was used to choose the colormap out of the ones we created.
    """
    if not os.path.exists("savedHeatmaps"):
        os.mkdir("savedHeatmaps")
    heatmaps = ["jet", "fall", "blues" "acidblues"]


    for name, heatmap in zip(names, heatmaps)[::-1]:
        for dataset in datasets:
            hm = "data/dumped/%s-10k_overlap.hm_corrected" % dataset
            plt.figure()
            data = np.loadtxt(hm)
            fillDiagonal(data, np.mean(np.diag(data, 2)) * 1.1, 1)
            fillDiagonal(data, np.mean(np.diag(data, 2)) * 1.1 , -1)
            fillDiagonal(data, np.mean(np.diag(data, 2)) * 1.2, 0)
            plt.imshow(data, origin="lower", cmap=heatmap, interpolation="none", vmin=0, vmax=0.035)
            plt.xticks([0, 100, 200, 300, 400], ["0", "1Mb", "2Mb", "3Mb", "4Mb"])
            plt.yticks([0, 100, 200, 300, 400], ["0", "1Mb", "2Mb", "3Mb", "4Mb"])
            plt.colorbar(orientation="vertical", ticks=[0, 0.01, 0.02, 0.03])
            ax = plt.gca()
            for i, line in enumerate(ax.get_xticklines() + ax.get_yticklines()):
                if i % 2 == 1:  # odd indices
                    line.set_visible(False)
            #plt.show()
            plt.savefig("/home/magus/Dropbox/Caulobacter-chromosome/heatmapsAllFromMax/%s_%s.pdf" % (dataset, name))


#a script to extract directionality ratio from heatmap
def directionalityRatio(dataset, size=20):
    heatmap = 1. * h5dict(hm(dataset))["heatmap"]  # extract heatmap

    #filling in the gaps in the heatmap. Not really needed as heatmaps are with overlaps,
    #so they have no gaps
    for _ in range(1):
        zeros = np.sum(heatmap, axis=0) == 0
        zeros = np.nonzero(zeros)[0]
        heatmap[zeros] = heatmap[zeros - 1]
        heatmap[:, zeros] = heatmap[:, zeros - 1]
    #Following regular IC protocol (see 033_....py)
    mirnylib.numutils.fillDiagonal(heatmap, 0, 0)
    mirnylib.numutils.fillDiagonal(heatmap, 0, 1)
    mirnylib.numutils.fillDiagonal(heatmap, 0, -1)
    heatmap = trunc(heatmap, low=0, high=0.0001)
    heatmap = ultracorrect(heatmap)
    diag2value = np.mean(np.diagonal(heatmap, 2))
    mirnylib.numutils.fillDiagonal(heatmap, 1.5 * diag2value, 0)
    mirnylib.numutils.fillDiagonal(heatmap, 1.2 * diag2value, 1)
    mirnylib.numutils.fillDiagonal(heatmap, 1.2 * diag2value, -1)
    heatmap /= np.mean(np.sum(heatmap, axis=0))

    #Put 9 copies of the heatmap in a huge square - Caulobacter is a ring.
    #this is a cheap-and-dirty way to account for that
    tiledHeatmap = np.hstack([heatmap, heatmap, heatmap])
    tiledHeatmap = np.vstack([tiledHeatmap, tiledHeatmap, tiledHeatmap])
    setExceptionHook()  # debug only
    start = len(heatmap)
    end = 2 * len(heatmap)
    ratios = []
    for mon in xrange(start, end):  #going through the central square
        upstream = tiledHeatmap[mon, mon:mon + size].sum()
        downstream = tiledHeatmap[mon - size:mon, mon].sum()
        #print upstream
        #print downstream
        ratios.append(upstream / (upstream + downstream))  #this is upstream/downstream ratio

    return ratios

def diamondScore(dataset, size=10):
    """
    Extract a so-called "diamond score" - inspired by  Suzana Hadjur talks
    see Sevil Sofueva, EMBO 2013 - Supp Figure 11
    (but this is a bit different from Supp Figure 11!!!)
    """
    heatmap = 1. * h5dict(hm(dataset))["heatmap"]
    for _ in range(1):
        zeros = np.sum(heatmap, axis=0) == 0
        zeros = np.nonzero(zeros)[0]
        heatmap[zeros] = heatmap[zeros - 1]
        heatmap[:, zeros] = heatmap[:, zeros - 1]
    mirnylib.numutils.fillDiagonal(heatmap, 0, 0)
    mirnylib.numutils.fillDiagonal(heatmap, 0, 1)
    mirnylib.numutils.fillDiagonal(heatmap, 0, -1)
    heatmap = trunc(heatmap, low=0, high=0.0001)
    heatmap = ultracorrect(heatmap)
    diag2value = np.mean(np.diagonal(heatmap, 2))
    mirnylib.numutils.fillDiagonal(heatmap, 1.5 * diag2value, 0)
    mirnylib.numutils.fillDiagonal(heatmap, 1.2 * diag2value, 1)
    mirnylib.numutils.fillDiagonal(heatmap, 1.2 * diag2value, -1)
    heatmap /= np.mean(np.sum(heatmap, axis=0))
    tiledHeatmap = np.hstack([heatmap, heatmap, heatmap])
    tiledHeatmap = np.vstack([tiledHeatmap, tiledHeatmap, tiledHeatmap])
    setExceptionHook()
    start = len(heatmap)
    end = 2 * len(heatmap)
    ratios = []
    for mon in xrange(start, end):
        diamond = tiledHeatmap[mon:mon + size, mon:mon - size:-1]
        inds = (np.arange(len(diamond))[:, None] + np.arange(len(diamond))[None, :]) < len(diamond)
        ratios.append(diamond[inds].sum())
    return np.array(ratios) - gaussian_filter(ratios, 30)

    return ratios




def allDirectionalityRatios(ratioFunction):
    """
    A simple plot which calculates all directionality ratios, plots them
    and puts lines at 20 top highly expressed genes (Supp figure from our paper)
    This is mostly matplotlib code.
    """
    if not os.path.exists("savedHeatmaps"):
        os.mkdir("savedHeatmaps")
    wildRatio = np.log(ratioFunction("Wildtype_0min_BglII_rep1"))
    for j, dataset in enumerate(datasets):
        ax = plt.subplot(len(datasets), 1, j + 1)
        curRatio = (ratioFunction(dataset))
        plt.title("{1},  r = {0:.2f}, p={2:.2e}".format(pearsonr(curRatio, wildRatio)[0], names[dataset],
                                                      pearsonr(curRatio, wildRatio)[1]), fontsize=10)
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.tick_params(axis='both', which='minor', labelsize=8)
        plt.plot(curRatio)
        plt.ylim((0.25, 0.75))
        plt.xlim((0, len(curRatio)))
        #plt.ylim((0, 1))
        plt.yticks((0.25, 0.5, 0.75))
        geneCoor = [1162773, 3509071, 1180887, 543099, 1953250, 2522439, 3328524, 1503879, 900483, 242693, 3677144, 3931680, 3677704, 3762707, 3480870, 3829656, 1424678, 901855, 1439056, 3678537]
        genePos = [i / 10000. for i in geneCoor]
        #genePos = []
        for lpos in genePos:
            plt.vlines(lpos , -.8, .8, alpha=0.2, linewidth=1, color="black")
        plt.xticks([0, 50, 100, 150, 200, 250, 300, 350, 400], ["" for i in xrange(9)], fontsize=98)
        removeAxes(ax=ax)
        plt.subplots_adjust(0.07, 0.05, 0.94, 0.95, 0.2, 0.5)



    plt.show()
    exit()
#allDirectionalityRatios(directionalityRatio)
#saveAllDatasets()
#exit()


#function which returns subplot configuration based on number of datasets
def asubplots(x):
    if x == 2: return 2, 1
    if x == 3: return 1, 3
    if x == 4: return 2, 2
    if (x == 5) or (x == 6): return 2, 3
    if (x == 7) or (x == 8): return 2, 4
    if x == 9: return 3, 3
    if x == 10: return 2, 5
    return 3, x / 3 + 1

def subplots(x):
    t = asubplots(x)
    return t[::-1]

def showAllDatasets():
    setExceptionHook()

    #plt.figure(figsize=(25, 15))
    fig = plt.figure()

    #size of the figure
    fw = fig.get_figwidth() * fig.get_dpi()
    fh = fig.get_figheight() * fig.get_dpi()

    #get subplot configuration
    sx, sy = subplots(len(datasets))

    for  j, dataset in enumerate(datasets):
        curPlot = plt.subplot(sx, sy, j + 1)
        heatmap = 1. * h5dict(hm(dataset), 'r')["heatmap"]

        #fill in gaps - obsolete, as heatmaps are with overlaps
        for _ in range(1):
            zeros = np.sum(heatmap, axis=0) == 0
            zeros = np.nonzero(zeros)[0]
            heatmap[zeros] = heatmap[zeros - 1]
            heatmap[:, zeros] = heatmap[:, zeros - 1]

        #regular IC protocol
        mirnylib.numutils.fillDiagonal(heatmap, 0, 0)
        mirnylib.numutils.fillDiagonal(heatmap, 0, 1)
        mirnylib.numutils.fillDiagonal(heatmap, 0, -1)
        heatmap = trunc(heatmap, low=0, high=0.0001)
        heatmap = ultracorrect(heatmap)
        diag2value = np.mean(np.diagonal(heatmap, 2))
        mirnylib.numutils.fillDiagonal(heatmap, 1.5 * diag2value, 0)
        mirnylib.numutils.fillDiagonal(heatmap, 1.2 * diag2value, 1)
        mirnylib.numutils.fillDiagonal(heatmap, 1.2 * diag2value, -1)
        newHeatmap = heatmap

        #Top highly expressed genes
        #genePos = [18, 56, 77, 117, 143, 215, 234, 256, 266, 286, 300, 326, 336, 367, 379]
        geneCoor = [1162773, 3509071, 1180887, 543099, 1953250, 2522439, 3328524, 1503879, 900483, 242693, 3677144, 3931680, 3677704, 3762707, 3480870, 3829656, 1424678, 901855, 1439056, 3678537]

        # here we commited to 10kb resolution - change below if you're not
        genePos = [i / 10000. for i in geneCoor]

        genePos = []

        #putting lines at highly expressed genes
        for lpos in genePos:
            plt.hlines(lpos , 0, 500, linewidth=0.7, color="black", alpha=0.2, zorder=1)
            plt.vlines(lpos , 0, 500, linewidth=0.7, color="black", alpha=0.2, zorder=1)
            pass

        #performing adaptive smoothing
        smoothedHeatmap = adaptiveSmoothing(newHeatmap, 20)
        smoothedHeatmap /= np.mean(np.sum(heatmap, axis=0))

        #print dataset, sum([np.diagonal(smoothedHeatmap, i).sum() for i in range(60, 140)])
        #maps = [[smoothedHeatmap, smoothedHeatmap[:30]],
        #         [smoothedHeatmap[:, :30], smoothedHeatmap[:30, :30]]]
        #smoothedHeatmap = np.hstack([np.vstack(i) for i in maps])

        allx = []
        ally = []

        plt.title(dataset, fontsize=10)
        plt.imshow((smoothedHeatmap), interpolation="none", vmax=0.035, cmap="acidblues", zorder=0)
        #plt.imshow((smoothedHeatmap), interpolation="nearest", vmin=0, vmax=np.exp(-4.5), cmap="fall", zorder=0)
        plt.xticks([])
        plt.yticks([])





        plt.subplots_adjust(left=0.05,  # the left side of the subplots of the figure
      right=0.95,  # the right side of the subplots of the figure
      bottom=0.05,  # the bottom of the subplots of the figure
      top=0.95 ,  # the top of the subplots of the figure
      wspace=0.1,  # the amount of width reserved for blank space between subplots
      hspace=0.2)
        #cPickle.dump(scaling, open(dataset.split("/")[-1] + "scaling", 'w'))
        #plt.ylim((400, 200))
        #plt.xlim((0, 200))

        #code below just puts the P(s) over the heatmap
        N = len(smoothedHeatmap)
        pts = np.array([[1, 0], [N, N], [N, 0]])
        p = Polygon(pts, closed=True, facecolor=(0.8, 0.8, 0.8), linewidth=0, alpha=0.7, zorder=2)
        ax = plt.gca()
        ax.add_patch(p)

        Bbox = matplotlib.transforms.Bbox.from_bounds(.55, .55, .35, .42)
        tBbox = matplotlib.transforms.TransformedBbox(Bbox, ax.transAxes).get_points()
        l, b, w, h = tBbox[0, 0] / fw, tBbox[0, 1] / fh, (tBbox[1, 0] - tBbox[0, 0]) / fw, (tBbox[1, 1] - tBbox[0, 1]) / fh
        axins = fig.add_axes([l, b, w, h], axisbg=(0, 0, 0, 0), xscale="log", yscale="log")
        removeAxes(ax=axins)
        for xlabel_i in axins.get_xticklabels(): xlabel_i.set_fontsize(6)
        for xlabel_i in axins.get_yticklabels(): xlabel_i.set_fontsize(6)

        N = len(smoothedHeatmap)
        st = int(0.05 * N)
        end = int(0.45 * N)
        st2 = int(0.55 * N)
        end2 = int(0.95 * N)
        axins.plot(*scaling(0.5 * (smoothedHeatmap[st:end, st:end] + smoothedHeatmap[st2:end2, st2:end2])), color="blue", label="intra-arm")
        if (dataset in ['Wildtype_0min_BglII_rep1', "ML2000_0hr"]):
            myscaling = scaling(0.5 * (smoothedHeatmap[st:end, st:end] + smoothedHeatmap[st2:end2, st2:end2]))
        #axins.plot(*scaling(smoothedHeatmap[st:end, end2:st2:-1]), color="green", label="inter-arm")
        axins.set_xlabel("kb", fontsize=6)
        axins.set_ylabel("Pc", fontsize=6)
        axins.grid()

        if "myscaling" in locals():
            axins.plot(*myscaling, color="grey")

        #axins.set_xticks([])
        #axins.set_yticks([])
        #axins.tick_params(color="red")

        #axins.set_xlabel("Mb")
        #axins.set_ylabel("Pc")
        for i, line in enumerate(axins.get_xticklines() + axins.get_yticklines()):
            if i % 2 == 1:  # odd indices
                line.set_visible(False)

        #if dataset != "Wildtype_0min_BglII_rep1":
        #    data = cPickle.load(open("scalings/{0}".format(dataset)))
        #    axins.plot(*data, color="blue")

        #axins.xscale("log")
        #axins.yscale("log")

        #end strange code





    plt.show()


showAllDatasets()

#plotScalingsForDifferentDatasets()


