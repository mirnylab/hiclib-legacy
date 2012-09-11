import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mirnylib.numutils import correct, ultracorrect
from mirnylib.plotting import removeAxes
from mirnylib.genome import Genome
from mirnylib import plotting
import cPickle
from scipy import stats

from hiclib.fragmentHiC import HiCdataset

def plotFigure2c():
    TR = HiCdataset()
    TR.load("GM-all.refined")
    hm = TR.buildHeatmap(1, 1, 1000000, False, False)
    TR.calculateWeights()
    TR.weights = np.ones(len(TR.weights), float)  # if you want to correct just by fragment density, not by length dependence
    hm2 = TR.buildHeatmap(1, 1, 1000000, False, weights=True)
    hm2[np.isnan(hm2)] = 0
    mask = np.sum(hm, axis=0) > 0
    """p1-6 are 6 lines to be plotted, below is plotting only"""
    p1 = np.sum(hm, axis=0)[mask]
    p3 = np.sum(correct(hm), axis=0)[mask]
    p5 = np.sum(ultracorrect(hm, 40), axis=0)[mask]
    p4 = np.sum(correct(hm2), axis=0)[mask]
    p2 = np.sum(hm2, axis=0)[mask]
    p6 = np.sum(ultracorrect(hm2, 40), axis=0)[mask]
    matplotlib.rcParams['font.sans-serif'] = 'Arial'
    dashstyle = (3, 3)
    plt.figure(figsize=(4, 4))

    ax = plt.subplot(2, 1, 1)
    plt.xlim((0, 80))
    plt.ylim((0, 2))
    plt.ylabel("Total coverage", fontsize=8)

    line21 = plt.plot(p1 / p1.mean(), "-", linewidth=1, color="#e5a826")[0]
    line22 = plt.plot(
        p3 / p3.mean(), "--", linewidth=1, color="#e5a826")[0]
    line22.set_dashes(dashstyle)
    line23 = plt.plot(p5 / p5.mean(), linewidth=1, color="grey")[0]

    for xlabel_i in ax.get_xticklabels():
        xlabel_i.set_fontsize(8)
    for xlabel_i in ax.get_yticklabels():
        xlabel_i.set_fontsize(8)
    legend = plt.legend([line21, line22, line23],
                        ["Raw data", "Single correction", "Iterative correction"], prop={"size": 6}, loc=1, handlelength=2)
    legend.draw_frame(False)
    removeAxes(shift=0, ax=ax)

    for i in ax.spines.values():
        i.set_color('none')
    ax.axhline(linewidth=1, color='black')
    ax.axvline(linewidth=1, color='black')

    ax2 = plt.subplot(2, 1, 2, sharex=ax)
    plt.xlim((0, 80))
    plt.ylim((0, 2))
    plt.xlabel("Position on chom 1 (MB)", fontsize=8)
    plt.ylabel("Total coverage", fontsize=8)

    line1 = plt.plot(p4 / p4.mean(), "--", color="#9b3811", linewidth=1)[0]
    line1.set_dashes(dashstyle)
    line2 = plt.plot(p2 / p2.mean(), "-", color="#9b3811", linewidth=1)[0]
    line3 = plt.plot(p6 / p6.mean(), linewidth=1, color="grey")[0]

    for xlabel_i in ax2.get_xticklabels():
        xlabel_i.set_fontsize(8)
    for xlabel_i in ax2.get_yticklabels():
        xlabel_i.set_fontsize(8)

    legend = plt.legend([line2, line1, line3],
                        ["HindIII corrected", "Single correction", "Iterative correction"], prop={"size": 6}, loc=1, handlelength=2)
    legend.draw_frame(False)
    removeAxes(shift=0, ax=ax2)
    plotting.niceShow()

#plotFigure2c()


def doSupplementaryCoveragePlot():
    TR = HiCdataset()
    TR.load("GM-all.refined")
    s1 = TR.fragmentSum(strands=1)
    TR.saveFragments()
    TR.maskFilter(TR.dists1 > 500)
    TR.originalFragments()
    s2 = TR.fragmentSum(strands=1)
    resolution = 1000000

    def coverage(s1, s2, TR):
        genome = Genome()
        genome.createMapping(resolution)
        label = genome.chromosomeStarts[TR.ufragments / TR.fragIDmult -
            1] + (TR.ufragments % TR.fragIDmult) / resolution
        counts = np.bincount(label, weights=s1)
        counts2 = np.bincount(label, weights=s2)
        data = cPickle.load(open("GC1M", 'rb'))
        eigenvector = np.zeros(genome.chromosomeEnds[-1], float)
        inds = np.argsort(counts)
        mask = inds[int(0.02 * len(inds)):]
        for chrom in range(1, 24):
            eigenvector[genome.chromosomeStarts[chrom - 1]:genome.chromosomeStarts[chrom - 1] + len(data[chrom - 1])] = data[chrom - 1]
        eigenvector[eigenvector < 35] = 35
        plt.scatter(counts[mask], counts2[mask], c=eigenvector[
            mask], s=6, linewidth=0)
        print stats.spearmanr(counts[mask], counts2[mask])
        plt.xlabel("Coverage from all reads")
        plt.xticks([0, 5000, 10000, 15000])
        plt.ylabel("Coverage from RBs")
        b = plt.colorbar()
        b.ax.set_xlabel("GC content")
    plt.subplot(121)
    plt.title("HinIII")
    coverage(s1, s2, TR)

    TR = HiCdataset()
    TR.load("GM-NcoI.refined")
    s1 = TR.fragmentSum(strands=1)
    TR.saveFragments()
    TR.maskFilter(TR.dists1 > 500)
    TR.originalFragments()
    s2 = TR.fragmentSum(strands=1)
    resolution = 1000000
    plt.subplot(122)
    plt.title("NcoI")
    coverage(s1, s2, TR)
    plt.show()

