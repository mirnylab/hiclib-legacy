import matplotlib.pyplot as plt
import functools
import os
from hiclib.fragmentHiC import HiCdataset
import numpy as np
from mirnylib.genome import Genome
import pandas as pd
from hiclib import hicShared
from mirnylib.h5dict import h5dict
from mirnylib.numutils import coarsegrain, observedOverExpected, completeIC
import joblib
import pickle
import cooler
from mirnylib.numutils import  zoomArray
from hiclib.hicShared import getResolution
from hiclib import binnedData
from mirnylib.systemutils import setExceptionHook


# defining genomes used here
mygenM = Genome("../../data/mm9", readChrms=["#", "X"])
mygenH = Genome("../../data/hg19", readChrms=["#", "X"])
genomDict = {"mm9":mygenM, "hg19":mygenH}

mem = joblib.Memory(".")




def getFrags():
    """
    A method that calculates a set of restriction fragments covered in all hg19, or in all mm9 datasets
     Used for the correct calculation of scalings
    """
    mouse = "/home/magus/HiC2011/DrosophilaSingleCell2015/alternativeFiltering/mm9/oocyte_combined_refined.frag"
    human = "/home/magus/HiC2011/DrosophilaSingleCell2015/alternativeFiltering/hg19/K562_combined_refined.frag"
    a = HiCdataset("bla", os.path.join("../../data/hg19"), enzymeName=1000, inMemory=True)
    a.load(human)
    humanFrags = a.fragmentSum() > 0

    a = HiCdataset("bla", os.path.join("../../data/mm9"), enzymeName=1000, inMemory=True)
    a.load(mouse)
    mouseFrags = a.fragmentSum() > 0
    return humanFrags, mouseFrags

getFrags = mem.cache(getFrags)
humanFrags, mouseFrags = getFrags()
fragdict = {"hg19":humanFrags,"mm9":mouseFrags}

def getLoopsDomains(gen):
    """
    A method that analyzes HICCUPS looplists from GEO entree for (Rao, 2014) paper
    """
    mygen = genomDict[gen]
    if gen == "mm9":
        loops = pd.read_csv("../GSE63525_mouse_lymphoblasts_HiCCUPS_looplist.txt", sep="\t")
        loops["chrms"] = [mygen.label2idx[i[3:]] for i in loops["chr1"].values]
        loops["starts"] = (loops["x1"] + loops["x2"]) / 2
        loops["ends"] = (loops["y1"] + loops["y2"]) / 2
        loops = loops[["chrms", "starts", "ends"]]
        loops = loops.sort(["chrms", "starts"])

        domains = pd.read_csv("../GSE63525_mouse_lymphoblasts_Arrowhead_domainlist.txt", sep="\t")
        domains["chrms"] = [mygen.label2idx[i[3:]] for i in domains["chr1"].values]
        domains["starts"] = domains["x1"]
        domains["ends"] = domains["x2"]
        domains = domains[["chrms", "starts", "ends"]]
        domains = domains.sort(["chrms", "starts"])

        CUTOFF = 80000
        domStartID = domains["chrms"].values * 500000000 + domains["starts"].values
        domEndID = domains["chrms"].values * 500000000 + domains["ends"].values
        loopStartID = loops["chrms"].values * 500000000 + loops["starts"].values
        loopEndID = loops["chrms"].values * 500000000 + loops["ends"].values
        mask = []
        for st, end in zip(loopStartID, loopEndID):
            if (np.abs(domStartID - st) + np.abs(domEndID - end)).min() < CUTOFF:
                mask.append(True)
            else:
                mask.append(False)

        domainsLoops = loops.ix[mask]
        return loops, domains, domainsLoops
    elif gen == "hg19":
        loopsH = pd.read_csv("../GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt", sep="\t")
        loopsH["chrms"] = [mygenH.label2idx[i] for i in loopsH["chr1"].values]
        loopsH["starts"] = (loopsH["x1"] + loopsH["x2"]) / 2
        loopsH["ends"] = (loopsH["y1"] + loopsH["y2"]) / 2
        loopsH = loopsH[["chrms", "starts", "ends"]]
        loopsH = loopsH.sort(["chrms", "starts"])
        loops = loopsH

        domains = pd.read_csv("../GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt", sep="\t")
        domains["chrms"] = [mygen.label2idx[i[0:]] for i in domains["chr1"].values]
        domains["starts"] = domains["x1"]
        domains["ends"] = domains["x2"]
        domains = domains[["chrms", "starts", "ends"]]
        domains = domains.sort(["chrms", "starts"])

        CUTOFF = 80000
        domStartID = domains["chrms"].values * 500000000 + domains["starts"].values
        domEndID = domains["chrms"].values * 500000000 + domains["ends"].values
        loopStartID = loops["chrms"].values * 500000000 + loops["starts"].values
        loopEndID = loops["chrms"].values * 500000000 + loops["ends"].values
        mask = []
        for st, end in zip(loopStartID, loopEndID):
            if (np.abs(domStartID - st) + np.abs(domEndID - end)).min() < CUTOFF:
                mask.append(True)
            else:
                mask.append(False)

        domainsLoops = loops.ix[mask]
        return loops, domains, domainsLoops

getLoopsDomains = mem.cache(getLoopsDomains)
loopDict = {i:getLoopsDomains(i) for i in ["hg19","mm9"]}  # dictionary of loops and domains

secondTypes = {"pronuc_male":["pronuc-male", "pronuc-w-o-inh-male"],
               "pronuc_female":["pronuc-female", "pronuc-w-o-inh-female"],
               "SN-Hoechst":["SN-Hoechst"],
               "NSN-Hoechst":["NSN-Hoechst"],
               "K562":["K562"]}  # data types used for plotting

df = pd.read_csv("statistics.csv", index_col=0)  # reading statistics
df = df.sort_values("cis", ascending=False)
df = df[~pd.isnull(df["filenames"])]
combinedDf = df[df["type"] == "combined"]
singleDf = df[df["type"] != "combined"]

# create a dictionary that lists files for each type, and for each combined type
# restricting analyses only to files with >40k reads
# will be used for several plots
mydf = df[df["cis"] > 40000]  # create dataframe with files with >40k reads
typeDict = {}
for i,grouped in pd.groupby(mydf, mydf["type"]):
    if i in ["combined"]:
        continue
    typeDict[i] = list(grouped.index.values)
secondTypeDict = {i:sum([typeDict[k] for k in j],[]) for i,j in secondTypes.items()}


def fakeData(filename, combinedFilename, results):
    """
    This is a method to create reshuffled data

    :param filename: filename of the data to create reshuffled data from
    :param combinedFilename: filename of combined data to sample reshuffled contacts from
    :param results: h5dict or dict-like to put reshuffled data into
    :return: updated h5dict or dict-like
    """
    h1 = h5dict(filename,'r')
    h2 = h5dict(combinedFilename,'r')
    for key in h1.keys():
        if hasattr(h1[key], "shape") and len(h1[key].shape) == 2:
            data1 = h1[key]
            dss = np.sum(data1, axis=1)
            ds = dss.sum()
            del data1
            data2 = h2[key]
            coverage = np.sum(data2, axis=1)
            marginals = 0.5 * dss  / (coverage + 0.00001)
            if (~np.isfinite(marginals)).sum() > 0:
                print("bad marginals")
            marginals[~np.isfinite(marginals)] = 0
            tocreate = data2 * marginals[:,None]
            del data2
            newdata = np.random.poisson(tocreate)
            del tocreate
            result = newdata + newdata.T
            print(key, ds, result.sum())
            del newdata
            results[key] = result
    return results



def controlLoops(df):
    """
    Creates "control" positions for loops

    :param df: loop (or domain) dataframe
    :return: five copies of loop dataframe shifted randomly by 100kb to 1100kb
    """
    dfs = []
    for i in range(20):
        df = df.copy()
        ran = (np.random.random(len(df)) + 0.1) * 1000000
        if i % 2 == 1:
            ran = -ran
        df["starts"] = df["starts"] + ran
        df["ends"] = df["ends"] + ran
        dfs.append(df)
    df = pd.concat(dfs)
    return df


def tadsToDf(tads, resolution):
    chrms = []
    starts = []
    ends = []
    for i,ar in enumerate(tads):
        chrms.append(np.zeros(len(ar)) + i)
        starts.append(ar[:,0])
        ends.append(ar[:,1])
    chrms, starts, ends = [np.concatenate(i) for i in [chrms, starts, ends]]
    starts = resolution * starts
    ends = resolution * ends
    df = pd.DataFrame({"chrms":chrms, "starts":starts, "ends":ends})
    return df

def averageLoops(loopPositions, filename, pad = 8):
    c = cooler.Cooler(filename)
    if "mm9" in filename:
        gen = "mm9"
    else:
        gen = "hg19"
    mygen = genomDict[gen]
    loops, domains, domainsLoops = loopDict[gen]
    resolution = hicShared.getResolution(filename)

    mymaps = []
    ind = os.path.split(filename)[-1].split("_")[0]
    try:
        int(ind[:5])  # workaround for some filename glitches
        ind = ind[:5]
    except:
        pass
    ind = ind.replace("-10k","")
    coverages = pickle.load(open("coverages/{0}".format(ind),'rb'))

    for mychr in range(mygen.chrmCount):
        mymap = np.zeros((2 * pad, 2 * pad), np.float64)
        totcov1 = np.zeros(2 * pad, np.float64)
        totcov2 = np.zeros(2 * pad, np.float64)
        # for mychr in [1, 2, 3]:
        print(mychr)


        #data = myd.get_dataset("{0} {0}".format(mychr))
        data = c.matrix(sparse=True, balance=False).fetch(mygen.idx2label[mychr]).tocsr()
        mycov = coverages[mychr]
        # plt.imshow(data[:2000, :2000])
        # plt.show()
        current = loopPositions.ix[loopPositions["chrms"] == mychr]
        assert len(current) > 0

        for st, end in zip(current["starts"].values, current["ends"].values):

            if abs(st - end) < 10 * resolution:
                continue
            stBin, endBin = int(np.floor(float(st) / resolution)), int(np.floor(float(end) / resolution))
            # print (stBin, endBin)
            if stBin - pad < 0:
                continue
            if endBin + pad > data.shape[0]:
                continue
            # plt.imshow(data[stBin - pad:stBin + pad, endBin - pad:endBin + pad])
            # plt.show()
            mymap = mymap + data[stBin - pad:stBin + pad, endBin - pad:endBin + pad].toarray()
            totcov1 = totcov1 + mycov[stBin-pad:stBin+pad]
            totcov2 = totcov2 + mycov[endBin-pad:endBin+pad]
        covmap = totcov1[:, None] * totcov2[None,:] + 1
        covmap = covmap / covmap.mean()
        assert mymap.shape == (2*pad, 2*pad)
        assert covmap.shape == (2*pad, 2*pad)
        mymaps.append(mymap / covmap)
    for i in mymaps:
        assert i.shape == (2 * pad, 2*pad)
    return mymaps


def averageLoopsWithControl(loopPositions, filename, cg=1, pad = 8 ):
    mymaps =  averageLoops(loopPositions, filename, pad = pad)
    mymaps2 = averageLoops(controlLoops(loopPositions), filename, pad = pad )
    "bla26"
    if cg != 1:
        mymaps = [coarsegrain(1. * i, cg) / (cg ** 2) for i in mymaps]
        mymaps2 = [coarsegrain(1. * i, cg) / (cg ** 2) for i in mymaps2]
    return mymaps, mymaps2



averageLoopsWithControl = mem.cache(averageLoopsWithControl)



def averageDomains( filename, domains=None, M = 30):
    c = cooler.Cooler(filename)
    minSize = 100000
    maxSize = 1000000
    resolution = hicShared.getResolution(filename)
    if "mm9" in filename:
        gen = "mm9"
    else:
        gen = "hg19"
    mygen = genomDict[gen]
    if domains is None:
        loopPositions = loopDict[gen][1]
    else:
        loopPositions=domains


    maps = []
    for mychr in range(mygen.chrmCount):
        futureMap = np.zeros((3 * M, 3 * M), float)
        # for mychr in [1, 2, 3]:
        print(mychr)

        data = c.matrix(sparse=True, balance=False).fetch(mygen.idx2label[mychr]).tocsc()
        current = loopPositions.ix[loopPositions["chrms"] == mychr]
        assert len(current) > 0

        for st, end in zip(current["starts"].values, current["ends"].values):
            if ((end - st) < minSize) or ((end - st) > maxSize):
                continue

            if abs(st - end) < 10 * resolution:
                continue
            stBin, endBin = int(np.rint(float(st) / resolution)), int(np.rint(float(end) / resolution))
            mylen = endBin - stBin + 1
            if stBin - mylen < 0:
                continue
            if endBin + mylen >= data.shape[0]:
                continue
            singleMap = data[stBin - mylen:endBin + mylen + 1, stBin - mylen:endBin + mylen + 1].toarray()
            # singleMap = data[stBin:endBin, stBin:endBin]
            # print singleMap.shape
            assert len(singleMap) % 3 == 0
            futureMap = futureMap + zoomArray(singleMap, futureMap.shape, order = 1)
        maps.append(futureMap)
    return maps

averageDomains = mem.cache(averageDomains)


def getCoverage(filename):
    c = cooler.Cooler(filename)
    if "mm9" in filename:
        gen = "mm9"
    else:
        gen = "hg19"
    mygen = genomDict[gen]

    myd = h5dict(filename, 'r')

    coverages = []
    for mychr in range(mygen.chrmCount):
        data = c.matrix(sparse=True, balance=False).fetch(mygen.idx2label[mychr])
        coverage = np.sum(data, axis=1)
        coverages.append(np.array(coverage)[:,0])
        assert len(coverages[-1]) == data.shape[0]
    return coverages

getCoverage = mem.cache(getCoverage)

def mysum(x):
    return functools.reduce(lambda x,y:x+y, x)


def doScaling(dataset):
    genome = os.path.split(os.path.split(dataset)[0])[-1]
    a = HiCdataset("bla", os.path.join("../../data/", genome), enzymeName=1000, inMemory=True)
    a.load(dataset)
    # a.load("../hadjurCohesin2012/mm9/AST-WT-AdCre-R1-Hi ndIII_refined.frag")
    a.maskFilter((a.chrms1 == a.chrms2) * (a.strands1 == a.strands2))

    sc = {}
    regions = []
    regions2 = []
    for chromosome in range(0, a.genome.chrmCount):
        cur = [(chromosome, 0, a.genome.cntrMids[chromosome]),
                   (chromosome, a.genome.cntrMids[chromosome],
                    a.genome.chrmLens[chromosome])]
        if chromosome % 2 == 0:
            regions += cur
        else:
            regions2 += cur
    frags = fragdict[genome]
    sc1 = a.plotScaling(excludeNeighbors=2,  mindist=6000, regions=regions,  plot=False, fragids1=frags, fragids2=frags)
    sc2 = a.plotScaling(excludeNeighbors=2,  mindist=6000, regions=regions2, plot=False, fragids1=frags, fragids2=frags)
    sc3 = a.plotScaling(excludeNeighbors=2,  mindist=6000, regions=regions + regions2, plot=False, fragids1=frags, fragids2=frags)
    return sc1, sc2, sc3

doScaling = mem.cache(doScaling)


def doEigenvector(filename, genome):
    if filename == "GC":
        gen = Genome("/home/magus/HiC2011/data/" + genome, readChrms=["#","X"])
        gen.setResolution(1000000)
        GC = np.concatenate(gen.GCBin)
        return GC
    resolution = getResolution(filename)
    BD = binnedData.binnedData(resolution, "/home/magus/HiC2011/data/" + genome, ["#","X"])

    BD.simpleLoad(filename, "bla")
    BD.removeDiagonal()

    BD.removeBySequencedCount(0.5)

    BD.removeCis()
    BD.truncTrans(high=0.0005)
    BD.removePoorRegions(cutoff=1)
    BD.fakeCis()
    BD.removeZeros()
    BD.doEig(numPCs=2)
    BD.restoreZeros(value=0)
    return BD.EigDict["bla"][0]

doEigenvector = mem.cache(doEigenvector)


eigmm = doEigenvector("/home/magus/HiC2011/BingRen2013Mouse/mm9/F123_ES_2014_Ren-all-HindIII-200k.hm", "mm9")
eighg = doEigenvector("/home/magus/HiC2011/Erez2014/hg19/K562_inSitu-all-MboI-200k.hm", "hg19")
eigenDict = {"hg19":eighg, "mm9":eigmm}

def doSaddle(filename, eig, gen):
    c = cooler.Cooler(filename)

    gen = Genome("/home/magus/HiC2011/data/" + gen, readChrms=["#", "X"])

    gen.setResolution(getResolution(filename))
    saddles = []
    for chrom in range(gen.chrmCount):
        saddle = np.zeros((5,5), dtype = float)
        st = gen.chrmStartsBinCont[chrom]
        end = gen.chrmEndsBinCont[chrom]
        cur = c.matrix(balance=False).fetch(gen.idx2label[chrom])
        cur = observedOverExpected(cur)
        mask = np.sum(cur , axis=0) > 0
        cur = cur [mask]
        cur = cur [:, mask]
        GC = eig[st:end]
        GC = GC[mask]
        if len(GC) > 5:
            for i in range(5):
                for j in range(5):
                    G1, G2 = np.percentile(GC, [20 * i, 20 * i + 20])
                    mask1 = (GC > G1) * (GC < G2)

                    G1, G2 = np.percentile(GC, [20 * j, 20 * j + 20])
                    mask2 = (GC > G1) * (GC < G2)
                    saddle[i, j] += cur[np.ix_(mask1, mask2)].mean()
        saddles.append(saddle)

    return saddles


doSaddle = mem.cache(doSaddle)


def doSaddleError(filename, eig, gen, correct=False):


    gen = Genome("/home/magus/HiC2011/data/" + gen, readChrms=["#", "X"])
    cur = 0
    data = h5dict(filename,'r')["heatmap"]
    if correct:
        data = completeIC(data)
    gen.setResolution(getResolution(filename))
    if eig == "GC":
        eig = np.concatenate(gen.GCBin)
    saddles = []
    permutted = []
    saddle = np.zeros((5,5), dtype = float)
    for i in range(100):
        permutted.append(np.zeros((5,5), dtype = float))

    for chrom in range(gen.chrmCount):
        st = gen.chrmStartsBinCont[chrom]
        end = gen.chrmEndsBinCont[chrom]
        cur = data[st:end, st:end]
        cur = observedOverExpected(cur)
        mask = np.sum(cur , axis=0) > 0
        cur = cur [mask]
        cur = cur [:, mask]
        GC = eig[st:end]
        GC = GC[mask]
        if len(GC) > 5:
            for i in range(5):
                for j in range(5):
                    G1, G2 = np.percentile(GC, [20 * i, 20 * i + 20])
                    mask1 = (GC > G1) * (GC < G2)
                    G1, G2 = np.percentile(GC, [20 * j, 20 * j + 20])
                    mask2 = (GC > G1) * (GC < G2)
                    addition = cur[np.ix_(mask1, mask2)]
                    addition = np.reshape(addition, (-1))
                    for k in range(100):
                        resampled = np.random.choice(addition, len(addition), replace=True)
                        permutted[k][i,j] += resampled.mean()
                    saddle[i, j] += addition.mean()
    return saddle, permutted


doSaddleError = mem.cache(doSaddleError)

