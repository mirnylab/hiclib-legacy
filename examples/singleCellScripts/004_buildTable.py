import glob
import os
import pandas as pd
import numpy as np
from hiclib.binnedData import binnedDataAnalysis
from mirnylib.genome import Genome
import mirnylib.numutils
import pickle
from mirnylib.h5dict import h5dict
from mirnylib.systemutils import  setExceptionHook
setExceptionHook()
sampleDict = pd.read_csv("samples.csv", index_col=0)


def getInterValues(filenames, genome):
    allValues = []
    for i in filenames:
        BD = binnedDataAnalysis(1000000, "/home/magus/HiC2011/data/{0}".format(genome), readChrms=["#", "X"])
        BD.simpleLoad("{1}/{0}-1000k.hm".format(i, genome), "test")
        BD.truncTrans(high=0.0001)
        allValues.append(BD.interchromosomalValues("test"))
    return allValues


def coveragePercent(ind, coarsegrainBy=1):
    try:
        covs =  pickle.load(open("coverages/{0}".format(ind), 'rb'))
    except:
        return np.NAN
    if coarsegrainBy != 1:
        covs = [mirnylib.numutils.coarsegrain(i, coarsegrainBy, extendEdge=True) for i in covs]
    totsum = np.sum([len(i) for i in covs])
    moreZero = np.sum([(i>0).sum() for i in covs])
    coverage = moreZero / totsum
    return coverage


def fractionCis20kb(filename):
    hd = h5dict(filename,'r')
    c1 = hd["chrms1"]
    c2 = hd["chrms2"]
    p1 = hd["cuts1"]
    p2 = hd["cuts2"]
    mask = c1 == c2
    cis = mask.sum()
    more20kb = (np.abs(p1[mask] - p2[mask]) > 20000).sum()
    return more20kb / cis




dfs = []

for genome  in ["mm9", 'hg19']:

    genomeObject = Genome("/home/magus/HiC2011/data/{0}".format(genome), readChrms = ["#", "X"])
    filenames = [i.replace(".1000000.cool", "") for i in os.listdir(genome) if
                 ".1000000.cool" in i and ("DpnII" in i or "combined" in i or "sperm" in i.lower() or "ad" in i.lower() )]


    sampleDict = pd.read_csv("samples.csv", index_col=0)
    #nsnDict = sampleDict["Stage"]

    ### Bug fix to read files only contained within sampleDict ###
    # RUN THIS THE FIRST TIME only
    #fnames = []; 
    #for f in list(sampleDict.index): 
    #    thisFile = [file for file in filenames if str(f) in file]
    #    if len(thisFile)>0:
    #        fnames.append(thisFile[0])
    #filenames = fnames
    #
    ### END bug fix # H.B. 23/10/2016
    """
    ## RUN THE SECOND TIME
    badfiles = [i for i in range(31141,31150+1)]
    badfiles.append(35737)
    badfiles.append(35738)
    [badfiles.append(i) for i in range(41669,41695+1)]
    fnames = []
    for f in filenames: 
        if not any(str(s) in f for s in badfiles):
            fnames.append(f)
    filenames = fnames
    # END SECOND RUN
    """
    indexArray = []
    readArray = []
    for i in filenames:
        curIndex = i.split("_")[0]
        try:
            b = curIndex.split("-")
            int(b[0])
            curIndex = b[0]
        except:
            pass
        indexArray.append(curIndex)
        try:
            numIndex = int(curIndex)
            dodo = True
        except:
            totCounts = 0 
            dodo = False
        if dodo:
        
            readCount = glob.glob("../mapped-*/{0}*/read_counts".format(numIndex))
            counts = [pickle.load(open(i, 'rb')) for i in readCount]
            try: 
                totCounts = sum(sum(counts, []))
                print(totCounts)
            except:
                totCounts = 0
        readArray.append(totCounts)



    df = pd.DataFrame({}, index=indexArray)
    df["raw reads"] = readArray
    allStats = []
    cis20kb = []
    for filename in filenames:
        if os.path.exists(os.path.join("statistics/{0}".format(genome), filename + ".stat")):
            fname = os.path.join("statistics/{0}".format(genome), filename + ".stat")
        else:
            fname = os.path.join("statistics/{0}".format(genome), filename + "_refined.stat")
        try:
            lines = open(fname,'r').readlines()
            lines = [i.strip() for i in lines]
            lines = [i.split(":") for i in lines]
            lines = [i for i in lines if len(i) == 2]
            lines = [[j.strip() for j in i] for i in lines]
            statdict = {i[0]:int(i[1]) for i in lines}
            allStats.append(statdict)
        except:
            print(fname)
            allStats.append({})
        refined = "{0}/{1}_refined.frag".format(genome, filename)
        cisfar = fractionCis20kb(refined)
        cis20kb.append(cisfar)
        #print(cisnew, cis20kb)


        #print(fname)a

    SSRat = []
    total = []
    cis = []
    allraw = []
    records = []
    for fn, stat in zip(filenames, allStats):
        try:
            raw = stat["010_MappedSide1"]
        except:
            raw = np.NAN

        allraw.append(raw)
        tot = stat.get("201_DS+SS","1")
        ss = stat.get("202_SSReadsRemoved","1")
        SSRat.append(float(ss) / float(tot))
        cis.append(int(stat.get("401_cisReads","1")))
        total.append(int(stat.get("400_readsAfterFiltering","1")))
        record = {}

        record["Duplicate-like reads removed"] = stat.get("321_quasiDuplicatesRemoved","")
        record["Removed reads in fragments with 8+ interactions"] = stat.get('360_removedMoreThan8readsPerFragment',"")
        record["Removed reads with start/end <500bp away"] = stat.get("210_sameFragmentReadsRemoved","") + stat.get('220_extraDandlingEndsRemoved',"")
        record["Side 1 mapped"] = stat.get("010_MappedSide1","")
        record["Side 2 mapped"] = stat.get("020_MappedSide2","")
        record["Both sides mapped"] = stat.get("200_totalDSReads","")

        #print(stat.keys())
        #exit()

        records.append(record)

    for i in records[0].keys():
        df[i] = [record[i] for record in records]

    covs40 = []
    covs1M = []

    for i in indexArray:
        covs40.append(coveragePercent(i,4))
        covs1M.append(coveragePercent(i, 100))
    df["cis more than 20kb"] = cis20kb
    df["coverage 40kb"] = covs40
    df["coverage 1M"] = covs1M
    df["filenames"] = filenames
    df["SS read fraction"] = SSRat
    df["raw reads (mapped at side 1)"] = allraw
    df["cis"] = cis
    df["total"] = total
    df["cis-to-total"] = np.array(cis) / (1. * np.array(total))


    #inters = getInterValues(filenames, genome)
    #x = np.percentile(inters, 5, axis=1) / np.median(inters, axis=1)
    #df["InterChrMaxToMedian"] = x

    df = df.sort("cis", ascending=False)


    df["genome"] = genome
    dfs.append(df)

df = pd.concat(dfs)
samples = pd.read_csv("samples.csv", index_col=0)
del samples["genome"]
samples.index = samples.index.astype(str)
print(samples.columns)
df = df.join(samples, how="left")
singleMask = []

for i in df.index.values:
    try:
        int(i)
        singleMask.append(True)
    except:
        singleMask.append(False)
df.set_value(~np.array(singleMask),"type","combined")
df = df.ix[~df["type"].isnull()]
df["isCombined"] = df["type"] == "combined"
df = df.sort_values(["isCombined", "genome", "type","cis"], ascending=[False,False,True,False])
#print(df)
#df["type"].fillna("combined", inplace=True)
print(df.columns)
df.to_csv("statistics.csv")

