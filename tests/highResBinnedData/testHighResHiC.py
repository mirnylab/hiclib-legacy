"""
This is an example of how HiResBinnedData behaves exactly like binnedData

It uses genome with 5 chromosomes because hiResBinnedData has a small per-chromosome
overhead, and therefore having 23*11 chromosome pairs would slow it down a bit.

Note the trick of converting a global heatmap map to a smaller map.

It tests only iterative correction now.
"""

from hiclib.highResBinnedData import HiResHiC
from mirnylib.genome import Genome
from hiclib.binnedData import binnedData
from mirnylib.h5dict import h5dict
import numpy as np
import os

genome = Genome("../../../data/hg19", readChrms=["1", "2", "3", "4", "5"])

a = HiResHiC(genome, 1000000, "hiResDict", mode='w')
a.loadData(dictLike="../fragmentHiC/test-1M-byChr.hm")
a.removeDiagonal()
a.removePoorRegions(2)
a.iterativeCorrection(1e-10)

b = binnedData(1000000, genome)

data = {"heatmap": h5dict("../fragmentHiC/test-1M.hm")["heatmap"]}
lim = b.genome.chrmEndsBinCont[-1]
data["heatmap"] = data["heatmap"][:lim, :lim]

b.simpleLoad(data, "data")
b.removeDiagonal()
b.removePoorRegions(cutoff=2)
b.iterativeCorrectWithoutSS(tolerance=1e-10)
a.export("testExport")

def compareData():
    dataHigh = a.getCombinedMatrix()
    dataLow = b.dataDict["data"]

    dataHigh /= dataHigh.mean()
    dataLow /= dataLow.mean()

    assert np.abs(dataHigh - dataLow).max() < 1e-6

compareData()



print "Tests finished successfully"

os.remove("hiResDict")
os.remove("testExport")


