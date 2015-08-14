import numpy as np 
from hiclib.binnedData import binnedDataAnalysis

import os
import sys


workingGenome = "hg18"
genomeFolder = "../../../data/hg19"
if not os.path.exists(genomeFolder):
    try:
        genomeFolder = sys.argv[1]
    except:
        raise StandardError("Please provide hg19 Genome folder in the code or as a first argument")


a = binnedDataAnalysis(1000000, genomeFolder)
a.simpleLoad("test-1M.hm", "test")
a.removeDiagonal()
a.removePoorRegions()
a.iterativeCorrectWithoutSS()
hm = a.dataDict["test"]
hsum = np.sum(hm, axis=0)
var = np.var(hsum[hsum>0])
print var, np.mean(hsum)
assert var < 0.0001
a.removeCis()
a.fakeCis()
a.removeZeros()
a.iterativeCorrectWithoutSS()
a.doEig()
a.export("test", "testHeatmap")
print "everything worked, but no verification of result was made, because we haven't written it yet..."

