from hiclib.fragmentHiC import HiCdataset
from mirnylib.h5dict import h5dict
import matplotlib.pyplot as plt 
import os
import sys
import numpy as np
from mirnylib.systemutils import setExceptionHook
from mirnylib.numutils import coarsegrain
from mirnylib.genome import Genome
setExceptionHook()

workingGenome = "hg19"

genomeFolder = sys.argv[1]
if not os.path.exists(genomeFolder):
    raise StandardError("Please provide hg19 Genome folder in the code or as a first argument")
mygenome = Genome(genomeFolder, readChrms = ["#","X"])

def source(ID):
    return os.path.join("%s-%s.hdf5" % (ID, workingGenome))  # determines path to the parsed file by ID

N = 2000000

chrms1 = np.random.randint(0,22,N)
chrms2 = chrms1.copy()
mask = np.random.random(N) < 0.5
chrms2[mask] = np.random.randint(0,22,mask.sum())
pos1 = np.array(np.array((0.1 + 0.8 * np.random.random(N)) * mygenome.chrmLens[chrms1]), dtype=int)
offset1 = np.exp(3 + np.random.random(N) * (np.log(1000000) - 3)) * (2 * (np.random.random(N)>0.5) - 1 )
pos2 = np.array(pos1 + offset1, dtype=int)

strands1 = np.random.random(N) > 0.5
strands2 = np.random.random(N) > 0.5



mydict = {"chrms1":chrms1, "chrms2":chrms2,"cuts1":pos1,"cuts2":pos2,"strands1":strands1,"strands2":strands2}

TR = HiCdataset("bla", genome=mygenome, enzymeName="MboI",maximumMoleculeLength=500, inMemory=True)
print "\nTesting loading new data without rsite information    "
TR.parseInputData(dictLike=mydict,
                  enzymeToFillRsites="MboI")
TR.filterLarge(cutlarge=50000, cutsmall=100)

sc = TR.plotScaling()

g = mygenome

print sc
plt.title("Scaling should be 1/x")
plt.plot(*sc)

plt.xscale("log")
plt.yscale("log")
plt.show()


