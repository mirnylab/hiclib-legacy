from hiclib.hicShared import byChrEig
from mirnylib.genome import Genome
import matplotlib.pyplot as plt
from mirnylib.systemutils import setExceptionHook
from mirnylib.plotting import nicePlot, maximumContrastList

setExceptionHook()

gen = Genome('../../../../hg19', readChrms=["#","X"])
mychroms = [0,2,5,13,20]
eigs = byChrEig("../fragmentHiC/test-1M-byChr.hm", gen, chromosomes = mychroms)

for j,chrom in enumerate(mychroms): 
    plt.scatter(eigs[j],gen.GCBin[chrom],  color = maximumContrastList[j], label = "Chr {0}".format(chrom+1))

plt.xlabel("eigenvector")
plt.ylabel("GC content")
nicePlot()

