"""
An example script to perform iterative correction of a raw heatmap.
"""

import os
from mirnylib.genome import Genome
from mirnylib.h5dict import h5dict
from hiclib.binnedData import binnedData
import numpy as np

inDatasets = {"HindIII": "../../../ErezPaperData/hg18/GM-HindIII-hg18-1M.hm",
              "NcoI": "../../../ErezPaperData/hg18/GM-NcoI-hg18-1M.hm"}

NumEigenvectors = 3  # number of eigenvectors to compute

# Write down filenames of output datasets
outFilenames = {}
for key, value in inDatasets.items():
    outFilenames[key] = value.rsplit(".", 1)[0] + ".eigenvectors"

genomeFolder = "../../../data/hg18"
readChrms = ["#",  # read all numbered chromosomes
             "X"]   # add X chromosome

for inDataset in inDatasets.values():
    if not os.path.exists(inDataset):
        raise IOError("Raw heatmap file does not exist: {}".format(inDataset))

if not os.path.isdir(genomeFolder):
    raise IOError("Genome folder does not exist")

# When you do this, be sure that readChrms used to save heatmap matches
# readChrms that you define here!
genome = Genome(genomeFolder, readChrms=readChrms)

# Read resolution from one of the datasets
sampleDataset = h5dict(inDatasets.values()[0], mode="r")  # random dataset
resolution = int(sampleDataset["resolution"])

# Define the binnedData object, load data
BD = binnedData(resolution, genome, readChrms)
for name, filename in inDatasets.items():
    BD.simpleLoad(filename, name)

BD.removeDiagonal()

# Remove bins with less than half of a bin sequenced
BD.removeBySequencedCount(0.5)

# We'll do iterative correction and Eigenvector expansion on trans data only!
# We want to remove cis, because later we want to remove poor regions in trans
BD.removeCis()

# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts)
# Do this before removing poor regions, because single blowouts may give
# lots of contacts to a region which does not have much contacts otehrwise.
BD.truncTrans(high=0.0005)

# Remove 1% of regions with low coverage
BD.removePoorRegions(cutoff=1)

# Fake cis counts. Data gets iteratively corrected during this process...
BD.fakeCis()

# Remove bins with zero counts for eigenvector analysis
BD.removeZeros()

# Perform eigenvector expansion. Eigenvectors are now in BD.EigDict
# Corresponding eigenvalues are in BD.eigEigenvalueDict
BD.doEig(numPCs=NumEigenvectors)

# BD.EigDict["HindIII"][0] is a first eigenvector of HindIII dataset
# BD.EigDict["NcoI"][2] is the third eigenvector of HindIII dataset
# etc.

# Now we restore regions with zero counts.
# We specify that we fill them in with zeros. By default it goes with NANs.
BD.restoreZeros(value=0)

# Export eigenvectors to files
for name, filename in outFilenames.items():
    eigenvectors = BD.EigDict[name]
    np.savetxt(filename, eigenvectors)

print "Finished!"
