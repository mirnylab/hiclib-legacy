"""
An example script to perform iterative correction of a raw heatmap.
"""

import os

from mirnylib.genome import Genome
from mirnylib.h5dict import h5dict
from hiclib.binnedData import binnedData


inDatasets = {"HindIII": "../../../ErezPaperData/hg18/GM-HindIII-hg18-1M.hm",
              "NcoI": "../../../ErezPaperData/hg18/GM-NcoI-hg18-1M.hm"}

# Write down filenames of output datasets
outFilenames = {}
for key, value in inDatasets.items():
    outFilenames[key] = value.rsplit(".", 1)[0] + "_ICed.hm"

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

# Remove 1% of regions with low coverage
BD.removePoorRegions(cutoff=1)

# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts)
BD.truncTrans(high=0.0005)

# Actually performe iterative correction
BD.iterativeCorrectWithoutSS()

# Export iteratively corrected datasets
for name, filename in outFilenames.items():
    BD.export(name, filename)

print "Finished!"
