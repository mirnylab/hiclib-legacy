import matplotlib.pyplot as plt
import numpy as np

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData

genome_db = genome.Genome('data/hg19', readChrms=['#', 'X'])

# Read resolution from the dataset.
raw_heatmap = h5dict.h5dict('data/heatmap-res-1M.hdf5', mode='r') 
resolution = int(raw_heatmap['resolution'])

# Create a binnedData object, load the data.
BD = binnedData.binnedData(resolution, genome_db)
BD.simpleLoad('data/heatmap-res-1M.hdf5', 'HindIII_GM_1')

# Remove the contacts between loci located within the same bin.
BD.removeDiagonal()

# Remove bins with less than half of a bin sequenced.
BD.removeBySequencedCount(0.5)

# Remove 1% of regions with low coverage.
BD.removePoorRegions(cutoff=1)

# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts).
BD.truncTrans(high=0.0005)

# Perform iterative correction.
BD.iterativeCorrectWithoutSS()

# Save the iteratively corrected heatmap.
BD.export('HindIII_GM_1', 'data/IC-heatmap-res-1M.hdf5')

# Plot the heatmap directly.
plotting.plot_matrix(np.log(BD.dataDict['HindIII_GM_1']))
plt.show()
