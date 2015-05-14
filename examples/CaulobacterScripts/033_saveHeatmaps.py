"""
A script which loads saved heatmaps and manuallly performs iterative correction on them. 
It then saves heatmaps with different degree of smoothing and correction. 
"""

from hiclib.fragmentHiC import HiCdataset
import os
import mirnylib
import numpy
np = numpy
from mirnylib import genome
from mirnylib.numutils import adaptiveSmoothing, trunc, ultracorrect
from mirnylib.h5dict import h5dict


genomeDb = genome.Genome('../data/caul', chrmFileTemplate="%s.fa", readChrms=[])

for expName in os.listdir("caul"):  # this directory contains folders with names of experiments. 
                                    # data will be loaded from different folders 
    TR = HiCdataset("bla", genome=genomeDb, inMemory=True,)   # inMemory, as data are small (<1e8 reads) 
    TR.load("data/" + expName + "_refined.frag")              # load filtered data 

    # Now save all heatmaps with different resolutions, etc. 
    TR.saveHeatmap("data/" + expName + "-5k_overlap.hm", 5000, useFragmentOverlap=True)
    TR.saveHeatmap("data/" + expName + "-10k_overlap.hm", 10000, useFragmentOverlap=True)
    TR.saveHeatmap("data/" + expName + "-20k_overlap.hm", 20000, useFragmentOverlap=True)
    TR.saveHeatmap("data/" + expName + "-50k.hm", 50000)
    TR.saveHeatmap("data/" + expName + "-20k.hm", 20000)
    TR.saveHeatmap("data/" + expName + "-10k.hm", 10000)  
    
    # Now iterate over all saved heatmaps... we can afford doing so as heatmaps are small - this scripts runs in minutes 
    for end in ["-5k_overlap.hm", "-10k_overlap.hm", "-20k_overlap.hm",
                "-10k.hm", "-20k.hm", "-50k.hm"]:
        mydict = h5dict("data/" + expName + end)  # manually loading data from h5dict
        origHeatmap = mydict["heatmap"]           # extracting heatmap 

        heatmap = np.array(origHeatmap)           # Converting to np.array (this string probably does nothing) 
        mirnylib.numutils.fillDiagonal(heatmap, 0, 0) # Fill main and second diagonals with zeros 
        mirnylib.numutils.fillDiagonal(heatmap, 0, 1)
        mirnylib.numutils.fillDiagonal(heatmap, 0, -1)
        heatmap = trunc(heatmap, low=0, high=0.0001)   # truncate heatmap at the  top 0.01% of bins (bin pairs). 
                                                       # This is one bin-pair per 10 bins, i.e. negligible. But will help to 
                                                       # get rid of PCR blowouts and other artifacts 
        # Before doing this step, I got sure that sum of reads in each row/column of the heatmap is more than 100 
        """
        If I was working with eucaryotes, or with highly reprtitive genomes, I would do this: 
        mask = np.sum(heatmap, axis=0) < 100 #take all bins with less than 100 reads
        heatmap[mask,:] = 0 
        heatmap[:,mask] = 0   #set these rows and columns to zero
        """
        heatmap = ultracorrect(heatmap)  # performing iterative correction 
        correctedHeatmap = np.array(heatmap)  # just making a copy (np.asarray does not make a copy, np.array does by default) 
        mask = np.sum(correctedHeatmap, axis=0) == 0  # record mask of rows/columns which have a sum of zero (e.g. no fragments).
        diag2value = np.mean(np.diagonal(correctedHeatmap, 2))  # get value from the third diagonal... the first one left 
        mirnylib.numutils.fillDiagonal(correctedHeatmap, 1.5 * diag2value, 0)   # Fill main diagonal with 1.5 times that
        mirnylib.numutils.fillDiagonal(correctedHeatmap, 1.2 * diag2value, 1)  # Fill second diagonals with 1.2 times that 
        mirnylib.numutils.fillDiagonal(correctedHeatmap, 1.2 * diag2value, -1)  # This is done not to confuse adaptiveSmoothing 
        correctedHeatmap[mask] = 0  # some rows and columns had no values. Now they have at the diagonal. 
        correctedHeatmap[:, mask] = 0  # Here we kill them so that adaptiveSmoothing does not get confused 
        # adaptiveSmoothing autodetects rows and columns with no reads and ignores them in a smart way 
        # It may be not as smart with respect to removed diagonals, though it should generally be OK. 
        smoothedHeatmap = adaptiveSmoothing(correctedHeatmap, 10)  # smooth heatmap to 
        mirnylib.numutils.fillDiagonal(smoothedHeatmap, 0, 0)  # Fill diagonals back with zeros 
        mirnylib.numutils.fillDiagonal(smoothedHeatmap, 0, 1)
        mirnylib.numutils.fillDiagonal(smoothedHeatmap, 0, -1)

        smoothedHeatmap /= np.mean(np.sum(smoothedHeatmap, axis=0))  # renormalize all heatmaps to have mean sum over rows/columns of 1 
        heatmap /= np.mean(np.sum(heatmap, axis=0))

        numpy.savetxt("data/dumped/" + expName + end, origHeatmap, fmt='%.4lf')   # saving heatmaps to .txt 
        numpy.savetxt("data/dumped/" + expName + end + "_corrected", heatmap, fmt='%.4lf')
        numpy.savetxt("data/dumped/" + expName + end + "_smoothed", smoothedHeatmap, fmt='%.4lf')





