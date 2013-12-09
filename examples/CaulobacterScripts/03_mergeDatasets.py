"""
This scripts is a rip-off of a large mergeDatasets script with certain adjustments. 
Follow comments along the text. 
"""

from hiclib.fragmentHiC import HiCdataset
import os

from mirnylib import genome

genomeDb = genome.Genome('../data/caul', chrmFileTemplate="%s.fa", readChrms=[])

for expName in os.listdir("caul"):

    TR = HiCdataset("bla", genome=genomeDb, inMemory=True,)   # inMemory, as files are probably small (less than hundreds mililon reads)
    TR.parseInputData("caul/" + expName, removeSS=True)       # We discard SS in our pipeline now
    TR.filterRsiteStart(offset=5)                             # We still do this filter to avoid strange "dangling end-like" molecules
    TR.filterDuplicates()
    #TR.save(out_file+".dat")
    TR.filterLarge(cutlarge=300000, cutsmall=100)             #Don't filter any large fragments. This was relevant for eucaryotes with
                                                              #their megabase-long stretches of repetitive or unmappable regions
    #TR.filterExtreme(cutH=0.0025, cutL=0)                    #ALl fragments in Caulobacter seemed to behave normally 
    TR.save("data/" + expName + "_refined.frag")               #Saving filtered dataset

#Below, saving all datasets at different resolutions. 
#Also, using new feature - fragment overlaps - which assins reads to all bins the fragment crosses. 
   
    TR.saveHeatmap("data/" + expName + "-5k_overlap.hm", 5000, useFragmentOverlap=True)
    TR.saveHeatmap("data/" + expName + "-10k_overlap.hm", 10000, useFragmentOverlap=True)
    TR.saveHeatmap("data/" + expName + "-20k_overlap.hm", 20000, useFragmentOverlap=True)
    TR.saveHeatmap("data/" + expName + "-50k.hm", 50000)
    TR.saveHeatmap("data/" + expName + "-20k.hm", 20000)
    TR.saveHeatmap("data/" + expName + "-10k.hm", 10000)

