"""
This is a draft of a script which corrects one or multiple Hi-C datasets. 
It uses BinnedData, meaning that you need to have resolution of 200kb or more
"""
from  mirnylib.h5dictUtils.h5dictToTxt import convertFile  # save txt
from hiclib import binnedData

fnames = ["fname1","fname2"]
names = ["dataset 1", "dataset 2"]
exportnames = ["fname1_corrected","fname2_corrected"] 
resolution = 500000
genFolder = "/folder/to/the/genome/files/and/gap.txt/file/according/to/the/mirnylib.genome/class"

#for one file it would be 
fnames = ["myfile.hm"]
resolution = 500000 
names = ["whatever"]
exportnames = ["filename_corrected"]
genFolder = "genomeFolder"


a = binnedData.binnedData(resolution,genFolder)    #folder should be openable by mirnylib.genome

for name,fname,exportname in zip(names,fnames,exportnames):
    a.simpleLoad(fname, name)

a.removeDiagonal()   #we never ever use diagonal

a.removeBySequencedCount()  # new filter: omit all bins with less than 0.5 coverage by sequenced bases (i.e. bases present in the genome)

a.removePoorRegions(cutoff = 0.5, coverage=True)  # remove .5% bins with the lowest number of records (i.e. non-zero entrees in the matrix)
# This filter was updated to remove bins which have zero contacts and one PCR blowout. Those bins would have many reads, but all reads will be with one or few other bins. 

a.removePoorRegions(cutoff = 0.5, coverage=False)  # standart filter. Cutoff reduced to 0.5 from 2. 
a.truncTrans() # remove PCR blowouts from trans data
a.iterativeCorrectWithoutSS()             #do iterative correction 
for name, exportname in names, exportnames: 
    a.export(name,exportname)

    convertFile(exportname, exportname+"_txt") #actually convert the file 
