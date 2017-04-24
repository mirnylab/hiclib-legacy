from mirnylib.h5dict import h5dict
import pandas as pd 
from mirnylib.systemutils import setExceptionHook
import os

"""
The goal of this file is to create a datasets.tsv file. 

The file has the following structure: 
Filename ExperimentName ReplicateName Genome RestrictionEnzyme

such as: 
mapped-hg19/SRR3141592/chunk0001.hdf5 K562 R1 hg19 HindIII
mapped-hg19/SRR3141592/chunk0002.hdf5 K562 R1 hg19 HindIII
mapped-hg19/SRR3141593/chunk0001.hdf5 K562 R1 hg19 HindIII
mapped-hg19/SRR3141593/chunk0002.hdf5 K562 R1 hg19 HindIII
mapped-hg19/SRR3141594/chunk0001.hdf5 K565 R2 hg19 HindIII
etc. 

Sometimes you will need more than 1000 lines in this file, and it may be a pain to create it manually. 
You can really do anything in any language to do this.

For (Rao, 2014) I had to write a parser of the spreadsheets which came with the paper, because 260 sequencing runs is too much to do manually. 
For many papers I use this template. 

This template will parse the file runs.tsv
That must be provided in the following format 

mapped-hg19/SRR3141592 K562 R1 hg19 HindIII
mapped-hg19/SRR3141593 K562 R1 hg19 HindIII
mapped-hg19/SRR3141594 K565 R2 hg19 HindIII
mapped-hg19/SRR3141595 HeLa R1 hg19 HindIII


etc.
"""


#If you want an easy way, replace PART 1 code with the following code 
SRRs = [i.split() for i in open("runs.tsv").readlines() if len(i) > 4]
        
newfile = open("datasets.tsv",'w')
for i in SRRs: 
    folder = i[0]
    print(folder) 
    filenames = [j for j in os.listdir(folder) if ("chunk" in j) and (j.endswith(".hdf5")) ]
    for fname in filenames:
        try: 
            mydict = h5dict(os.path.join(folder,fname),'r')
        except:
            pass 
        if "strands1" not in mydict:
            raise
        if len(mydict.get_dataset("strands1")) < 10000:
            raise
        newfile.write("{0} {1} {2} {3} {4}\n".format(os.path.join(folder,fname),i[1],
                                                    i[2], i[3], i[4]))
        
    
     
                       
