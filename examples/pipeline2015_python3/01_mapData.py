# (c) 2015 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Maxim Imakaev (imakaev@mit.edu)

import atexit
import glob
import os
import sys
import gzip

import logging
from hiclib import mapping
from mirnylib import h5dict, genome
import pickle
import numpy as np
logging.basicConfig(level=logging.DEBUG)
from mirnylib.systemutils import setExceptionHook

setExceptionHook()

from defineGenome import getGenome

chunkSize = 10000000

#---------------------------
# This code is not, by any means, an "essential" component of the library
# It is a script which I eventually converged to when I need to map and process a ton of Hi-C data
# It was modified to allow mapping huge datasets, such as (Rao, 2014), with minimal storage/RAM/etc.
# It can be run on a single computer, or on the cluster.
#
# This code is smart enough to be run in parallel on the cluster.
# Just spawn many instances of this code; better with ~5s delays between starting two instances (to avoid collisions with lockfiles).
# It may be very benefitial to assign threads to the # of cores, and spawn two instances on each machine.
# Then if one instance is downloading files, or filtering the .fastq files (both using 1 or several threads), the other will be mapping (using all threads)
#
# The code will automatically detect that some other instance is working on a given ID, and will distribute tasks accordingly
# If a program fails in a brutal way (segfault, etc.), it will leave a lock file in "mapped-hg19(mm9,etc)" folder. Delete stray lockfiles and restart the program.
# If program is ctrl+C'ed or fails normally, it will clear the temporary files and lockfiles.

# To run the code, you must change the parameters below:
# -- provide path to genome files (see mirnylib.genome.Genome for how to do this).
# -- provide path to bowtie and bowtie index
# -- provide a list of SRR numbers you want to map
# -- put fastq-dump in the folder with the code

# To map the .fastq.gz files, do the following:
# 1. Separate your fastq files into two strands
# 2. Put your .fastq.gz files (or symlinks to them) in one folder ("fastq" by default, can be change in fastqFolder variable))
# 3. Be sure your files end with _1.fastq.gz and _2.fastq.gz 
#       Or that you have a unique prefix before .fastq.gz (like side1_data.fastq.gz, not side1_SRR12345.fastq.gz)
#       Then just adjust "sidePrefixes" variable below 

# The code below is moderately complicated.
# For .sra, it automatically downloads .sra files from GEO. Replace the wget call with something else if you have already downloaded them.
# It then splits the file in chunks; using fastq-dump for .sra, and plain Python for .fastq.gz
# It then maps each chunk separately
# Once all chunks are mapped, it marks given dataset as completed


mode = "sra"
#mode = "fastq"

# -------------------- Parameter definitions ------------
inFastqDir = "fastq3"  # for mode="fastq" only

#sidePrefixes = ("side1", "side2")   # a version for naming ....side1.fastq.gz
sidePrefixes = ("_1","_2")  # a prefix preceeding .fastq.gz, which will be used to distinguish side 1 and side 2
# If your files are named "run32167_something_side1_somethingElse.fastq.gz", then "side1_somethingElse" should be the prefix. 

threads = 4
tmpDir = "/tmp"  # this will contain up to 3X the size of the largest input .sra file (256GB is enough for (Rao, 2014), but 128 is not)
# Make sure your system drive (where /tmp usually is) has enough space. If it is a small SSD, it may not.
# Also, there is a lot of IO through the tmpDir. Put it on a local drive, not on a network drive, if you can.

genomeName = "cb10"
genome_db = getGenome(genomeName)

bowtiePath = "../bin/bowtie2/bowtie2"
bowtieIndex = "../bin/bowtie2/index/{0}".format(genomeName)  # change this if your index is named differently from the genome
bowtieFlags = "--very-sensitive"

"IDs from GEO (SRR numbers)"
GEOids = list(range(1665087,1665096))
# Set this for for mapping .sra files
# You can do it like this:
# GEOids = range(1658523,1658540) + [398318, 398920,398921]  #taken from an actual study

seqSkipStart = 0  # skip first 2 bp of the read, if you want
minMapLen = 20  # start mapping at this length
# This will adjust iterative mapping automatically
# Other iterative mapping parameters are in "calculateStep" function defenition

# -------------------- End parameter definitions ------------
if not os.path.exists(bowtiePath): raise

def cleanFile(filename):
    if os.path.exists(filename):
        os.remove(filename)


def cleanDirectory(dirName):
    for i in os.listdir(dirName):
        os.remove(os.path.join(dirName, i))


def calculateStep(length, minlen, approxStep=10, maxSteps=4):
    """returns minimum length and step based on the
    length of sequence and proposed minimum length"""

    actualDif = length - minlen
    if actualDif < approxStep * 0.6:
        return length, 100

    numIter = np.array(np.around(actualDif / float(approxStep)), dtype=int)
    if numIter == 0:
        numIter = 1
    if numIter > maxSteps:
        numIter = maxSteps
    actualStep = actualDif / numIter

    minlen = length - actualStep * numIter

    return minlen, actualStep


if mode == "sra":
    iterList = 2 * GEOids
elif mode == "fastq":
    iterList = 2 * sorted(os.listdir(inFastqDir))


for i in  iterList:
    if mode == "sra":
        sraNum = i
        expName = "SRR{0}".format(i)
        i = expName
        num = int(i[3:])
    if mode == "fastq":
        if not i.endswith(sidePrefixes[0] + ".fastq.gz"):
            if not i.endswith(sidePrefixes[1] + ".fastq.gz"):
                print("ignoring file", i, "does not end with fastq.gz")
            continue
        expName = i.replace(sidePrefixes[0] + ".fastq.gz", "")

    fastqFolder = os.path.join(tmpDir, "{0}-fastq-{1}".format(expName, genomeName))
    samFolder = os.path.join(tmpDir, "{0}-sams-{1}".format(expName, genomeName))

    savePath = "mapped-{0}".format(genomeName)
    saveFolder = os.path.join(savePath, expName)


    for folder in [samFolder, saveFolder, fastqFolder]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    lockName = saveFolder + ".lock"
    completedName = os.path.join(saveFolder, "completed")

    if os.path.exists(completedName) and not os.path.exists(lockName):
        print("skipping", expName)
        continue

    if os.path.exists(lockName):
        print("someone is working on", expName)
        continue
    if os.path.exists(completedName) and os.path.exists(lockName):
        raise

    cleanDirectory(fastqFolder)
    cleanDirectory(samFolder)
    atexit.register(cleanDirectory, fastqFolder)
    atexit.register(cleanDirectory, samFolder)


    lock = open(lockName, "w")
    lock.close()
    atexit.register(cleanFile, lockName)

    if mode == "sra":
        sraName = os.path.join(tmpDir, expName + "_original.sra")
        atexit.register(cleanFile, sraName)
        ret = os.system("wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR{2}/{0}/{0}.sra -O {1}".format(expName, sraName, str(sraNum)[:3]))
        if ret != 0:
            raise RuntimeError("Download failed with return code {0}".format(ret))

        counters = mapping.splitSRA(sraName, os.path.join(fastqFolder, expName + "_chunk{0:04d}_fhtagn_side{1}.fastq.gz"), chunkSize)  # creating unique IDs requires creativity
        os.remove(sraName)
        pickle.dump(counters, open(os.path.join(saveFolder, "read_counts"), 'wb'))
    elif mode == "fastq":

        firstSide = os.path.join(inFastqDir, expName + sidePrefixes[0] + ".fastq.gz")
        secondSide = os.path.join(inFastqDir, expName + sidePrefixes[1] + ".fastq.gz")
        counter = mapping.splitSingleFastq(firstSide, os.path.join(fastqFolder, expName + "_chunk{0:04d}_fhtagn_side1.fastq.gz"), chunkSize)
        mapping.splitSingleFastq(secondSide, os.path.join(fastqFolder, expName + "_chunk{0:04d}_fhtagn_side2.fastq.gz"), chunkSize)


    inFiles = [os.path.join(fastqFolder, i) for i in os.listdir(fastqFolder)]
    for i in inFiles:
        atexit.register(cleanFile, i)

    inFiles1 = sorted([i for i in inFiles if "fhtagn_side1" in i])
    inFiles2 = sorted([i for i in inFiles if "fhtagn_side2" in i])
    assert len(inFiles1) == len(inFiles2)
    outFiles = [os.path.join(saveFolder, "chunk{0:04d}.hdf5".format(i + 1)) for i in range(len(inFiles1))]


    def doOne(inData):
        file1, file2, outfile = inData
        print("Mapping {0} and {1} into {2}".format(*inData))


        for onefile in file1, file2:
            a = gzip.open(onefile, 'r')
            a.readline()
            length = len(a.readline()) - 1
            if length < 10:
                raise ValueError("Length of your sequence is {0}. Something is wrong".format(length))
            minlen, step = calculateStep(length - seqSkipStart, minMapLen)

            mapping.iterative_mapping(
                bowtie_path=bowtiePath,
                bowtie_index_path=bowtieIndex,
                fastq_path=onefile,
                out_sam_path=os.path.join(samFolder, os.path.split(onefile)[1] + ".sam"),
                seq_start=seqSkipStart,
                min_seq_len=minlen,  # for bacteria mimimal mappable length is 15 bp, so I start with something slightly longer
                len_step=step,  # and go with a usualy step
                nthreads=threads,  # on intel corei7 CPUs 4 threads are as fast as
                             # 8, but leave some room for you other applications
                # max_reads_per_chunk = 10000000,  #optional, on low-memory machines
                temp_dir=tmpDir,
                bowtie_flags=bowtieFlags,
                )

        os.remove(file1)
        os.remove(file2)

        # Second step. Parse the mapped sequences into a Python data structure,
        #    assign the ultra-sonic fragments to restriction fragments.
        mapped_reads = h5dict.h5dict(outfile)
        sf1, sf2 = [os.path.join(samFolder, os.path.split(onefile)[1] + ".sam") for onefile in [file1, file2]]
        mapping.parse_sam(sam_basename1=sf1, sam_basename2=sf2,
            out_dict=mapped_reads, genome_db=genome_db, save_seqs=False, maxReads=int(chunkSize*1.6), IDLen=50)
        for i in os.listdir(samFolder):
            if (os.path.split(file1)[1] in i) or (os.path.split(file2)[1] in i):
                print("deleting", i)
                os.remove(os.path.join(samFolder, i))

    if len(inFiles1) == 0:
        raise ValueError("No files supplied")
    list(map(doOne, list(zip(inFiles1, inFiles2, outFiles))))
    a = open(completedName, 'w')
    a.close()

    os.remove(lockName)

