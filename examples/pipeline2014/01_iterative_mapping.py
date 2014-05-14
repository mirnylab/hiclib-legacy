"""
This scripts takes fastq files from fastq directory, maps them to the genome and
saves them to .hdf5 files in a directory "genomeName-mapped".
Please follow comments along the text.
"""

import glob
import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome

import numpy as np
logging.basicConfig(level=logging.DEBUG)


genomeName = "mm10"
threads = 10
bowtiePath = "../bin/bowtie2/bowtie2"
if not os.path.exists(bowtiePath): raise
fastqDir = "fastq"
bowtieIndex = "../bin/bowtie2/index/{0}".format(genomeName)
tmpDir = "/tmp"
samFolder = "sams-{0}".format(genomeName)

savePath = "mapped-{0}".format(genomeName)

if not os.path.exists(samFolder):
    os.mkdir(samFolder)
else:
    os.system("rm -rf {0}/*".format(samFolder))

if not os.path.exists(savePath):
    os.mkdir(savePath)

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



for i in sorted(os.listdir("fastq")):
    expName = i
    print i
    file1 = os.path.join(fastqDir, expName)
    if not os.path.exists(file1): raise
    if not os.path.exists(file1): raise
    lengthFile = os.path.join("lengths", expName)
    length = (int(open(lengthFile).readlines()[0]) - 1) / 2
    print length
    minlen, step = calculateStep(length, 25)


    finalName = '%s/%s.hdf5' % (savePath, expName.replace(".sra", ""))
    print finalName
    if os.path.exists(finalName):
        print "skipping", finalName
        continue

# A. Map the reads iteratively.
    mapping.iterative_mapping(
        bowtie_path=bowtiePath,
        bowtie_index_path=bowtieIndex,
        fastq_path=file1,
        out_sam_path='{0}/{1}_1.bam'.format(samFolder, expName),
        min_seq_len=minlen,  # for bacteria mimimal mappable length is 15 bp, so I start with something slightly longer
        len_step=step,  # and go with a usualy step
        nthreads=threads,  # on intel corei7 CPUs 4 threads are as fast as
                     # 8, but leave some room for you other applications
        #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
        temp_dir=tmpDir,
        seq_start=0,
        seq_end=length,
        bash_reader="fastq-dump -Z",
        bowtie_flags=" --very-sensitive ",
        )

    mapping.iterative_mapping(
        bowtie_path=bowtiePath,
        bowtie_index_path=bowtieIndex,
        fastq_path=file1,
        out_sam_path='{0}/{1}_2.bam'.format(samFolder, expName),
        min_seq_len=minlen,
        len_step=step,
        nthreads=threads,  # on intel corei7 CPUs 4 threads are as fast as
                     # 8, but leave some room for you other applications
        #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
        temp_dir=tmpDir,
        seq_start=length,
        seq_end=2 * length,
        bash_reader="fastq-dump -Z",
        bowtie_flags=" --very-sensitive ",
        )

    # B. Parse the mapped sequences into a Python data structure,
    #    assign the ultra-sonic fragments to restriction fragments.
    mapped_reads = h5dict.h5dict(finalName)
    genome_db = genome.Genome('../data/{0}'.format(genomeName), readChrms=["#", "X"])

    mapping.parse_sam(
        sam_basename1='{0}/{1}_1.bam'.format(samFolder,expName),
        sam_basename2='{0}/{1}_2.bam'.format(samFolder,expName),
        out_dict=mapped_reads,
        genome_db=genome_db,
	save_seqs=False)

