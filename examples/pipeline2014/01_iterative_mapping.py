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

logging.basicConfig(level=logging.DEBUG)


genomeName="hg19"

bowtiePath = "../bin/bowtie2/bowtie2"
if not os.path.exists(bowtiePath): raise
fastqDir = "fastq"
bowtieIndex = "../bin/bowtie2/index/{0}".format(genomeName)

savePath = "{0}-mapped".format(genomeName)

if not os.path.exists("sams"):
    os.mkdir("sams")
else:
    os.system("rm -rf sams/*")

if not os.path.exists(savePath):
    os.mkdir(savePath)

for i in sorted(os.listdir("fastq")):
    expName = i
    print i 
    file1 = os.path.join(fastqDir, expName)
    if not os.path.exists(file1): raise
    if not os.path.exists(file1): raise
    lengthFile = os.path.join("lengths", expName)
    length = (int(open(lengthFile).readlines()[0]) - 1)/2
    print length
    finalName = '%s/%s' % (savePath,expName)
    if os.path.exists(finalName):
        print "skipping", finalName
        continue

# A. Map the reads iteratively.
    mapping.iterative_mapping(
        bowtie_path=bowtiePath,
        bowtie_index_path=bowtieIndex,
        fastq_path=file1,
        out_sam_path='sams/%s_1.bam' % expName,
        min_seq_len=30,  # for bacteria mimimal mappable length is 15 bp, so I start with something slightly longer
        len_step=10,  # and go with a usualy step
        nthreads=6,  # on intel corei7 CPUs 4 threads are as fast as
                     # 8, but leave some room for you other applications
        #max_reads_per_chunk = 10000000,  #optional, on low-memory machines

        seq_start=0,
        seq_end=length,
        bash_reader = "fastq-dump -Z",
        )

    mapping.iterative_mapping(
        bowtie_path=bowtiePath,
        bowtie_index_path=bowtieIndex,
        fastq_path=file1,
        out_sam_path='sams/%s_2.bam' % expName,
        min_seq_len=30,
        len_step=10,
        nthreads=6,  # on intel corei7 CPUs 4 threads are as fast as
                     # 8, but leave some room for you other applications
        #max_reads_per_chunk = 10000000,  #optional, on low-memory machines

        seq_start=length,
        seq_end = 2 * length,
        bash_reader = "fastq-dump -Z",
        )

    # B. Parse the mapped sequences into a Python data structure,
    #    assign the ultra-sonic fragments to restriction fragments.
    mapped_reads = h5dict.h5dict(finalName)
    genome_db = genome.Genome('../data/{0}'.format(genomeName), readChrms=["#","X"])

    mapping.parse_sam(
        sam_basename1='sams/%s_1.bam' % expName,
        sam_basename2='sams/%s_2.bam' % expName,
        out_dict=mapped_reads,
        genome_db=genome_db,
        enzyme_name='HindIII',
	save_seqs=False)

