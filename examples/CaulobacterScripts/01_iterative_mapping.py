"""
This scripts takes fastq files from fastq directory, maps them to the genome and 
saves them to .hdf5 files in a directory "caul". 
Please follow comments along the text. 
"""

import glob
import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome

logging.basicConfig(level=logging.DEBUG)


bowtiePath="../bin/bowtie2-2.1.0/bowtie2"
if not os.path.exists(bowtiePath): raise
fastqDir = "fastq"
bowtieIndex = "../bin/bowtie2-2.1.0/index/caul"

for i in sorted(os.listdir(fastqDir)):
    expName = i
    folder = os.path.join(fastqDir, expName)
    file1 = glob.glob(folder+"/*1.fastq")[0]
    file2 = glob.glob(folder+"/*2.fastq")[0]
    if not os.path.exists(file1): raise
    if not os.path.exists(file2): raise
    
# A. Map the reads iteratively.
    mapping.iterative_mapping(
        bowtie_path=bowtiePath,
        bowtie_index_path=bowtieIndex,
        fastq_path=file1,
        out_sam_path='sams/%s_1.bam' % expName,
        min_seq_len=10,   # for bacteria mimimal mappable length is slightly over 10bp, so I start with 10bp 
        len_step=3,       # and go with a smaller step
        seq_start=0,
        seq_end=40,
        nthreads=4,  # on intel corei7 CPUs 4 threads are as fast as
                     # 8, but leave some room for you other applications
        #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
        temp_dir='tmp',  # optional, keep temporary files here
        bowtie_flags='--very-sensitive')

    mapping.iterative_mapping(
        bowtie_path=bowtiePath,
        bowtie_index_path=bowtieIndex,
        fastq_path=file2,
        out_sam_path='sams/%s_2.bam' % expName,
        min_seq_len=10,
        len_step=3,
        seq_start=0,
        seq_end=40,
        nthreads=4,  # on intel corei7 CPUs 4 threads are as fast as
                     # 8, but leave some room for you other applications
        #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
        temp_dir='tmp',  # optional, keep temporary files here
        bowtie_flags='--very-sensitive')

    # B. Parse the mapped sequences into a Python data structure,
    #    assign the ultra-sonic fragments to restriction fragments.
    mapped_reads = h5dict.h5dict('caul/%s' % expName)
    genome_db    = genome.Genome('../data/caul', chrmFileTemplate="%s.fa", readChrms=[])

    mapping.parse_sam(
        sam_basename1='sams/%s_1.bam' % expName,
        sam_basename2='sams/%s_2.bam' % expName,
        out_dict=mapped_reads,
        genome_db=genome_db, 
        enzyme_name='BglII')

