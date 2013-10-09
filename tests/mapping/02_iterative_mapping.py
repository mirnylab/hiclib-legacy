import os
import glob
import logging

import numpy as np

import mirnylib.genome
import hiclib.mapping
import mirnylib.h5dict

logging.basicConfig(level=logging.DEBUG)

if not os.path.exists('tmp'):
    os.mkdir('tmp')

bowtie_bin_path = '../../bin/bowtie2/bowtie2'
bowtie_index_path = '../../bin/bowtie2/index/sacCer3'
fasta_path = '../../fasta/sacCer3'

genome_db = mirnylib.genome.Genome(fasta_path)

print "Map the insilico generated Hi-C reads..."
hiclib.mapping.iterative_mapping(
    bowtie_path=bowtie_bin_path,
    bowtie_index_path=bowtie_index_path,
    fastq_path='./tmp/insilico_1.fastq',
    out_sam_path='./tmp/insilico_1.bam',
    min_seq_len=25,
    len_step=5,
    nthreads=4,
    tmp_dir='./tmp',
    bowtie_flags='--very-sensitive')

hiclib.mapping.iterative_mapping(
    bowtie_path=bowtie_bin_path,
    bowtie_index_path=bowtie_index_path,
    fastq_path='./tmp/insilico_2.fastq',
    out_sam_path='./tmp/insilico_2.bam',
    min_seq_len=25,
    len_step=5,
    nthreads=4,
    tmp_dir='./tmp',
    bowtie_flags='--very-sensitive')
print "Done!"

print 'Parse the generated BAMs...'
lib = mirnylib.h5dict.h5dict('./tmp/insilico_mapped_reads.hdf5')
hiclib.mapping.parse_sam(
    sam_basename1='./tmp/insilico_1.bam',
    sam_basename2='./tmp/insilico_2.bam',
    out_dict=lib,
    genome_db=genome_db,
    keep_ids=True,
    enzyme_name='HindIII')
print 'Done!'

