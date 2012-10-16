import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome

logging.basicConfig(level=logging.DEBUG)

if not os.path.exists('../../data/sample/tmp/'):
    os.mkdir('../../data/sample/tmp/')

# A. Map the reads iteratively.
mapping.iterative_mapping(
    bowtie_path='../../bin/bowtie2/bowtie2',
    bowtie_index_path='../../bin/bowtie2/index/hg19',
    fastq_path='../../data/sample/SRR027956.sra',
    out_sam_path='../../data/sample/SRR027056_1.bam',
    min_seq_len=25,
    len_step=5,
    seq_start=0,
    seq_end=75,
    nthreads=4,  # on intel corei7 CPUs 4 threads are as fast as
                 # 8, but leave some room for you other applications
    #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
    temp_dir='../../data/sample/tmp',  # optional, keep temporary files here
    bowtie_flags='--very-sensitive',
    bash_reader='../../bin/sra/bin/fastq-dump -Z')

mapping.iterative_mapping(
    bowtie_path='../../bin/bowtie2/bowtie2',
    bowtie_index_path='../../bin/bowtie2/index/hg19',
    fastq_path='../../data/sample/SRR027956.sra',
    out_sam_path='../../data/sample/SRR027056_2.bam',
    min_seq_len=25,
    len_step=5,
    seq_start=76,
    seq_end=151,
    nthreads=4,  
    #max_reads_per_chunk = 10000000, 
    temp_dir='../../data/sample/tmp',  
    bowtie_flags='--very-sensitive',
    bash_reader='../../bin/sra/bin/fastq-dump -Z')

# B. Parse the mapped sequences into a Python data structure,
#    assign the ultra-sonic fragments to restriction fragments.
mapped_reads = h5dict.h5dict('../../data/sample/mapped_reads.hdf5')
genome_db    = genome.Genome('../../fasta/hg19', readChrms=['#', 'X'])

mapping.parse_sam(
    sam_basename1='../../data/sample/SRR027056_1.bam',
    sam_basename2='../../data/sample/SRR027056_2.bam',
    out_dict=mapped_reads,
    genome_db=genome_db, 
    enzyme_name='HindIII')

