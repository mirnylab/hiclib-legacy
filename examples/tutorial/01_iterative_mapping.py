from hiclib import mapping
from mirnylib.h5dict import h5dict

# A. Map the reads iteratively.
mapping.iterative_mapping(
    bowtie_path='bin/bowtie2-2.0.0-beta7/bowtie2',
    bowtie_index_path='bin/bowtie2-2.0.0-beta7/index/hg19',
    fastq_path='data/SRR027956.sra',
    out_sam_path='SRR027056_1.bam',
    min_seq_len=25,
    len_step=5,
    seq_start=0,
    seq_end=75,
    nthreads=4,  # on intel corei7 CPUs 4 threads are as fast as
                 # 8, but leave some room for you other applications
    #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
    temp_dir='data/tmp',  # optional, keep temporary files here
    bowtie_flags='--very-sensitive',
    bash_reader='bin/sratoolkit.2.1.16-ubuntu32/bin/fastq-dump')

mapping.iterative_mapping(
    bowtie_path='bin/bowtie2-2.0.0-beta7/bowtie2',
    bowtie_index_path='bin/bowtie2-2.0.0-beta7/index/hg19',
    fastq_path='data/SRR027956.sra',
    out_sam_path='SRR027056_2.bam',
    min_seq_len=25,
    len_step=5,
    seq_start=76,
    seq_end=151,
    nthreads=4,  
    #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
    temp_dir='data/tmp',  # optional, keep temporary files here
    bowtie_flags='--very-sensitive',
    bash_reader='bin/sratoolkit.2.1.16-ubuntu32/bin/fastq-dump')

# B. Parse the mapped sequences into a Python data structure.
mapped_reads = h5dict('data/mapped_reads.hdf5')

mapping.parse_sam(
    sam_basename1='data/SRR027956_1.bam',
    sam_basename2='data/SRR027956_2.bam',
    out_dict=mapped_reads,
    genome_db='../data/hg19')

# C. Assign the ultra-sonic fragments to restriction fragments.
mapping.fill_rsites(
    lib=mapped_reads,
    genome_db='data/hg19',
    enzyme_name='HindIII')


