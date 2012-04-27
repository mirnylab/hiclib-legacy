from mirnylab.hic import mapping
from mirnylab.h5dict import h5dict

# A. Map the reads iteratively.
mapping.iterative_mapping(
    bowtie_path='~/bin/bowtie2/bowtie2',
    bowtie_index_path='~/bin/bowtie2/index/hg19',
    fastq_path='~/data/hic/1.fastq',
    out_sam_path='~/data/hic/1.bam',
    min_seq_len=25,
    len_step=5,
    nthreads=8,
    max_reads_per_chunk = 10000000,  #optional, to split reads into smaller groups
    temp_dir = "~/data/tmp",         #optional, keep temporary files here
    bowtie_flags='--very-sensitive --score-min L,-0.6,-0.2')

mapping.iterative_mapping(
    bowtie_path='~/bin/bowtie2/bowtie2',
    bowtie_index_path='~/bin/bowtie2/index/hg19',
    fastq_path='~/data/hic/2.fastq',
    out_sam_path='~/data/hic/2.bam',
    min_seq_len=25,
    len_step=5,
    nthreads=8,
    max_reads_per_chunk = 10000000,   
    temp_dir = "~/data/tmp",          
    bowtie_flags='--very-sensitive --score-min L,-0.6,-0.2')

# B. Parse the mapped sequences into a Python data structure.
lib = h5dict('../hic/hic_lib.hdf5')

mapping.parse_sam(
    sam_basename1='~/data/hic/1.sam',
    sam_basename2='~/data/hic/2.sam',
    out_dict=lib,
    genome_db='../data/hg19')

# C. Assign the ultra-sonic fragments to restriction fragments.
mapping.fill_rsites(
    lib=lib, 
    genome_db='../data/hg19',
    enzyme_name='HindIII')
