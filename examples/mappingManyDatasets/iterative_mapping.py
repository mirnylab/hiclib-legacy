from mirnylib.hic import mapping
from mirnylib.h5dict import h5dict

# A. Map the reads iteratively.

names = open("datasets.txt").readlines()
names = [i.split() for i in names if len(i) > 5]
assert False not in [len(i) == 4 for i in names]


for filename, length, genome, enzyme in names:
    print "started side 1 of ", filename, genome, enzyme
    base = filename.split(".")[0]
    length = int(length)
    mapping.iterative_mapping(
        bowtie_path='/home/magus/HiC2011/bin/bowtie2-2.0.0-beta5/bowtie2',
        bowtie_index_path='/home/magus/HiC2011/bin/'
        'bowtie2-2.0.0-beta5/index/%s' % genome,
        fastq_path='%s.fastq' % base,
        out_sam_path='tmp/%s_1.bam' % base,
        min_seq_len=26,
        len_step=5,
        nthreads=4,
        temp_dir="tmp",
        seq_start=0,
        seq_end=length,
        bowtie_flags='--fast')
    print "started side 2 of ", filename, genome, enzyme
    mapping.iterative_mapping(
        bowtie_path='/home/magus/HiC2011/bin/bowtie2-2.0.0-beta5/bowtie2',
        bowtie_index_path='/home/magus/HiC2011/bin/'\
        'bowtie2-2.0.0-beta5/index/%s' % genome,
        fastq_path='%s.fastq' % base,
        out_sam_path='tmp/%s_2.bam' % base,
        min_seq_len=26,
        len_step=5,
        nthreads=4,
        temp_dir="tmp",
        seq_start=length,
        seq_end=2 * length,

        bowtie_flags='--fast')

    # B. Parse the mapped sequences into a Python data structure.
    lib = h5dict('data/%s.hdf5' % filename)
    print "parsing sams"
    mapping.parse_sam(
        sam_basename1='tmp/%s_1.bam' % base,
        sam_basename2='tmp/%s_2.bam' % base,
        out_dict=lib,
        genome_db='/home/magus/HiC2011/data/%s' % genome)

    # C. Assign the ultra-sonic fragments to restriction fragments.
    mapping.fill_rsites(
        lib=lib,
        genome_db='/home/magus/HiC2011/data/%s' % genome,
        enzyme_name=enzyme)
    print "Finished"
