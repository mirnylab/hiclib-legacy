from hiclib import mapping
from mirnylib.h5dict import h5dict
import os
import glob
import mirnylib.genome

if not os.path.exists("tmp"):
    os.mkdir("tmp")


bowtieBinPath = "../../../bin/bowtie2-2.0.0-beta5/bowtie2"
bowtieIndexPath = "../../../bin/bowtie2-2.0.0-beta5/index/hg19"
genomeFastaFiles = "../../../data/hg19"

if not os.path.isfile(bowtieBinPath):
    print "Please provide path to bowtie binary, not to the folder"
    exit()

if not os.path.isdir(genomeFastaFiles):
    print "Please provide path to folder with fasta files"
    exit()

print "Trying to create a Genome from provided directory..",
myGenome = mirnylib.genome.Genome(genomeFastaFiles)
print "Successfull"
del myGenome

print "Starting mapping.."

# A. Map the reads iteratively.
mapping.iterative_mapping(
    bowtie_path=bowtieBinPath,
    bowtie_index_path=bowtieIndexPath,
    fastq_path="test.fastq",
    out_sam_path='1.bam',
        min_seq_len=25,
        len_step=5,
        seq_start=0,
        seq_end=40,
        nthreads=4,
                # on intel corei7 CPUs 4 threads are as fast as
                #8, but leave some room for you other applications
    #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
    temp_dir="tmp", # optional, keep temporary files here
    bowtie_flags='--very-sensitive')

mapping.iterative_mapping(
    bowtie_path=bowtieBinPath,
    bowtie_index_path=bowtieIndexPath,
    fastq_path='test.fastq',
    out_sam_path='2.bam',
    min_seq_len=25,
    len_step=5,
    seq_start=40,
    seq_end=80,
    nthreads=2,
    temp_dir="tmp",
    bowtie_flags='--very-sensitive')

print "Mapping successfull!"

print "Starting read parsing..."

# B. Parse the mapped sequences into a Python data structure.
lib = h5dict('mappedReads.hdf5')

mapping.parse_sam(
    sam_basename1='1.bam',
    sam_basename2='2.bam',
    out_dict=lib,
    genome_db=genomeFastaFiles)

print "Reads parsed successfully"

print "Starting restriction site search.."
# C. Assign the ultrar8-sonic fragments to restriction fragments.
mapping.fill_rsites(
    lib=lib,
    genome_db=genomeFastaFiles,
    enzyme_name='HindIII')

print "Restriction sites filled in!"
print
print "Because mapping is probabilistic, we cannot verify it precisely."
print
print "Loading resulting library..",
mylib = h5dict("mappedReads.hdf5")
print "Loaded!"

print "Getting number of mapped reads.."
side1 = mylib["chrms1"] > 0
side2 = mylib["chrms2"] > 0
print "Loaded!"

print "First side mapped: ", side1.sum()
print "Second side mapped: ", side2.sum()
print "Both sides mapped: ", (side1 * side2).sum()

if abs(side1.sum() - 4080) > 500:
    print "Number of first side mapped reads is too bad (should be 4080)"
    raise
if abs(side2.sum() - 3487) > 500:
    print "Number of second side mapped reads is too bad (should be 3487)"
    raise

print "Numbers of mapped reads seems to be consistent"
print
print "Mapping finished successfully!"

os.remove("mappedReads.hdf5")
for i in os.listdir("tmp"):
    os.remove(os.path.join("tmp", i))
os.rmdir("tmp")
for i in glob.glob("*.bam.*"):
    os.remove(i)
print "Cleaned up!"
