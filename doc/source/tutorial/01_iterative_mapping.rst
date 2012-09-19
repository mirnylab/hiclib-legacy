1. Map reads to the genome
==========================

Align the raw reads with bowtie2
--------------------------------

At the first step, the raw sequences of both ends of Hi-C molecules are mapped
separately to the genome database using the Bowtie2 aligning software. 
The DNA molecules produced in a Hi-C experiment are composed of two 
segments coming from different genome regions. The exact position of the joint
between the segments is not known and may enter the sequenced region and 
introduce a mapping error. This is why we need to trim the sequences to a 
minimal mappable length and increase it for those sequences that cannot be mapped.
The exact details of the mapping algorithm can be found in 
`the original hiclib publication <http://dx.doi.org/10.1038/nmeth.2148>`_.

Usually, the raw sequences from new-generation sequencers are presented in two 
files, for the left and right sides of the DNA molecules, respectively. These 
files should be mapped with separate commands and the results are stored in 
separate files. In our case, the raw sequences are stored in one file, with 
the sequences of the left and right sides of a molecule stored in the first 76
and next 76 symbols, respectively. We also omit the last base pair as it is 
known to have an inferior quality at least in some of the sequencers.

Finally, we supply the full path to fastq-dump binary to the bash_reader
parameter to convert the .sra data file into the .fastq format.

.. literalinclude:: ../../../examples/tutorial/01_iterative_mapping.py
   :language: python
   :lines: 1-38

This code will generate a series of files with names 'data/SRR027956_1.bam.25',
'data/SRR027956_1.bam.30', etc with the terminating numbers corresponding to the
truncation length. For the further details, please consult 
:doc:`the mapping API </mapping>`.

Parse the bowtie2 output
------------------------

Next, the output of bowtie2 have to be parsed, the separate sides of the reads
combined into pairs and the result stored in a dict-like object. All of this
is done with the :func:`hiclib.mapping.parse_sam` function.

As an input for this function, you have to supply a mirnylib.genome.Genome 
object corresponding to the target organism. The Genome class is a Python 
interface to a FASTA database that facilitates calculation and caching of 
auxiliary genomic properties, such as the positions of the restriction sites
or the GC content averaged over a given resolution. The full description of this
class may be found on :doc:`the corresponding API chapter </genome>`. Note, that 
when you create a Genome object, you specify the chromosomes that you want to 
include. In our case, it is all numbered chromosomes ('#') and X. The reads 
mapped to the chromosomes not included into the Genome object are ignored.

Finally, you also have to specify the name of the restriction enzyme used
in the Hi-C experiment. It is used to assign the Hi-C molecules to the
restriction fragments.

.. literalinclude:: ../../../examples/tutorial/01_iterative_mapping.py
   :language: python
   :lines: 39-49

The preferable data type for Hi-C data storage is :class:`mirnylib.h5dict.h5dict`
that is included into the library. h5dict mimics the interface of the standard
Python dict, but stores the data on the hard drive in the highly efficient HDF5 format.

:download:`Download the source code for this chapter. <../../../examples/tutorial/01_iterative_mapping.py>`
