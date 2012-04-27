Mapping
=======

The mapping procedure consists of three simple steps.

A. Iterative mapping.
---------------------

At the first step, the raw sequences of both sides of Hi-C molecules are mapped
separately to the genome database using Bowtie2 aligner.

.. literalinclude:: iterative_mapping.py
   :language: python
   :lines: 1-21

The description of the iterative_mapping function and its arguments
This code will generate a series of files with names '~/data/hic/1.bam.25',
'~/data/hic/1.bam.30', etc with the terminating numbers corresponding to the
truncation length.

B. Parsing of mapped reads.
---------------------------

Then, the mapped sequences of both sides should be parsed and combined.

.. literalinclude:: iterative_mapping.py
   :language: python
   :lines: 24-30

The preferable data type for Hi-C data storage is h5dict included into the 
library. h5dict mimics the interface of the standard Python dict, but stores
the data on the hard drive in highly efficient HDF5 format.

C. Restriction fragments assignment.
------------------------------------

Finally, the ultra-sonic fragments mapped onto the genome are assigned to
the restriction fragments.

.. literalinclude:: iterative_mapping.py
   :language: python
   :lines: 33-36


