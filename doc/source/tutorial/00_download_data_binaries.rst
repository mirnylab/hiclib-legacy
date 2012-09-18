0. Download software and data
=============================


Before we start, we need to make sure we have every tool and dataset that we 
may need. The :doc:`installation </index>` section of the intro describes how to
obtain the source codes of the library and gives the list of the dependencies.

Software
--------

First, we need to obtain the bowtie2 aligning software and the hg19 compiled 
index that it requires:

.. literalinclude:: ../../../examples/tutorial/00_download_data_binaries.sh
   :language: bash
   :lines: 1-18

Besides the aligning software, we also need the 
`NCBI SRA toolkit <http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software>`_.
The `original dataset <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199>`_
was published in the SRA format that has to be converted into FASTQ before 
bowtie2 could read it.

.. literalinclude:: ../../../examples/tutorial/00_download_data_binaries.sh
   :language: bash
   :lines: 20-24

Data
----

In this tutorial we are going to process only one file from the original dataset.
This dataset is the 1st replica of the Hi-C experiment on the GM cell line, the 
restriction enzyme is HindIII. The whole original dataset may be found on 
`GEO <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199>`_.

.. literalinclude:: ../../../examples/tutorial/00_download_data_binaries.sh
   :language: bash
   :lines: 26-31

hg19 genome
-----------

Finally, we will need the FASTA files for 
`the latest version of the human genome <http://hgdownload.cse.ucsc.edu/downloads.html#human>`_
and the gap file containing the genomic coordinates of the centromeres:

.. literalinclude:: ../../../examples/tutorial/00_download_data_binaries.sh
   :language: bash
   :lines: 33-42
