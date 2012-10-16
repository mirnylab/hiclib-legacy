0. Download software and data
=============================

Before we start, we need to make sure we have every tool and dataset that we 
may need. The :doc:`installation </index>` section of the intro describes how to
obtain and install the library and gives the list of the dependencies.

In order to obtain the required software and data, complete the following steps:

- Install bowtie2 as described in the :doc:`installation </index>` section;
- Download and index the human genomic sequence as described in the 
  :doc:`installation </index>` section.
- Install the
  `NCBI SRA toolkit <http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software>`_ 
  via the download_sra.sh script in the hiclib folder.
  The `original dataset <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199>`_
  was published in the SRA format and we will need the toolkit to convert it 
  into FASTQ before bowtie2 could read it.
- Execute the download_sample_data.sh script in the hiclib folder.
  In this tutorial we are going to process only one file from the original dataset.
  This dataset is the 1st replica of the Hi-C experiment on the GM cell line, the 
  restriction enzyme is HindIII. The whole original dataset may be found on 
  `GEO <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199>`_.

