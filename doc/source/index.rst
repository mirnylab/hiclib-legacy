.. hic-lib documentation master file, created by
   sphinx-quickstart on Mon Mar 26 21:38:46 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
 
Welcome to hic-lib's documentation!
===================================

Installation
------------
This library will work only on Linux. 
To install it please run install_linux.py in the hiclib/src directory. 
This will make hiclib module accessible from any folder on the machine 
(don't forget to restart bash in order to apply the changes made by install_linux.py). 
You can then import hiclib sub-modules as "import hiclib.binnedData", etc. 

Requirements
------------

You will need mirnylab library to use this library. That library is installed in a similar way.  

This library is not python 2.5-compatible, but works fine with python 2.6. 
However, it might require a newer versions of python packages, than those available with python 2.6-equipped linux distributions. 
It does not work with python 3.x 

Hiclib requires the following python libraries: joblib, h5py, pysam, numpy, scipy, matplotlib, numexpr, biopython, bx-python (preferably from bitbucket repo)

Hiclib requires following non-python binaries to be installed: samtools. 

Ubuntu 11.04 with default numpy-scipy-matplotlib combination works well. However, for ubuntu 10.04 you'll need to update numpy, scipy and biopython to newer versions.

If you're updating numpy/scipy/matplotlib using pip (pip install --upgrade numpy), be sure to delete/replace the original package. 
Often pip installs new package to a different location, and python still loads the old copy as it looks there first.
You might need to specifically delete files corresponding to the original package, as running "apt-get remove python-numpy" might take down a significant number of packages. Locations of files can be determined from the python console by >>>print numpy.__file__, and version by >>>print numpy.__version__. 

         
Hardwre requirements
--------------------

Fragment-based analysis uses HDD to store all the information about the reads. 
However, at each point in time a few individual tracks are loaded to RAM to perform certain analysis. 
Memory requirements for this part can be estimated as 30-50 bytes per read, dependent on the number of filters used. 
For example, a 400 mln read Hi-C dataset will require no more than 20GB RAM. 

Binned data analysis uses RAM to store all the heatmaps. Application of certain filters will create a copy of a heatmap, one at a time, even if multiple datasets are loaded. 
For example, working with 3 datasets at 200-kb resolution will require: 
(3 GB/200kb)^2 * (8 bytes per float64) * (3 datasets + 1 extra) = 8 GB RAM

Timewise, fragment-based correction of a dataset is relatively quick, being about 1 hour for 500 million raw reads. 
Binned data analysis performance depends on resolution only, usually being of the order of seconds for 1M-resolution analysis, and scaling as (1/resolution^2), 
reaching few minutes for 200kb resolution, and tens minutes for 100kb. 

Overall, time-limiting step is mapping, and it's the only part of the code that is well-paralellized. 
Multi-threading is only utilized by bowtie when mapping reads, and takes about a few days to map 500-million-read dataset on 
a high-level desktop, or a day on a high-performance workstation. 

Storage concept
---------------

mirnylab.h5dict is a default storage used by all parts of the library. 
H5dict is a persistent python dictionary with all the data stored on the HDD. 
It uses HDF5 to store and access the data, and can be viewed using external HDF5 viewers or python. 
It stores data in a compressed format, and can be possibly modified to hold data in RAM only. 

H5dict can be converted to txt format using a small utility called h5dictToTxt.py

Pipeline structure
------------------

Two sides of the read are mapped to the genome independently using iterative mapping function. 
Acceptable read format is fastq, but .sra format can be incorporated using a bash_reader parameter 
of iterative_mapping and a fastq-dump tool from SRA toolkit. 
Iterative mapping outputs a combination of sam/bam files (autodetect by filename), which can be then parsed and written to an h5dict containing read data.

Fragment-level analysis inputs an h5dict, or any dictionary-like interface, containing at minimum information about read positions (directions are highly recommended). 
It stores all the data in an h5dict, and modifies it as the analysis goes. 
An example script shows how  fragment-level analysis can be used to merge multiple datasets together, filter them and save heatmaps to an h5dict. 

Heatmaps and vectors of SS reads are then passed to a binnedData class, that can be used to further filter the Hi-C data,
perform multiple types of iterative correction and eigenvector expansion, and compare results with different genomic tracks.


Example scripts
---------------



Contents
--------

.. toctree::
   :maxdepth: 3
   

   glossary   
   genome
   mapping_tutorial
   fragment 
   binneddata


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

