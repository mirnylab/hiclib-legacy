.. hic-lib documentation master file, created by
   sphinx-quickstart on Mon Mar 26 21:38:46 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
 
Documentation for the Hi-C data analysis library by Leonid Mirny lab
====================================================================

Code is available here: https://bitbucket.org/mirnylab/hiclib

Overview
--------

Here we present a collection of tools to map, filter and analyze Hi-C data. 
The libarary is written in Python, an easy-to-learn human-friendly programming language. 

Installation
------------
This library will work only on Linux/Unix based systems. 

.. note:: Windows support is not currently possible because of the usage of weave.inline function to inline C++ code.
	MacOS was reported to work with weave.inline. It requires gcc and additional compiler directives
	for each call of weave.inline in both mirnylib and hiclib.

To install the library, please run install_linux.py in the hiclib/src directory. 
This will add hiclib to PYTHONPATH (in .bashrc & .bash_profile), what will make 
hiclib library accessible from any folder of the computer. 
You can then import hiclib sub-modules as "import hiclib.binnedData", etc. 
(don't forget to restart bash in order to apply the changes made by install_linux.py).

Requirements
------------

You will need mirnylib library to use this library. 
Mirnylib library is publicly available at https://bitbucket.org/mirnylab/mirnylib. 

Both libraries are library not compatible with python 2.5 and python 3.x, but work fine with python 2.6, 2.7.

However, it might require a newer versions of python packages (specifically numpy and scipy), than those available with python 2.6-equipped linux distributions (e.g. Ubuntu 10.04 LTS, ubuntu 11.04+ works fine!).

Hiclib requires the following python libraries: joblib, h5py, pysam, numpy, scipy, matplotlib, numexpr, biopython, bx-python (preferably from bitbucket repo by james_taylor)

Hiclib requires following non-python binaries to be installed: samtools. 

If you're upgrading numpy/scipy/matplotlib using pip (pip install --upgrade numpy), be sure to delete/replace the original package. 
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

Timewise, fragment-based correction of a dataset is relatively quick, being about an hour for 500 million raw reads.
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
It stores data in a compressed format. For small datasets, it can be initialized to hold data in memory only.

.. warning:: H5dict is not equivalent to dict, as it does not provide direct references to arrays on an HDD. 
	For example ``>>> a = h5dict["myarray"] >>> a[0] = 1`` will not change the array in h5dict. You would need to run ``>>> h5dict["myarray"] = a`` to commit the change to the hard disk.  

H5dict file can be converted to txt format or matlab format using a small utilities
h5dictToTxt.py and h5dictToMat.py. Both utilities are a part of mirnylib repository. 

Pipeline structure
------------------

Two sides of the read are mapped to the genome independently using iterative mapping function. 
Acceptable read format is fastq, but .sra format can be incorporated using a bash_reader parameter 
of iterative_mapping and a fastq-dump tool from SRA toolkit. 
Iterative mapping outputs a combination of sam/bam files (autodetect by filename), which can be then parsed and written to an h5dict containing read data using parse_sams function. 

Fragment-level analysis inputs an h5dict, or any dictionary-like interface, containing at minimum information about read positions (however, strand information is highly recommended). 
It stores all the data in an h5dict, and modifies it as the analysis goes. 

.. warning:: fragmentHiC module modifies the data file when any filter is applied. 

An example script shows how  fragment-level analysis can be used to merge multiple datasets together, filter them and save heatmaps to an h5dict.

Heatmaps and vectors of SS reads are then passed to a binnedData class, that can be used to further filter the Hi-C data,
perform multiple types of iterative correction and eigenvector expansion, and compare results with different genomic tracks.


Example scripts
---------------

Our pipeline also presents a few example scripts that can be used to analyze the data. 
These scripts can be applied to any Hi-C data with minimal modification. 

To perform iterative mapping, use examples/iterative_mapping.py scipt. 

To combine data from multiple sequencing lanes, filter all lanes together (including filtering duplicates), and write out heatmaps, use examples/filteringManyDataests. You can modify genome files location inside the script, and provide mapping output files in the datasets.tsv file. By default, this script saves heatmaps for all replicas of an experiments, then for each experiment combines all replicas together and saves combined heatmaps.


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

