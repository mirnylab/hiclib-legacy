.. hic-lib documentation master file, created by
   sphinx-quickstart on Mon Mar 26 21:38:46 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to hic-lib's documentation!
===================================

Installation
------------
This library will work only on Linux. 
To install it please run install_linux.py in the hiclib/src/hiclib directory. 
You can then import functions as "import hiclib.binnedData". 

Requirements
------------

You will need mirnylab library to use this library. That library is installed in a similar way.  

Python2.6 holds now for every part of the code. However, we advice the use of python 2.7 if possible, as it has numerous bug-fixes relative to 2.6. 

Other python libraries required are: joblib, h5py, pysam, numpy, scipy, matplotlib, numexpr, biopython, bx-python (preferably from bitbucket repo)

Non-python libraries are: samtools, 

Ubuntu 11.04 with default numpy-scipy-matplotlib combination works well. However, for ubuntu 10.04 you'll need to update numpy, scipy and biopython to a newer version.

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


 
Contents:

.. toctree::
   :maxdepth: 3

   fragment 
   binneddata


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

