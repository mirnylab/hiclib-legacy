.. openmm-lib documentation master file, created by
   sphinx-quickstart on Mon Mar 26 21:38:46 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to openmm-lib's documentation!
======================================

Installation
------------
This library will work only on Linux. 
To install it please run install_linux.py in the hiclib/src/hiclib directory. 
You can then import functions as "import hiclib.binnedData". 

Requirements
------------

You will need mirnylab library to use this library. That library is installed in a similar way.  

Python2.6 compatibility should hold for most parts of the pipeline, though things were not extensively tested yet (mapping wasn't tested) 

Other libraries required are: joblib, h5py, pysam, numpy, scipy, matplotlib, numexpr, biopython



 
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

