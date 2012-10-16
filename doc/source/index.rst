.. raw:: html

    <script type="text/javascript">
    
    var _gaq = _gaq || [];
    
    _gaq.push(['_setAccount', 'UA-33583600-1']);
    
    _gaq.push(['_trackPageview']);
    
    (function() {
    
    var ga = document.createElement('script'); ga.type =
    
    'text/javascript'; ga.async = true;
    
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' :
    
    'http://www') + '.google-analytics.com/ga.js';
    
    var s = document.getElementsByTagName('script')[0];
    
    s.parentNode.insertBefore(ga, s);
    
    })();
    
    </script>

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
The library is written in Python, an easy-to-learn human-friendly programming language. 

This page contains the general information about the library.
Detailed documentation for each module can be found in API documentation below. 

Requirements
------------

The hiclib was designed and tested only under Ubuntu environment, but
may theoretically be ported to the other Linux distributions and OS. 

Hiclib requires the following python libraries: 

- cython (0.16+)
- numpy(1.6+) 
- matplotlib
- scipy
- biopython
- joblib (0.6.3+)
- h5py
- pysam
- bx-python (preferably from 
  `the bitbucket repo <https://bitbucket.org/james_taylor/bx-python/wiki/Home>`_ 
  by james_taylor)
- `mirnylib <https://bitbucket.org/mirnylab/mirnylib>`_ (covered
  in the `installation`_ section below).

hiclib and mirnylib are not compatible with python 2.5 and python 3.x, but work
fine with python 2.6, 2.7.

Besides that, you will need the following binaries and libraries (the names 
are given for the corresponding Ubuntu packages):

- samtools
- hdf5-tools
- libhdf5-serial-dev

Finally, you will also need the bowtie2 mapping software which can be downloaded
manually or installed via the download_bowtie.sh script in the hiclib folder.

Installation
------------

To install the library in Linux, do the following procedure:

1. Install all the required python packages, binaries and libraries.

2. Download the latest version of the 
   `hiclib <https://bitbucket.org/mirnylab/hiclib/get/tip.zip>`_
   and unpack it into an empty folder. This folder will be called 
   "the hiclib folder" all throughout the document.
3. Run install_linux.py in the hiclib folder. This script modifies the .bashrc 
   and .bash_profile scripts and adds the the hiclib folder to the PYTHONPATH 
   environment variable. After restarting the terminal, python will be able to 
   locate the hiclib library.
4. Download `mirnylib <https://bitbucket.org/mirnylab/mirnylib/get/tip.zip>`_,
   unpack it into a separate folder and run install_linux.py.
   Failure during installation of mirnylib often means that cython is not up to date. 
5. Run download_bowtie.sh script to download and install the bowtie2 mapping
   software.
6. Run ./make_hg19.sh and ./make_sacCer3.sh scripts to download and index 
   the human and baker yeast genome. These scripts download the fasta files 
   with the nucleotide sequences and a gap file containing the positions of 
   the centromeres. The scripts for the other genomes are not 
   supplied and have to be created by a user. In most cases, it should be enough
   to change the GENOME_NAME variable of the make_hg19.sh script to the name 
   of the target genome, i.e. 'mm10' or 'bosTau7'.
7. Test the installation using the scripts from the /tests folders.

.. note:: If you're upgrading a package (numpy/scipy/matplotlib/cython/etc...) 
          using the pip Python package manager, never use "pip install PACKAGENAME",
          only "pip install --upgrade PACKAGENAME"! 

          If you have another version of the same package installed via the
          standard system package manager, e.g. apt-get or aptitude,
          the second command will most likely create a second copy of the package.
          This causes numerous problem that are very hard to detect.
          If 'pip install" says that package is already installed, remove it
          via the system package manager and only then try to install it again.
          However, if this asks you to remove a ton of dependent packages, 
          you will need to specifically delete files corresponding to the 
          original package. Current version and location of package files 
          can be determined from the python console by executing,
          e.g., >>>print numpy.__version__. and >>>print numpy.__file__
          
          Please refer to troubleshooting guide for a list of known errors! 
         
Hardware requirements
---------------------

Fragment-based analysis uses HDD to store all the information about the reads. 
However, at each point in time a few individual tracks are loaded to RAM to perform certain analysis. 
Memory requirements for this part can be estimated as 20-30 bytes per read, dependent on the number of filters used. 
For example, a 500 mln read Hi-C dataset will require no more than 16GB RAM. 
We're working on rewriting the core of the library to consume less memory. 

Binned data analysis uses RAM to store all the heatmaps.
Application of certain filters will create a copy of a heatmap, one at a time, even if multiple datasets are loaded. 
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
	For example ``>>> a = h5dict["myarray"] >>> a[0] = 1`` will not change the array on the hdd. You would need to run ``>>> h5dict["myarray"] = a`` to commit the change to the hard disk.  If you want to edit a part of an array, you can get a direct h5py dataset using ``>>> h5dict.get_dataset``

H5dict file can be converted to txt format or matlab format using a small utilities
h5dictToTxt.py and h5dictToMat.py. Both utilities are a part of mirnylib repository.

Pipeline structure
------------------

Two sides of the read are mapped to the genome independently using iterative mapping function. 
Acceptable read format is fastq (including fastq.gz), but .sra format can be incorporated using a bash_reader parameter 
of iterative_mapping and a fastq-dump tool from SRA toolkit. 
Iterative mapping outputs a combination of sam/bam files (autodetect by filename), which can be then parsed and written to an h5dict containing read data using parse_sams function. 

Fragment-level analysis inputs an h5dict, or any dictionary-like interface, containing at minimum information about read positions (however, strand information is highly recommended). 
It stores all the data in an h5dict, and modifies it as the analysis goes. 

.. warning:: fragmentHiC module modifies the data file when any filter is applied. 

An example script shows how  fragment-level analysis can be used to merge multiple datasets together, filter them and save heatmaps to an h5dict.

Heatmaps and vectors of SS reads are then exported by fragment-level analysis to an h5dict, that can be passe to a binnedData class. 
BinnedData class will load this h5dict, or any other dict-like object,
perform multiple types of iterative correction and eigenvector expansion, and compare results with different genomic tracks.


Example scripts
---------------

Our pipeline also presents a few example scripts that can be used to analyze the data. 
These scripts can be applied to any Hi-C data with minimal modification. 

To perform iterative mapping, use examples/iterative_mapping.py scipt. 

To combine data from multiple sequencing lanes, filter all lanes together (including filtering duplicates), and write out heatmaps. 
You can do it using examples/filteringManyDataests. 
You can modify genome files location inside the script, and provide mapping output files in the datasets.tsv file. 
By default, this script saves heatmaps for all replicas of an experiments, then for each experiment combines all replicas together and saves combined heatmaps.

To perform iterative correction and eigenvector expansion, use example scripts (example/iterativeCorrectionEigenvectorExpansion)


API documentation
-----------------

.. toctree::
   :maxdepth: 2
   
   glossary   
   genome
   mapping
   fragment 
   binneddata
   tutorial
   troubleshooting


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. raw:: html

    <hr width=50 size=10>


