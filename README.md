hiclib
======
 

Notice 6/14/19
------------
This repository is not actively maintained. 

Please see cooler (https://github.com/mirnylab/cooler/), cooltools (https://github.com/mirnylab/cooltools), and distiller (https://github.com/mirnylab/distiller-nf) for actively maintained & developed Hi-C analysis tools.
 
 
Read the [documentation](https://mirnylab.bitbucket.io/hiclib/).

Installation
------------
The easiest way is to use pip.

First, install [mirnylib](https://bitbucket.org/mirnylab/mirnylib), e.g:

`$ pip install https://bitbucket.org/mirnylab/mirnylib/get/tip.tar.gz`

Then install hiclib:

`$ pip install https://bitbucket.org/mirnylab/hiclib/get/tip.tar.gz`

Or clone the repo and install from there. Alternatively, there is a install script that doesn't use pip (Linux only).

Installation requirements
-------------------------

### Python dependencies

Python 2.7/3.3+

Getting the basic Scientific Python stack can be trickier on some platforms than others. For more details, see the [instructions on scipy.org](http://www.scipy.org/install.html). You should already have these dependencies installed and running correctly before attempting to install this package.

- numpy (1.6+)
- scipy
- matplotlib
- h5py
- cython
- numexpr
- statsmodels

We highly recommend using the [conda](http://conda.pydata.org/miniconda.html) package/environment manager if you have trouble building the core scientific Python packages.

`$ conda install numpy scipy matplotlib h5py cython ...`

Other deps:

- biopython
- pysam

Optional:

- Sphinx (to generate documentation)

### Non-python dependencies

- Mirnylib provides a HDF5 dictionary type based on the Python h5py package. For the h5py package to install properly, you need to have the shared library and development headers for HDF5 1.8.4 or newer installed (`libhdf5-dev` or similar). See the [h5py docs](http://docs.h5py.org/en/latest/build.html) for more information.

- You will also need the `bowtie2` mapping software which can be downloaded manually or installed via the download_bowtie.sh script in the hiclib folder.