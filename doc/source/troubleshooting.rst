HiC-Lib troubleshooting guide
=============================

.. warning :: This document is quite wild. 

Problems with cython (numutils_new in mirnylib) 
-----------------------------------------------

**Error**  

File "<stdin>", line 1, in <module>

ImportError: ./numutils_new.so: undefined symbol: PyCapsule_New

**Solution** Please rebuild numuitls_new in mirnylib (using install_linux.py)

**Error** (while running install_linux.py )

Error converting Pyrex file to C:

cdef extern from "stdlib.h": 

    long c_libc_random "random"()
     
ctypedef fused my_type:

**Solution** 

upgrade cython (pip install --upgrade cython)


Problems with other packages
----------------------------

**problem** 
import of h5py prints this: 
ValueError: "numpy.dtype does not appear to be the correct type object"

Most likely arises on old python2.6 distributions

**solution** 
For me this was caused by a stray copy of numpy 1.3.1 and an installation of h5py built against that copy. 
I solved it by carefully uninstalling all copies of numpy, scipy, cython and h5py to the point when 
"import numpy" or "import scipy" or "import h5py" will all say "no module found..". 
Then I installed all packages using pip, and everything worked. 

**problem** 
installation of h5py using pip prints a lot of errors such as h5py/defs.c:25784: error: ‘__pyx_f_4h5py_4defs_H5Screate_simple’ undeclared (first use in this function)

**solution**
make sure that libhdf5-serial-dev is installed 

