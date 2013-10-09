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

Problems with h5dict
--------------------

A recently discovered problem is caused by new version of h5py (that is 2.1, released Oct 9). 
We're trying to fix this... in the meantime just don't upgrade to it. 

> Checking for mirnylib.h5dict install..
> Traceback (most recent call last):
>   File "install_linux.py", line 52, in <module>
>     a["numpy"] = b
>   File
> "/net/gs/vol3/software/modules-sw/python/2.7.2/Linux/RHEL6/x86_64/lib/python2.7/site-packages/mirnylib/h5dict.py",
> line 168, in __setitem__
>     self.__self_dump__()
>   File
> "/net/gs/vol3/software/modules-sw/python/2.7.2/Linux/RHEL6/x86_64/lib/python2.7/site-packages/mirnylib/h5dict.py",
> line 98, in __self_dump__
>     data=cPickle.dumps(data, protocol= -1))
>   File "build/bdist.linux-x86_64/egg/h5py/_hl/group.py", line 71, in
> create_dataset
>   File "build/bdist.linux-x86_64/egg/h5py/_hl/dataset.py", line 94, in
> make_new_dset
>   File "h5d.pyx", line 221, in h5py.h5d.DatasetID.write (h5py/h5d.c:2866)
>   File "_proxy.pyx", line 152, in h5py._proxy.dset_rw (h5py/_proxy.c:1766)
>   File "defs.pyx", line 1780, in h5py.defs.H5Tconvert (h5py/defs.c:18439)
>   File "_errors.pyx", line 91, in h5py._errors.set_exception
> (h5py/_errors.c:711)
>   File "_conv.pyx", line 430, in h5py._conv.str2vlen (h5py/_conv.c:3382)
>   File "_conv.pyx", line 116, in h5py._conv.generic_converter
> (h5py/_conv.c:1302)
>   File "_conv.pyx", line 253, in h5py._conv.conv_str2vlen
> (h5py/_conv.c:2158)
> ValueError: VLEN strings do not support embedded NULLs
>




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

