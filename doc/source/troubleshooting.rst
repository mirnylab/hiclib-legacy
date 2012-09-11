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


