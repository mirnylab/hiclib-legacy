from libc.stdint cimport int32_t, int64_t, uint64_t

import numpy as np
cimport numpy as np 

cdef extern from "mycode.h":
    void int64BinarySearch(uint64_t* data, long int N, uint64_t* keys, long int M, int32_t* result,uint64_t mult)
    


def binarySearch(_data, _keys):
    maxValue = int(_keys.max())
    a = 2**63
    cdef uint64_t mymultiplyer  = a / (maxValue + 1000) - 1000 
    
    
    cdef np.ndarray[np.uint64_t, ndim=1, mode = 'c'] data = np.ascontiguousarray(_data, dtype = np.uint64)
    cdef np.uint64_t* data_pointer = <uint64_t*> data.data
    
    cdef np.ndarray[np.uint64_t, ndim=1, mode = 'c'] keys = np.ascontiguousarray(_keys, dtype = np.uint64)
    cdef np.uint64_t* keys_pointer = <uint64_t*> keys.data

    
    cdef np.ndarray[np.int32_t, ndim=1, mode = 'c'] result = np.empty(len(_data), dtype = np.int32)
    cdef np.int32_t* result_pointer = <int32_t*> result.data
    int64BinarySearch(data_pointer, len(_data), keys_pointer, len(keys), result_pointer,mymultiplyer)
    return result

"""
def binarySearch(_data, _keys):
    
    cdef np.ndarray[np.uint64_t, ndim=1, mode = 'c'] data = np.ascontiguousarray(_data, dtype = np.uint64)
    cdef np.uint64_t* data_pointer = <uint64_t*> data.data

    cdef np.ndarray[np.uint64_t, ndim=1, mode = 'c'] keys = np.ascontiguousarray(_keys, dtype = np.uint64)
    cdef np.uint64_t* keys_pointer = <uint64_t*> keys.data
    
    cdef np.ndarray[np.int32_t, ndim=1, mode = 'c'] result = np.empty(len(_data), dtype = np.int32)
    cdef np.int32_t* result_pointer = <int32_t*> result.data
    uint64BinarySearch(data_pointer, len(_data), keys_pointer, len(keys), result_pointer, mymultiplyer)
    return result
"""
