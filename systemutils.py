import os, sys, cPickle, pdb 
import traceback

import ctypes 
from copy import copy






def _exceptionHook(infoType, value, tb):
    "sets exception hook to pdb"
    
    traceback.print_exception(infoType, value, tb)
    print
    pdb.post_mortem(tb) 

def setExceptionHook():
    sys.excepthook = _exceptionHook

def run_in_separate_process(func, *args, **kwds):
    pread, pwrite = os.pipe()
    pid = os.fork()
    if pid > 0:
        os.close(pwrite)
        with os.fdopen(pread, 'rb') as f:
            status, result = cPickle.load(f)
        os.waitpid(pid, 0)
        if status == 0:
            return result
        else:
            raise result
    else: 
        os.close(pread)
        try:
            result = func(*args, **kwds)
            status = 0
        except Exception, exc:
            result = exc
            status = 1
        with os.fdopen(pwrite, 'wb') as f:
            try:
                cPickle.dump((status,result), f, cPickle.HIGHEST_PROTOCOL)
            except cPickle.PicklingError, exc:
                cPickle.dump((2,exc), f, cPickle.HIGHEST_PROTOCOL)
        os._exit(0)


def _nprocessors():
    try:
        try:
            # Mac OS
            libc=ctypes.cdll.LoadLibrary(ctypes.util.find_library('libc'))
            v=ctypes.c_int(0)
            size=ctypes.c_size_t(ctypes.sizeof(v))
            libc.sysctlbyname('hw.ncpu', ctypes.c_voidp(ctypes.addressof(v)), ctypes.addressof(size), None, 0)
            return v.value
        except:
            # Cygwin (Windows) and Linuxes
            # Could try sysconf(_SC_NPROCESSORS_ONLN) (LSB) next.  Instead, count processors in cpuinfo.
            s = open('/proc/cpuinfo', 'r').read()
            return s.replace(' ', '').replace('\t', '').count('processor:')
    except:
        return 1

nproc = _nprocessors()

def fmap(f, *a, **kw):
    import struct
    builtin_map = map 
    """
    forkmap.map(..., n=nprocessors), same as map(...).
    n must be a keyword arg; default n is number of physical processors.
    """
    def writeobj(pipe, obj):
        s = cPickle.dumps(obj)
        s = struct.pack('l', -len(s)) + s
        os.write(pipe, s)

    def readobj(pipe):
        n = struct.unpack('l', os.read(pipe, 8))[0]
        s = ''
        an = abs(n)
        while len(s) < an:
            s += os.read(pipe, min(65536, an-len(s)))
        return cPickle.loads(s)

    n = kw.get('n', nproc)
    if n == 1:
        return builtin_map(f, *a)
    
    if len(a) == 1:
        L = a[0]
    else:
        L = zip(*a)
    try:
        len(L)
    except TypeError:
        L = list(L)
    n = min(n, len(L))
    
    ans = [None] * len(L)
    pipes = [os.pipe() for i in range(n-1)]
    
    for i in range(n):
        if i < n-1 and not os.fork():
        # Child, and not last processor
            try:
                try:
                    if len(a) == 1:
                        obj = builtin_map(f, L[i*len(L)//n:(i+1)*len(L)//n])
                    else:
                        obj = [f(*x) for x in L[i*len(L)//n:(i+1)*len(L)//n]]
                except Exception, obj:
                    pass
                writeobj(pipes[i][1], obj)
            except:
                traceback.print_exc()
            finally:
                os._exit(0)
        elif i == n-1:
            try: 
            
                if len(a) == 1:
                    ans[i*len(L)//n:] = builtin_map(f, L[i*len(L)//n:])
                else:
                    ans[i*len(L)//n:] = [f(*x) for x in L[i*len(L)//n:]]
                for k in range(n-1):
                    obj = readobj(pipes[k][0])
                    if isinstance(obj, Exception):
                        raise obj
                    ans[k*len(L)//n:(k+1)*len(L)//n] = obj
            finally:
                for j in range(n-1):
                    os.close(pipes[j][0])
                    os.close(pipes[j][1])
                    os.wait()
    return ans



def _fmapredcount(function,data,reduction = lambda x,y:x+y, n=4,exceptionList = [IOError]):
    """fork-map-reduce
    Performs fork-map of function on data, automatically reducing the data inside each worker. 
    If evaluation throws the exception from exceptionList, this results are simply ignored 
    """
    def funsum(x,y):
        """reduces two x[0],y[0], keeping track of # of
        successful evaluations that were made
        Also keeps track of None's that can occur if evaluation failed"""
        if x==None:
            if y==None:
                return None
            else:
                return y
        else:
            if y ==None:
                return x
            else:
                return (reduction(x[0],y[0]),x[1]+y[1])            
    
    def newfunction(x):
        "if function is evaluated, it was evaluated one time"
        return function(x),1
    
    if len(data) < n:
        n = len(data) 
    datas = []
    
    for i in xrange(n): 
        datas.append(copy(data[i::n]))   #split like that if beginning and end of the array have different evaluation time
    
    def worker(dataList):        
        try:
            dataList[0] = newfunction(dataList[0])
            return reduce(lambda z,y:funsum(z,newfunction(y)),dataList)  #reducing newfunction with our new reduction algorithm
        except tuple(exceptionList):
            return None
                     
    reduced = fmap(worker,datas,n=n)
    return reduce(funsum,reduced)

def fmapred(function,data,reduction = lambda x,y:x+y, n=4,exceptionList = [IOError]):
    """reduces two x[0],y[0], keeping track of # of
    successful evaluations that were made
    Also ignores failed evaluations with exceptions from exceptionList. 
    
    Parameters
    ----------
    function : function
        function to be applied to the data
    data : iterable
        input data
    reduction : function, optional
        Reduction function. By default - sum
    n : int, optional
        number of CPUs
    exceptionList : list, optional 
        list of exceptions to be ignored during reduction. By default, only IOError is ignored. 
    """
    return _fmapredcount(function,data,reduction = reduction, n=n,exceptionList = exceptionList)[0]

def fmapav(function,data,reduction = lambda x,y:x+y, n=4,exceptionList = [IOError]):
    """Calculates averate of [fucntion(i) for i in data]
    Also ignores failed evaluations with exceptions from exceptionList.
    
    Parameters
    ----------
    function : function
        function to be applied to the data
    data : iterable
        input data
    reduction : function, optional
        Reduction function. By default - sum
    n : int, optional
        number of CPUs
    exceptionList : list, optional 
        list of exceptions to be ignored during reduction. By default, only IOError is ignored. 
    """

    a = _fmapredcount(function,data,reduction = reduction, n=n,exceptionList = exceptionList)
    return a[0]/float(a[1])
    


