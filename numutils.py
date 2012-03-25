import numpy
na = numpy.array 
import  scipy.weave,scipy.sparse.linalg  
from scipy import weave 
import scipy.linalg

#-------------------------
"Mathematical utilities first" 
#-------------------------

def rank(x):
    "Returns rank of an array"
    tmp = numpy.argsort(x)
    return na(numpy.arange(len(x)),float)[tmp.argsort()]

def trunk(x,low=0.005,high = 0.005):
    "Truncates top 'low' fraction and top 'high' fraction of an array "    
    lowValue, highValue = numpy.percentile(x,[low*100.,(1-high)*100.])     
    return numpy.clip(x,a_min = lowValue, a_max = highValue)    
    
def arraySearch(array,tosearch):
    "returns location of tosearch in array; -->> assumes that elements exist!!! <--- " 
    inds = numpy.argsort(array)
    arSorted = array[inds]
    newinds = numpy.searchsorted(arSorted[:-1],tosearch)    
    return inds[newinds]

def arrayInArray(array,filterarray):    
    """gives you boolean array of indices of elements in array that are contained in filterarray
    a faster version of  [(i in filterarray) for i in array]"""           #sorted array
    array = numpy.asarray(array)  
    mask = numpy.zeros(len(array),'bool')   
    args = numpy.argsort(array)  
    arsort = array[args]
    diffs = numpy.r_[0,numpy.nonzero(numpy.diff(arsort) > 0.5)[0]+1,len(arsort)]  #places where sorted values of an array are changing
    values = arsort[diffs[:-1]]  #values at that places
    allinds = numpy.searchsorted(values[:-1],filterarray)   
    exist = values[allinds] == filterarray                 #check that value in filterarray exists in array
    N = len(allinds)    
    N,args,exist  #Eclipse warning remover 
    code = r"""
    #line 54 "numutils"
    using namespace std;
    for (int i = 0; i < N; i++)
    {    
        if (exist[i] == 0) continue;
        for (int j=diffs[allinds[i]];j<diffs[allinds[i]+1];j++)
        {
            mask[args[j]] = true;
        }
    } 
    """
    weave.inline(code, ['allinds', 'diffs' , 'mask' ,'args','N','exist'], 
                 extra_compile_args=['-march=native -malign-double -O3'],
                 support_code = r"#include <math.h>" )
    return mask

    
            
def arraySumByArray(array,filterarray,meanarray):
    "faster [sum(meanrrray [array == i]) for i in filterarray]"
    array = numpy.asarray(array,dtype = float)
    meanarray = numpy.asarray(meanarray,dtype = float)    
    if len(array) != len(meanarray): raise ValueError    
    args = numpy.argsort(array)
    arsort = array[args]
    diffs = numpy.r_[0,numpy.nonzero(numpy.diff(arsort) > 0.5)[0]+1,len(arsort)]
    values = arsort[diffs[:-1]]
    allinds = numpy.searchsorted(values[:-1],filterarray)
    exist = values[allinds] == filterarray
    N = len(allinds)
    args,exist,N #Eclipse warning removal 
    ret = numpy.zeros(len(allinds),meanarray.dtype)    
    code = """
    #line 50 "binary_search.py"
    using namespace std;
    for (int i = 0; i < N; i++)
    {
        if (exist[i] == 0) continue; 
        for (int j=diffs[allinds[i]];j<diffs[allinds[i]+1];j++)
        {
            ret[i] += meanarray[args[j]];
        }
    } 
    """
    support = """
    #include <math.h>  
    """
    weave.inline(code, ['allinds', 'diffs' , 'args' , 'ret','N','meanarray','exist'], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
    return ret


def sumByArray(array,filterarray, dtype = None):
    "faster [sum(array == i) for i in filterarray]"            
    arsort = numpy.sort(array)    
    diffs = numpy.r_[0,numpy.nonzero(numpy.diff(arsort) > 0.5)[0]+1,len(arsort)]
    if dtype != None: diffs = numpy.array(diffs, dtype = dtype)
    values = arsort[diffs[:-1]]
    del arsort        
    allinds = numpy.searchsorted(values[:-1],filterarray)
    notexist = values[allinds] != filterarray    
    del values
    c = diffs[allinds + 1] - diffs[allinds]
    c[notexist] = 0 
    return c


def trimZeros(x):
    "trims leading and trailing zeros of a 1D/2D array"
    if len(x.shape) == 1:
        nz = numpy.nonzero(x)[0]
        return x[nz.min():nz.max()+1]         
    ax1 = numpy.nonzero(numpy.sum(x,axis = 0))[0]    
    ax2 = numpy.nonzero(numpy.sum(x,axis = 1))[0]
    return x[ax1.min():ax1.max()+1, ax2.min() : ax2.max()+1 ]
    
def zoomOut(x,shape):
    "rescales an array preserving the structure and total sum"
    M1 = shape[0]
    M2 = shape[1]
    N1 = x.shape[0]
    N2 = x.shape[1]
    if (N1 < M1) or (N2 < M2):
        d1 = M1/N1 + 1
        d2 = M2/N2 + 1
        d = max(d1,d2)
        newx = numpy.zeros((d*N1,d*N2))        
        for i in xrange(d):
            for j in xrange(d):
                newx[i::d,j::d] = x/(1. * d**2)
        x = newx
        M1,M2,N1,N2 = M1*d,M2*d,N1*d,N2*d   #array is bigger now
    
    shift1 = N1/float(M1) + 0.000000001
    shift2 = N2/float(M2) + 0.000000001
    x = numpy.array(x,numpy.double,order = "C")
    tempres = numpy.zeros((M1,N2),float)    
    for i in xrange(N1):
        beg = (i/shift1)
        end = ((i+1)/shift1)
        if int(beg) == int(end): 
            tempres[beg,:] += x[i,:]
        else:
            tempres[beg,:] += x[i,:] * (int(end) - beg) / (end - beg)
            tempres[beg+1,:] += x[i,:] * (end - int(end)) / (end - beg)
    res = numpy.zeros((M1,M2),float)
    for i in xrange(N2):
        beg = (i/shift2)
        end = ((i+1)/shift2)
        if int(beg) == int(end): 
            res[:,beg] += tempres[:,i]
        else:
            res[:,beg] += tempres[:,i] * (int(end) - beg) / (end - beg)
            res[:,beg+1] += tempres[:,i] * (end - int(end)) / (end - beg)
    return res

smartZoomOut = zoomOut  #backwards compatibility

def coarsegrain(array,size):
    "coarsegrains array by summing values in sizeXsize squares; truncates the last square"
    if len(array.shape) == 2:
        N = len(array) - len(array) % size 
        array = array[:N,:N]
        a = numpy.zeros((N/size,N/size),float)
        for i in xrange(size):
            for j in xrange(size):
                a += array[i::size,j::size]
        return a
    if len(array.shape) == 1:
        array = array[:(len(array) / size) * size]
        narray = numpy.zeros(len(array)/size,float)
        for i in xrange(size):
            narray += array[i::size]
        return narray

#-----------------------------------
"Iterative correction, PCA and other Hi-C related things"
#-----------------------------------

def maskPCA(A,mask):
    "attempts to perform PCA-like analysis of an array with a masked part"    
    from numpy import linalg
    mask = numpy.array(mask, int)
    bmask = mask == 0
    A[bmask] = 0  
    sums = numpy.sum(A,axis = 0)
    means = sums / (1. * numpy.sum(mask, axis = 0))     
    M = (A - means).T
    M[bmask]  = 0
    covs = numpy.zeros((len(M),len(M)),float)
    for i in xrange(len(M)):
        myvector = M[i]
        mymask = mask[i]
        allmask = mask * mymask[None,:]
        tocov = myvector[None,:] * M 
        covsums = tocov.sum(axis = 1)
        masksums = allmask.sum(axis = 1)
        covs[i] = covsums/masksums
    [latent,coeff] = linalg.eig(covs)
    print latent[:4]
    return coeff


def PCA(A):
    "performs PCA analysis, and returns 6 best principal components"
    A = numpy.array(A,float)
    M = (A-numpy.mean(A.T,axis=1)).T 
    covM = numpy.dot(M,M.T)
    [latent,coeff] =  scipy.sparse.linalg.eigsh(covM,6)
    print latent
    return coeff[:,::-1]


def EIG(A):
    "Performs mean-centered engenvector expansion"
    A = numpy.array(A,float)    
    M = (A - numpy.mean(A)) # subtract the mean (along columns)
    if (M -M.T).var() < numpy.var(M[::10,::10]) * 0.000001:
        [latent,coeff] = scipy.sparse.linalg.eigsh(M,3)
        print "herm"   #matrix is hermitian
    else: 
        [latent,coeff] = scipy.sparse.linalg.eigs(M,3)
        print "norm"   #Matrix is normal
    alatent = numpy.argsort(numpy.abs(latent)) 
    print latent[:4]
    coeff = coeff[:,alatent]
    return coeff[:,::-1]

def project(data,vector):
    "project data on a single vector"
    dot = numpy.sum((data*vector[:,None]),axis = 0)     
    den = (vector * vector).sum()
    return vector[:,None] * (dot / den)[None,:]


def projectOnEigenvalues(data,N=1):
    "projects symmetric data on the first N eigenvalues"
    #TODO: rewrite properly for both symmetric and non-symmetric case 
    meanOfData = numpy.mean(data)
    mdata = data - meanOfData
    symData = 0.5*(mdata + mdata.T)    
    values,vectors = scipy.linalg.eig(symData)
    ndata = 0
    for i in xrange(N):
        ndata += values[i] * vectors[:,i][:,None] * vectors[:,i][None,:]
    return ndata + meanOfData 
    

def correct(y):
    "Correct non-symmetric or symmetirc data once"
    x = numpy.array(y,float)        
    s = numpy.sum(x,axis = 1)
    s /= numpy.mean(s[s!=0])    
    s[s==0] = 1     
    s2 = numpy.sum(x,axis = 0)        
    s2 /= numpy.mean(s2[s2!=0])
    s2[s2==0] = 1
    return x / (s2[None,:] * s[:,None])

def correctInPlace(x):
    "works for non-symmetric and symmetric data"            
    s = numpy.sum(x,axis = 1)
    s /= numpy.mean(s[s!=0])    
    s[s==0] = 1     
    s2 = numpy.sum(x,axis = 0)        
    s2 /= numpy.mean(s2[s2!=0])
    s2[s2==0] = 1    
    x /= (s2[None,:] * s[:,None])


def ultracorrectSymmetricWithVector(x,v = None,M=50,diag = -1):
    """Main method for correcting DS and SS read data. Possibly excludes diagonal.
    By default does iterative correction, but can perform an M-time correction"""    
    totalBias = numpy.ones(len(x),float)
    code = """
    #line 288 "numutils" 
    using namespace std;
    for (int i = 0; i < N; i++)    
    {    
        for (int j = 0; j<N; j++)
        {
        x[N * i + j] = x [ N * i + j] / (s[i] * s[j]);
        }
    } 
    """
    if v == None: v = numpy.zeros(len(x),float)  #single-sided reads
    x = numpy.array(x,float,order = 'C')
    v = numpy.array(v,float,order = "C")
    N = len(x)
    N #Eclipse warning remover     
    for _ in xrange(M):         
        s0 = numpy.sum(x,axis = 1)         
        mask = [s0 == 0]            
        v[mask] = 0   #no SS reads if there are no DS reads here        
        nv = v / (totalBias * (totalBias[mask==False]).mean())
        s = s0 + nv
        for dd in xrange(diag + 1):   #excluding the diagonal 
            if dd == 0:
                s -= numpy.diagonal(x)
            else:
                dia = numpy.diagonal(x,dd)
                #print dia
                s[dd:] -= dia
                s[:-dd] -= dia 
        s /= numpy.mean(s[s0!=0])
        s[s0==0] = 1
        totalBias *= s
        scipy.weave.inline(code, ['x','s','N'], extra_compile_args=['-march=native -malign-double -O3'])  #performing a correction
    corr = totalBias[s0!=0].mean()  #mean correction factor
    x *= corr * corr #renormalizing everything
    totalBias /= corr
    return x,v/totalBias 



def ultracorrectSymmetricByMask(x,mask,M = 50):
    """performs iterative correction excluding some regions of a heatmap from consideration.
    These regions are still corrected, but don't contribute to the sums 
    """
    code = """
    #line 333 "numutils.py"
    using namespace std;    
    for (int i=0;i<N;i++)
    {
        for (int j = 0;j<N;j++)
        {
        if (mask[N * i + j] == 1)
            {
            sums[i] += x[N * i + j];
            counts[i] += 1;             
            }
        }
    }
    float ss = 0;
    float count = 0; 
    for (int i = 0;i<N;i++)
    {    
        if ((counts[i] > 0) && (sums[i] > 0))
        {
            sums[i] = sums[i] / counts[i];
            ss+= sums[i];             
            count+= 1;           
        }
        else
        {
            sums[i] = 1; 
        }
    }
    for (int i = 0;i<N;i++)
    {
        sums[i] /= (ss/count);
    }
    for (int i = 0; i < N; i++)    
    {    
        for (int j = 0; j<N; j++)
        {
        x[N * i + j] = x [ N * i + j] / (sums[i] * sums[j]);
        }
    } 
    """
    support = """
    #include <math.h>  
    """    
    x = numpy.asarray(x,float,order = 'C')
    N = len(x)
    allsums = numpy.zeros(N,float)
    for i in xrange(M):
        sums = numpy.zeros(N,float)
        counts = numpy.zeros(N,float)
        i,counts  #Eclipse warning removal     
        scipy.weave.inline(code, ['x','N','sums','counts','mask'], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
        allsums *= sums
    return x,allsums


def ultracorrect(x,M=20):
    "just iterative correction of symmetric matrix"
    x = numpy.array(x,float)
    print numpy.mean(x),
    newx = numpy.array(x)
    for _ in xrange(M):         
        correctInPlace(newx)
    print numpy.mean(newx),
    newx /= (1. * numpy.mean(newx)/numpy.mean(x))
    print numpy.mean(newx)
    return newx

def correctBias(y):
    "performs single correction and returns data + bias"
    x = numpy.asarray(y,dtype=float)        
    s = numpy.sum(x,axis = 1)
    s /= numpy.mean(s[s!=0])    
    s[s==0] = 1     
    return x / (s[None,:] * s[:,None]),s 

def ultracorrectBiasReturn(x,M=20):
    "performs iterative correction and returns bias"
    x = numpy.array(x,float)
    print numpy.mean(x),
    newx = numpy.array(x)
    ball = numpy.ones(len(x),float)
    for _ in xrange(M):
        newx,b = correctBias(newx)
        ball *= b
    print numpy.mean(newx),
    newx /= (1. * numpy.mean(newx)/numpy.mean(x))
    print numpy.mean(newx)
    return newx,ball

def create_regions(a):
    "creates array of nonzero regions"    
    a = numpy.array(a,int)
    a = numpy.concatenate([numpy.array([0],int),a,numpy.array([0],int)])
    a1 = numpy.nonzero(a[1:] * (1-a[:-1]))[0]
    a2 = numpy.nonzero(a[:-1] * (1-a[1:]))[0]    
    return numpy.transpose(numpy.array([a1,a2]))
    
