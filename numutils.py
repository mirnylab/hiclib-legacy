import base   
from math import sin,cos,exp,log,sqrt 
import numpy
na = numpy.array
from numpy import r_
import scipy.sparse, scipy.weave  
 
import hashlib 
import random
 
#import mdp 


  
    

def rank(x):
    x = numpy.array(x)
    tmp = x.argsort()
    return na(numpy.arange(len(x)),float)[tmp.argsort()]

def trunk(x,low=0.005,high = 0.005):
    x  = numpy.array(x)
    sortedArray = numpy.sort(x.ravel())
    NS = len(sortedArray)
    lowValue = sortedArray[int(NS* (low))]
    highValue = sortedArray[int(NS* (1 - high))]
    x[x>highValue] = highValue
    x[x<lowValue] = lowValue 
    return x

def maskPCA(A,mask):    
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
        #mat_img(allmask)
        tocov = myvector[None,:] * M 
        #mat_img(tocov)
        covsums = tocov.sum(axis = 1)
        masksums = allmask.sum(axis = 1)
        #print covsums
        #print masksums
        covs[i] = covsums/masksums
    [latent,coeff] = linalg.eig(covs)
    print latent[:4]
    return coeff
    





def PCA(A):

    A = numpy.array(A,float)
    """ performs principal components analysis 
        (PCA) on the n-by-p data matrix A
        Rows of A correspond to observations, columns to variables. 
    
    Returns :  
     coeff :
       is a p-by-p matrix, each column containing coefficients 
       for one principal component.
     score : 
       the principal component scores; that is, the representation 
       of A in the principal component space. Rows of SCORE 
       correspond to observations, columns to components.
    
     latent : 
       a vector containing the eigenvalues 
       of the covariance matrix of A.
    """
# computing eigenvalues and eigenvectors of covariance matrix
    M = (A-numpy.mean(A.T,axis=1)).T # subtract the mean (along columns)
    covM = numpy.dot(M,M.T)
    [latent,coeff] =  scipy.sparse.linalg.eigsh(covM,6)
    print latent
    return coeff[:,::-1]


def EIG(A):

    A = numpy.array(A,float)
# computing eigenvalues and eigenvectors of covariance matrix
    M = (A - numpy.mean(A)) # subtract the mean (along columns)
    if (M -M.T).var() < numpy.var(M[::10,::10]) * 0.000001:
        [latent,coeff] = scipy.sparse.linalg.eigsh(M,3)
        print "herm"
    else: 
        [latent,coeff] = scipy.sparse.linalg.eigs(M,3)
        print "norm"
    alatent = numpy.argsort(numpy.abs(latent)) 
    print latent[:4]
    coeff = coeff[:,alatent]
    return coeff[:,::-1]





#performICA(data,4)

def projectOnEigenvalues(data,N=1):
    meanOfData = numpy.mean(data)
    mdata = data - meanOfData
    symData = 0.5*(mdata + mdata.T)
    import scipy.linalg
    values,vectors = scipy.linalg.eig(symData)
    ndata = 0
    for i in xrange(N):
        ndata += values[i] * vectors[:,i][:,None] * vectors[:,i][None,:]
    return ndata + meanOfData 
    
    
    

    


def correct(y):
    x = numpy.array(y,float)        
    s = numpy.sum(x,axis = 1)
    s /= numpy.mean(s[s!=0])    
    s[s==0] = 1     
    s2 = numpy.sum(x,axis = 0)        
    s2 /= numpy.mean(s2[s2!=0])
    s2[s2==0] = 1
    return x / (s2[None,:] * s[:,None])

def correctInPlace(x):            
    s = numpy.sum(x,axis = 1)
    s /= numpy.mean(s[s!=0])    
    s[s==0] = 1     
    s2 = numpy.sum(x,axis = 0)        
    s2 /= numpy.mean(s2[s2!=0])
    s2[s2==0] = 1    
    x /= (s2[None,:] * s[:,None])
    




def ultracorrectSymmetricWithVector(x,v = None,M=50,chromosomes = None,diag = -1):
    totalBias = numpy.ones(len(x),float)
    code = """
    #line 50 "binary_search.py"
    using namespace std;
    for (int i = 0; i < N; i++)    
    {    
        for (int j = 0; j<N; j++)
        {
        x[N * i + j] = x [ N * i + j] / (s[i] * s[j]);
        }
    } 
    """
    support = """
    #include <math.h>  
    """
    
    x = numpy.array(x,float,order = 'C')
    if v == None: v = numpy.zeros(len(x),float)

    N = len(x)
    N #Eclipse warning remover 
    v = numpy.array(v,float,order = "C")
    if chromosomes != None: 
        chromosomes = r_[chromosomes,len(v)]
         
    for _ in xrange(M):
         
        s0 = numpy.sum(x,axis = 1)
         
        mask = [s0 == 0]            
        v [s0 == 0] = 0
        
        nv = v / (totalBias * (totalBias[mask==False]).mean())
        
         
        
        
         
        s = s0 + nv
          
    
        for dd in xrange(diag + 1):
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
        #print s 
        scipy.weave.inline(code, ['x','s','N'], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
        
        

    #x = x / (s1/s_in)
    #v = v / (s1/s_in)     

    corr = totalBias[s0!=0].mean()
    x *= corr * corr
    totalBias /= corr
    return x,v/totalBias 


#a = numpy.array([[3,2,1],[2,4,5],[1,5,6]])
#b = numpy.array([7,8,9])
#biases = numpy.random.random(3) + 0.5
#biases /= biases.mean()
#ab = a * biases[:,None] * biases[None,:]
#bb = b * biases
#
#
#ac,bc = ultracorrectSymmetricWithVector(ab, bb,100)
#nbiases =  bb/bc
#print nbiases
#
#print ab   - ac * ((nbiases[:,None] * nbiases[None,:]))
#print bb   - bc * nbiases * nbiases.mean()  

 



def ultracorrectSymmetricByMask(x,mask,M = 50):
    code = """
    #line 50 "binary_search.py"
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
    
    x = numpy.array(x,float,order = 'C')
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
    x = numpy.array(y,float)        
    s = numpy.sum(x,axis = 1)
    s /= numpy.mean(s[s!=0])    
    s[s==0] = 1     
    return x / (s[None,:] * s[:,None]),s 

def ultracorrectBiasReturn(x,M=20):
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



  

def project(data,vector):
    dot = numpy.sum((data*vector[:,None]),axis = 0)
     
    den = (vector * vector).sum()
    return vector[:,None] * (dot / den)[None,:]









def coolAverage(*data):
    if len(data) ==1: data = data[0]    
    data = [numpy.array(i) for i in data]
    
        
    shapes = [i.shape for i in data]
    for i in shapes:
        if len(i) != len(shapes[0]):
            print "shape mismatch"
            return -1                  
    nshapes = numpy.array([numpy.array(i) for i in shapes])
    mshape = numpy.max(nshapes,axis = 0)
    newshape = tuple([len(data)] + list(mshape))    
    newarray = numpy.zeros(newshape,float)
    mask = numpy.zeros(newshape,"bool")
    for i in xrange(len(data)):
        myslice = [slice(0,end,None) for end in shapes[i]]
        newarray[i][myslice] = data[i]
        mask[i][myslice] = True
    sums  = numpy.sum(newarray,axis = 0)
    dividers = numpy.sum(mask,axis = 0)
    return sums / dividers


def arraySearch(array,tosearch):
    "returns location of tosearch in array" 
    inds = numpy.argsort(array)
    arSorted = array[inds]
    newinds = numpy.searchsorted(arSorted[:-1],tosearch)    
    return inds[newinds]







def arrayInArray(array,filterarray):
    "gives you indices in array that are contained in filterarray"
    
        
    from scipy import weave
    arsort = numpy.sort(array)
    mask = numpy.zeros(len(array),'bool')
    args = numpy.argsort(array)
    diffs = numpy.r_[0,numpy.nonzero(numpy.diff(arsort) > 0.5)[0]+1,len(arsort)]
    #print arsort[diffs[0]:diffs[1]]
    values = arsort[diffs[:-1]]
    allinds = numpy.searchsorted(values[:-1],filterarray)
    exist = values[allinds] == filterarray
    N = len(allinds)    
    N,args,exist  #used below, warning remover 
    code = """
    #line 50 "binary_search.py"
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
    support = """
    #include <math.h>  
    """
    weave.inline(code, ['allinds', 'diffs' , 'mask' ,'args','N','exist'], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
    return mask
        


        
def trimZeros(x):
    if len(x.shape) == 1:
        nz = numpy.nonzero(x)[0]
        a = nz.min()
        b = nz.max() 
        return x[a:b+1] 
    print x.shape
    s1 = numpy.sum(x,axis = 0) > 0
    a = numpy.nonzero(s1 == True)[0]
    print numpy.min(a),numpy.max(a),
    s1[numpy.min(a) : numpy.max(a)] = True    
    s2 = numpy.sum(x,axis = 1) > 0
    a = numpy.nonzero(s2 == True)[0]
    print numpy.min(a),numpy.max(a)
    s2[numpy.min(a) : numpy.max(a)] = True
        
     
    
    y = x[s2][:,s1]
    return y 
    
def zoomOut(x,shape):
    M1 = shape[0]
    M2 = shape[1]
    N1 = x.shape[0]
    N2 = x.shape[1]
    if (N1 < M1) or (N2 < M2): raise
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

def smartZoomOut(x,shape):
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
        return zoomOut(newx,shape)
    else: return zoomOut(x,shape)
    
    

    

def arraySumByArray(array,filterarray,meanarray):
    "return [ sum(meanrrray [array == i]) for i in filterarray]"
    array = numpy.array(array,dtype = float)
    meanarray = numpy.array(meanarray,dtype = float)
    if len(array) != len(meanarray): raise 
    from scipy import weave
    arsort = numpy.sort(array)
    args = numpy.argsort(array)
    diffs = numpy.r_[0,numpy.nonzero(numpy.diff(arsort) > 0.5)[0]+1,len(arsort)]
    # print arsort[diffs[0]:diffs[1]]
    values = arsort[diffs[:-1]]
    means = []
    allinds = numpy.searchsorted(values[:-1],filterarray)
    exist = values[allinds] == filterarray

    #print values[allinds],filterarray
    N = len(allinds)
    args,means,exist,N #used below, warning remover 
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
    "return [ sum(array == i) for i in filterarray]"            
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



    


def create_regions(a):
    "creates array of nonzero regions"
    
    a = numpy.array(a,int)
    a = numpy.concatenate([numpy.array([0],int),a,numpy.array([0],int)])
    a1 = numpy.nonzero(a[1:] * (1-a[:-1]))[0]
    a2 = numpy.nonzero(a[:-1] * (1-a[1:]))[0]
    
    return numpy.transpose(numpy.array([a1,a2]))
    
    
def fill_sphere(r):
    points = []
    def positions(x):
        N = int(x)
        return [x * i / (N+1) for i in xrange(N+1)]
    OneLen = numpy.pi * r
    heights = positions(OneLen)
    angles = [i/r for i in heights]
    print angles
    for i in angles:
        len2 = 2 * numpy.pi * sin(i) * r
        if abs(len2) < 0.5: angles2 = [0]
        else: angles2 = [k / (r * sin(i)) for k in positions(len2)]
        
        for j in angles2:
            z = r * cos(i)
            xy = r * sin(i)
            x = xy * cos(j)
            y = xy * sin(j)
            points.append((x,y,z))
    return numpy.array(points)
    
            
        


 
def exprandom(a,b):
    A = exp(a)
    B = exp(b)
    ra = numpy.random.random()
    ra2 = ra*(B-A) + A 
    return log(ra2)

   
class harray(numpy.ndarray):
    "hashable numpy array, do not change it on the way"
    def __new__(self,subtype, data):
        subarr = numpy.array(data)
        subarr = subarr.view(subtype)
        return subarr

        
    def __hash__(self):
        try:
            self.hashed
            
        except:
            self.hashed = int(hashlib.md5(self.tostring()).hexdigest(),16)
        return self.hashed
    def __eq__(self,a):
        return self.__hash__() == a.__hash__()

def realContacts(a):
    a = a>0
    while True:
        b = numpy.dot(a,a)
        b = b>0       
        if numpy.sum(b) == numpy.sum(a):
            return a
        a = b
        
def contactPower(a,N): 
    a = a>0
    for _ in xrange(N):
       
        b = numpy.dot(a,a)
        b = b>0       
        a = b
    return a 


def coarsegrain(array,size):
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



#from inout import img_show
#img_show(coarsegrain(numpy.random.random((1000,1000)),13))    

def corr2d(x):
    x = numpy.array(x)
    t = numpy.fft.fft2(x)
    return numpy.real(numpy.fft.ifft2(t*numpy.conjugate(t)))





def logbins(a, b, pace, N_in=0):
    "create log-spaced bins"
    beg = log(a)
    end = log(b - 1)
    pace = log(pace)
    N = int((end - beg) / pace)
    if N_in != 0: N = N_in  
    pace = (end - beg) / N
    mas = numpy.arange(beg, end + 0.000000001, pace)
    ret = numpy.exp(mas)
    ret = numpy.array([int(i) for i in ret])
    ret[-1] = b 
    for i in xrange(len(ret) - 1):
        if ret[i + 1] <= ret[i]:
            ret[i + 1: - 1] += 1
    return [int(i) for i in ret]



def alogbins(a,b,pace,N_in=0):
    "create floating point log bins"
    beg = log(a)
    end = log(b - 1)
    pace = log(pace)
    N = int((end - beg) / pace)
    if N_in != 0: N = N_in  
    pace = (end - beg) / N
    mas = numpy.arange(beg, end + 0.000000001, pace)
    ret = numpy.exp(mas)
    return ret

def rescale(data):
    "rescales array to zero mean unit variance"
    if type(data) !=numpy.ndarray: data = numpy.array(data,'float')
    L = len(data)
    av = numpy.sum(data)/float(L)
    fun = lambda x:(x-av)**2
    dis = sqrt(numpy.sum(fun(data))/float(L))+0.000001
    f = lambda x:(x-av)/dis
    return numpy.array(f(data),'float')
    
def autocorr(x):
    x = rescale(x)
    result = numpy.correlate(x, x, mode='full')
    return result[result.size/2:]

def rotationMatrix(theta):
    tx,ty,tz = theta
    
    Rx = numpy.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
    Ry = numpy.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0, cos(ty)]])
    Rz = numpy.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0,0,1]])
    
    return numpy.dot(Rx, numpy.dot(Ry, Rz))


def random_on_sphere(r=1):
    while True:        
        a = numpy.random.random(2);
        x1 = 2*a[0]-1
        x2 = 2*a[1]-1
        if x1**2 + x2**2 >1: continue
        t = sqrt(1 - x1**2 - x2**2)
        x = r*2*x1* t
        y = r*2* x2 * t
        z = r*(1 - 2*(x1**2 + x2**2))
        return (x,y,z)

def random_in_sphere(r=1):
    while True:
        a = numpy.random.random(3)*2-1
        if numpy.sum(a**2) < 1:
            return r*a
randomInSphere = random_in_sphere



class node():
    "representation of a tree"
    def __init__(self,parent = None ,generation = 0 ):
        self.parent = parent
        self.children = []
        self.generation = generation
        self.childCount = 0
    def createChildren(self,N):
        
        for _ in xrange(N):
            child = node(self)
            child.generation = self.generation + 1 
            self.children.append(child)
        return self.children
    def add_child(self):        
        self.children.append(node(self))
        self.children[-1].generation = self.generation + 1 
    def giveAllChildren(self):
        children  = sum([i.giveAllChildren()    for i in self.children],[]) + [self]
        self.childCount = len(children)
        return children
    def giveCOM(self):
        coms  = [i.giveCOM()    for i in self.children]
        com,com2 = self.pos.copy(),self.pos**2
        
        for i in coms:
            com += i[0]
            com2 += i[1]
        #print com,com2,self.pos
        self.rg = sqrt(((-com**2/float(self.childCount) + com2)/float(self.childCount)).sum())
        return  com,com2
    def createLambdaChildren(self,L,maxgen):
        if (self.generation) == maxgen:
            return
        x = numpy.random.poisson(L)        
        self.createChildren(x)
        for i in self.children:
            i.createLambdaChildren(L,maxgen)
def percentileMask(array,low=0.01,high=0.01):
    sortarray = numpy.sort(array)
    N = len(array)
    low = sortarray[low * N]
    high = sortarray[(1-high)*N]
    return (array > low) * (array < high)
class coolnode(node):
    def add_child(self):
        try:self.children.append(coolnode(self))
        except: 
            print self.generation
            return -1 
            
        self.children[-1].generation = self.generation + 1 
        
    def cool_add_children(self,fun,children):
        for i in children:
            if self.add_child() != -1: 
                self.children[-1].pos = i
            else:
                return  
        for i in self.children:
            i.cool_add_children(fun,fun(i.pos))
    def rewire(self,parent):
        
        if parent ==-1:
            if self.parent == None: self.parent = parent 
            self.children = [i for i in (self.children + [self.parent] ) if i != parent]
            self.generation = 0
            self.parent = None 
            
            
        elif self.parent != None: 
            self.children = [i for i in (self.children + [self.parent] ) if i != parent]
            self.generation = parent.generation + 1
            self.parent = parent 
        else: 
            self.children = [i for i in (self.children ) if i != parent]
         
            self.generation = parent.generation + 1
            self.parent = parent
             
        for i in self.children: i.rewire(self)
    def branchCount(self):
        if type(self.parent) == type(self):
            return len(self.children) + 1
        return len(self.children)

class nodes:
    def __init__(self,parent):
        self.parent = parent
        assert isinstance(parent,(node,coolnode))
        self.nodes = parent.giveAllChildren()
        
        
    def give_generations(self):
        return [i.generation for i in self.nodes]
        
    def rewire(self,new_parent):
        self.parent = new_parent
        self.parent.rewire(-1)
        
    def max_rewire(self):
        
            
        node = random.choice(self.nodes)
        node.generation = -1 
         
        self.rewire(node)
        generations1 = na(self.give_generations())
        node = self.nodes[numpy.argmax(generations1)]
        self.rewire(node)
        generations1 = na(self.give_generations())
        node = self.nodes[numpy.argmax(generations1)]
        self.rewire(node)        
        na(self.give_generations())   

            

        

def createLambdaTree(N,l):
    "creates tree with branching probability l"
    mynodes = [coolnode() for i in xrange(N)]
    for j,i in  enumerate(mynodes):
        r = numpy.random.random() 
        if r<l:
            i.tri = 3
        else:
            i.tri = 2
        i.ID = j 
    
    while True:
        #print [i.ID for i in nodes]
        #time.sleep(0.1) 
        a,b = random.choice(mynodes),random.choice(mynodes)
        if a.ID == b.ID:
            continue
        if (a.branchCount() < a.tri) and (b.branchCount() < b.tri):
            if (a.parent != b) and type(a.parent) == type(a):
                a.children.append(a.parent) 
            a.parent = b
            b.children.append(a)
            a.rewire(-1)
            
            myid = a.ID
            t = a.giveAllChildren()
            #print len(t)  
            if len(t) == N:
                allNodes  = nodes(a)
                allNodes.max_rewire()
                return allNodes.parent.giveAllChildren()  
            for i in t: i.ID = myid
    
    

#bins = logbins(10,15000,4)

def onetree(treeBin):
    def mytree(x):
        numpy.random.seed(x)
        random.seed(x) 
        return numpy.mean([i.generation for i in createLambdaTree(treeBin,1)])    
    return numpy.mean(na(base.fmap(mytree,range(4))))
#from inout import pointplot    
#pointplot(3,[bins,map(onetree,bins)])            
        
        
        
            






 



def expsmeer(profile,smeer_length):
    data = numpy.array(profile)
    if smeer_length <=1: return data

    trailing = numpy.zeros(smeer_length * 4,float)+numpy.mean(data[:smeer_length])
    finishing = numpy.zeros(smeer_length * 4,float)+numpy.mean(data[-smeer_length:])
    ndata = numpy.concatenate([trailing,data,finishing])
    

    data = ndata
    
    L = len(data)
    a = numpy.array(range(4*smeer_length+1),float) - 2*smeer_length
    add_profile = numpy.exp(-a*a/smeer_length**2)
    smeered = numpy.zeros(len(data),float)
    integral = numpy.sum(add_profile)
    for i in xrange(2*smeer_length+1,L-2*smeer_length):
        smeered[i-2*smeer_length:i+2*smeer_length+1] += data[i]*add_profile
        
    ndata =  smeered[4*smeer_length:-4*smeer_length]
    return ndata/integral
  






        



def test_lambda(l,ngen):
    
    weights = []
    while len(weights) < 3000000:
        while True:
            x = node()
            x.createLambdaChildren(l,ngen)
            if len(x.children) != 0:
                break
        allElements = x.giveAllChildren()
        for element in allElements:
            #element = random.choice(all)
            weights.append(len(element.giveAllChildren()))
    weights = na(weights)
    bins0 = logbins(2,max(weights),1.4)
    bins = [(bins0[i],bins0[i+1]) for i in xrange(len(bins0) - 1)]
    counts = [numpy.sum((i[0] < weights) * (weights < i[1])) for i in bins]
    possible = [i[1] - i[0] for i in bins]
            
    return [[(i[0]+i[1])/2. for i in bins],na(counts)/na(possible)]

def cross(a, b):
    return numpy.array((a[1] * b[2] - a[2] * b[1], a[2] * b[0] - b[2] * a[0], a[0] * b[1] - b[0] * a[1]))
def dot(a, b):
    return numpy.sum(a * b)
def SameSide(p1, p2, a, b):
    cp1 = cross(b - a, p1 - a)
    cp2 = cross(b - a, p2 - a)
    if dot(cp1, cp2) >= 0: return True
    else: return False
def PointInTriangle(p, a, b, c):
    if (SameSide(p, a, b, c) and SameSide(p, b, a, c) and SameSide(p, c, a, b)): return True
    else: return False

def intersect(t1, t2, t3, r1, r2):
    "intersection of a line and a triangle"
    t1 = numpy.array(t1)
    t2 = numpy.array(t2)
    t3 = numpy.array(t3)
    r1 = numpy.array(r1)
    r2 = numpy.array(r2)
    t2 -= t1
    t3 -= t1
    r1 -= t1
    r2 -= t1
    n = cross(t2, t3)
    c1 = dot(r1, n)
    c2 = dot(r2, n)
    if c1 * c2 > 0: return 0
    alpha = c1 / (c1 - c2)
    a = r1 + ((r2 - r1) * alpha)
    if not PointInTriangle(a, numpy.array((0, 0, 0)), t2, t3): return 0
    elif c1 > 0: return 1
    else: return - 1
#print cross(numpy.array((0,1,0)),numpy.array((0,0,1)))
#print intersect((1,1,1),(2,2,2),(1000,0,0),(10,2,1),(20,3,0))
def newintersect(t1, t2, t3, r1, r2):
    
    
    A = numpy.transpose(numpy.array([r1 - r2, t2 - t1, t3 - t1]))        
    D = numpy.array(r1 - t1)
    B = numpy.linalg.inv(A)
    C = numpy.dot(B, D)
    if (0 < C[0] < 1 and 0 < C[1] < 1 and 0 < C[2] < 1 and 0 < C[2] + C[1] < 1): return True
    else: return False
#print newintersect(numpy.array((1,1,1)),numpy.array((2,2,2)),numpy.array((1000,0,0)),numpy.array((10,2,1)),numpy.array((10,1,2)))
#raw_input()

def superintersect(t1, t2, t3, r1, r2):
    B = (t2[0] - t1[0], t2[1] - t1[1], t2[2] - t1[2])                         
    C = (t3[0] - t1[0], t3[1] - t1[1], t3[2] - t1[2])
    D = (r2[0] - t1[0], r2[1] - t1[1], r2[2] - t1[2])
    A = (r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2])
    det = A[0]*(B[1]*C[2]-C[1]*B[2])-A[1]*(B[0]*C[2]-B[2]*C[0])+A[2]*(B[0]*C[1]-C[0]*B[1])
    if det==0: return 0
    if det>0:
        t = D[0]*(B[1]*C[2]-C[1]*B[2])-D[1]*(B[0]*C[2]-B[2]*C[0])+D[2]*(B[0]*C[1]-C[0]*B[1])
        if t<0 or t>det: return 0
        u = A[0]*(D[1]*C[2]-C[1]*D[2])-A[1]*(D[0]*C[2]-D[2]*C[0])+A[2]*(D[0]*C[1]-C[0]*D[1])
        if u<0 or u>det: return 0
        v = A[0]*(B[1]*D[2]-D[1]*B[2])-A[1]*(B[0]*D[2]-B[2]*D[0])+A[2]*(B[0]*D[1]-D[0]*B[1])
        if v<0 or v>det: return 0
        if u+v>det: return 0        
        return 1
    if det<0:
        t = D[0]*(B[1]*C[2]-C[1]*B[2])-D[1]*(B[0]*C[2]-B[2]*C[0])+D[2]*(B[0]*C[1]-C[0]*B[1])
        if t>0 or t<det: return 0
        u = A[0]*(D[1]*C[2]-C[1]*D[2])-A[1]*(D[0]*C[2]-D[2]*C[0])+A[2]*(D[0]*C[1]-C[0]*D[1])
        if u>0 or u<det: return 0
        v = A[0]*(B[1]*D[2]-D[1]*B[2])-A[1]*(B[0]*D[2]-B[2]*D[0])+A[2]*(B[0]*D[1]-D[0]*B[1])
        if v>0 or v<det: return 0
        if u+v<det: return 0
        return 1
        
    
def createRW(N):
    data = numpy.zeros((N,3),float)
    for i in xrange(N-1):
        data[i+1,:] = data[i,:] + na(random_on_sphere(1))
    return data
def boundedRW(N):
    rad = (N**(1/3.) )
    print rad
    data = numpy.zeros((N,3),float)
    for i in xrange(N-1):
        data[i+1,:] = data[i,:] + na(random_on_sphere(1))
        while numpy.sqrt((data[i+1]**2).sum()) > rad:  
            data[i+1,:] = data[i,:] + na(random_on_sphere(1))
    return data
    
    

def writeRW():
    for size in [4000,32000,128000]:
        def dowrite(ind):
            t = boundedRW(size)
            a = open("/home/magus/evo/trajectories/boundedRWs/%d/rw%d.dat" % (size,ind),'w')
            a.write("%d\n" % size)
            for i in t:
                a.write("%f\t%f\t%f\n" % tuple(i))
            a.close()
        base.fmap(dowrite,range(200),n=4)

def get_distribution(x,num = 200):
    num  = float(num)
    x = numpy.array(x,float)
    import scipy.stats as st
    est = st.gaussian_kde(x)
    datamin,datamax = numpy.min(x),numpy.max(x)
    toevaluate = numpy.arange(datamin,datamax,(datamax - datamin)/num)
    return (toevaluate,est.evaluate(toevaluate))
                 
def getLinkingNumber(data1,data2):
    if len(data1) == 3: 
        data1 = numpy.array(data1.T)
    if len(data2) ==3: 
        data2 = numpy.array(data2.T)
    if len(data1[0]) !=3: raise ValueError
    if len(data2[0]) !=3: raise ValueError
    
    olddata = numpy.concatenate([data1,data2],axis = 0)
    olddata = numpy.array(olddata,dtype = float, order = "C")    
    M = len(data1)
    N = len(olddata)
    returnArray = numpy.array([0])
    
    support = r"""
#include <stdio.h>
#include <stdlib.h>

double *cross(double *v1, double *v2) {
    double *v1xv2 = new double[3];
    v1xv2[0]=-v1[2]*v2[1] + v1[1]*v2[2];
    v1xv2[1]=v1[2]*v2[0] - v1[0]*v2[2];
    v1xv2[2]=-v1[1]*v2[0] + v1[0]*v2[1];    
    return v1xv2;
}

double *linearCombo(double *v1, double *v2, double s1, double s2) {
    double *c = new double[3];
    c[0]=s1*v1[0]+s2*v2[0];
    c[1]=s1*v1[1]+s2*v2[1];
    c[2]=s1*v1[2]+s2*v2[2];
        return c;
}

int intersectValue(double *p1, double *v1, double *p2, double *v2) {
    int x=0;
    double *v2xp2 = cross(v2,p2), *v2xp1 = cross(v2,p1), *v2xv1 = cross(v2,v1);
    double *v1xp1 = cross(v1,p1), *v1xp2 = cross(v1,p2), *v1xv2 = cross(v1,v2);
    double t1 = (v2xp2[2]-v2xp1[2])/v2xv1[2];
    double t2 = (v1xp1[2]-v1xp2[2])/v1xv2[2];
    if(t1<0 || t1>1 || t2<0 || t2>1) {
        free(v2xp2);free(v2xp1);free(v2xv1);free(v1xp1);free(v1xp2);free(v1xv2);
        return 0;
    }
    else {
        if(v1xv2[2]>=0) x=1;
        else x=-1;
    }
    double *inter1 = linearCombo(p1,v1,1,t1), *inter2 = linearCombo(p2,v2,1,t2);
    double z1 = inter1[2];
    double z2 = inter2[2];
    
    free(v2xp2);free(v2xp1);free(v2xv1);free(v1xp1);free(v1xp2);free(v1xv2);free(inter1);free(inter2);
    if(z1>=z2) return x;
    else return -x;
}
    """
    
    
    code = r"""
    #line 1149 "numutils.py"
    double **data = new double*[N];
    int i,j;
    for(i=0;i<N;i++) {

        data[i] = new double[3];
        data[i][0]=olddata[3*i];
        data[i][1]=olddata[3*i+1];
        data[i][2]=olddata[3*i + 2];
    }
    
    int L = 0;
        for(i=0;i<M;i++) {
            for(j=M;j<N;j++) {
                double *v1, *v2;
                if(i<M-1) v1 = linearCombo(data[i+1],data[i],1,-1);
                else v1 = linearCombo(data[0],data[M-1],1,-1);
                
                if(j<N-1) v2 = linearCombo(data[j+1],data[j],1,-1);
                else v2 = linearCombo(data[M],data[N-1],1,-1);                
                L+=intersectValue(data[i],v1,data[j],v2);
                free(v1);free(v2);
            }
        }
        
    returnArray[0] =  L;

"""
    M,N #Eclipse warning removal
    scipy.weave.inline(code, ['M','olddata','N',"returnArray"], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
    return returnArray[0]
    

    

