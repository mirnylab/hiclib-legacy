import base 
base #Eclipse warning damper 
import numpy 

from copy import copy 
from math import sqrt 



#from inout import *

class Statistics:
    def __init__(self,data):
        #print data
        self.N = len(data)
        N = self.N
        self.data = [numpy.array(i,float) for i in data]
        self.lengthes = []
        for i in xrange(N):
            self.lengthes.append(len(data[i]))
        self.SortedData = []
        for i in xrange(N):
            self.SortedData.append(copy(self.data[i])) 
        for i in xrange(N):
            self.SortedData[i].sort()
        self.expcutoff = 0.3
        
        
        self.RocN = 30
        #Number of ROC poitns
        
    def giveROC(self,a,b,n):
        "first is prediction, second is experiment"
        cutoff = float(n)/float(self.RocN)
        lena = self.lengthes[a]
        lenb = self.lengthes[b]
        lenmax = min(lena,lenb)
        data1 = self.data[a][:lenmax]
        data2 = self.data[b][:lenmax]
        med1 = self.SortedData[a][int(lenmax*cutoff)]
        med2 = self.SortedData[b][int(lenmax*self.expcutoff)]
        tp = numpy.sum((data1>med1)&(data2>med2))
        #both positive
        fp =   numpy.sum((data1>med1) & (data2<med2))
        tn = numpy.sum((data1<med1) & (data2<med2))
        fn = numpy.sum((data1<med1) & (data2>med2))
        #< changed to >
        #our is positive, real is negative)    
        return tuple([float(i) for i in [tp,fp,tn,fn]])
    
    def binROC(self,a,b,n):
        "first is prediction, second is binar experiment; 1 is match"
        cutoff = float(n)/float(self.RocN)
        lena = self.lengthes[a]
        lenb = self.lengthes[b]
        lenmax = min(lena,lenb)
        data1 = self.data[a][:lenmax]
        data2 = numpy.array(self.data[b][:lenmax],bool)
        data2n = numpy.invert(data2)
        med1 = self.SortedData[a][int(lenmax*cutoff)]
        #med2 = self.SortedData[b][int(lenmax*self.expcutoff)]
        tp = numpy.sum((data1>med1)& data2)
        #both positive
        fp =   numpy.sum((data1>med1) & data2n)
        tn = numpy.sum((data1<med1) & data2n)
        fn = numpy.sum((data1<med1) & data2)
        #< changed to >
        #our is positive, real is negative)    
        return tuple([float(i) for i in [tp,fp,tn,fn]])
        
    def plotROC(self,a,b):
        ROCx = numpy.array([0. for _ in xrange(self.RocN)],float)
        ROCy = numpy.array([0. for _ in xrange(self.RocN)],float)
        for i in xrange(1,self.RocN):
            myroc = self.binROC(a,b,i)
            #print myroc
            ROCx[i] = myroc[1]/(myroc[1]+myroc[2])
            ROCy[i] = myroc[0]/(myroc[0] + myroc[3])
        #print ROCx
        #print ROCy
        return (ROCx,ROCy)
    def area(self,a,b):
        self.ROCN = 20
        a = self.plotROC(a,b)
        x = list(copy(a[0]))
        y = list(copy(a[1]))
        y[0] = 1
        x[0] = 1
        y.append(0)
        x.append(0)
        mySum = 0.
        for i in xrange(len(x) -1):
            sum += (x[i]-x[i+1])*(y[i]+y[i+1])/2
        return mySum
            
        
        pass
        
def ROC_area(exp,data2):
    stat = Statistics([exp,data2])
    return stat.area(0,1)

def give_av_plot(data,method = "median"):
    def search(array,num):
        a = 0
        b = len(array) -1
        for _ in xrange(100):
            half = (a+b)/2            
            if array[half] > num:
                if array[half-2] < num:
                    return half-1
                b = half
            else:
                a = half
        return (a+b)/2
    
    def sigma(data,av):
        if len(data) == 0: return 0
        a = sqrt(numpy.sum((data - av)**2)/float(len(data)))
        
        
        return a
        
    minlen = 10
    per=5    
    N = len(data)    
    working = False
    cut = N
    for i in xrange(N):
        if len(data[i]) == 0: 
            data[i] = [0]
            if working == False:
                working = True
                cut = i
        else:
            working = False
             
        data[i] = numpy.array(data[i],float)
        data[i] += 0.000001*numpy.random.random(len(data[i]))
        data[i].sort()    
    N = cut
    print N
    if method == "average":
        av = [0 for i in xrange(N)]
        errtop = [0 for i in xrange(N)]
        errbot = [0 for i in xrange(N)]
        for i in xrange(N):            
            av[i] = numpy.sum(data[i])/len(data[i])
            mid = search(data[i],av[i])
            errtop[i] =  av[i] + sigma(data[i][mid:],av[i])
            errbot[i] =  av[i] - sigma(data[i][:mid],av[i])
        ret = (av,errbot,errtop)
    if method == "median1":
        points = [0.5,0.25,0.75]
        ret = [[0. for i in xrange(N)] for j in xrange(3)]
        for i in xrange(N):
            for j in xrange(3):
                ret[j][i] = data[i][int(len(data[i])*points[j])]
        ret
    if method == "median2":
        points = [0.125,0.25,0.5,0.75,0.875]
        ret = [[0. for i in xrange(N)] for j in xrange(5)]
        for i in xrange(N):
            for j in xrange(5):
                ret[j][i] = data[i][int(len(data[i])*points[j])]
    points = [[],[]]
    ret2 = [[] for i in xrange(len(ret)-1)]
    i = 0
    while i<N:
        if len(data[i]) >  minlen and i%per == 0:
       
            points[1].append(ret[0][i])
            points[0].append(i)
            for j in xrange(1,len(ret)):
                ret2[j-1].append(ret[j][i])
            i+=1
        else:            
            points[1].append(0)
            points[0].append(0)
            for j in xrange(1,len(ret)):
                ret2[j-1].append(0)
            
            count = 0            
            while ((count < minlen or i%per !=1)   and i<N): 
                count += len(data[i])
                l = len(data[i])
                points[1][-1] += ret[0][i] * l
                points[0][-1] += i*l
                for j in xrange(1,len(ret)):
                    ret2[j-1][-1] += ret[j][i] * l
                i+=1
            try: points[1][-1] = points[1][-1]/count
            except: print i
            points[0][-1] = points[0][-1]/count
            for j in xrange(len(ret2)):
                ret2[j][-1] = ret2[j][-1] / count
                
                
            
            
     
    return points,ret2
    
                
                
        
            
#a = give_av_plot([numpy.random.normal(1,1,5000),2*numpy.random.random(1000),[1],[1,2,3,4,5,6,7],[1],[2,3,4,5],[1],[2],[3],[]],"average")           
#error_script(a[0],a[1],"fuck me")
            
            
            
            
        
def point_distribution(a):
    a = numpy.array(a,int)
    low,high = numpy.min(a),numpy.max(a)
    ret = numpy.zeros(high-low+1,int)
    for i in a:
        ret[i-low] += 1
    return [range(low,high+1),ret]




#inout.pointplot(1,*filter_distribution(range(10000),range(10000),[100,3000,5000,8000],14))
            
         
            
         
#a = Statistics([[1,2,3,4,5,6,7,3,4,5,54,63,6,5],[7,6,5,4,3,2,1,6,3,4,5,54,63,6,5]])
#print a.plotROC(0,1)            