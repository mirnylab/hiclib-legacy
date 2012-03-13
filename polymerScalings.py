import base
from base import fmap, fmapav 

import dnaUtils
import plotting 
import numutils
from numutils import logbins  

from dnaUtils import load, Cload, give_contacts 
import cPickle  
from math import exp,sqrt,log 

import numpy

import matplotlib.pyplot as plt 

from copy import copy  




def give_endpoint_scaling(base,shift=300, start = 100, end = 200,pace = 3):
    
    bins = numpy.arange(0.9,55,0.3)
    mul = 0.85,1    
    def dist(i, j,k):
        tx = datax[j:k] - datax[i]
        ty = datay[j:k] - datay[i]
        tz = dataz[j:k] - dataz[i]
        dists = numpy.sqrt(tx**2 + ty**2 + tz**2)
        return dists
    
    all_hists = numpy.zeros(len(bins)-1,float) 
    count = 0
    for i in xrange(start,end,pace):

            
        try:a = load(base % i)
        except:a = load(base)
        datax,datay,dataz = a[0],a[1],a[2]

      

                
        t = []
        for j in xrange(1,len(datax) - shift,shift/600+1):
            
                        
            t.append(dist(j,j+int(shift*mul[0]),j+int(shift*mul[1])))
            if len(t)*len(t[0]) > 300000:
                t = numpy.concatenate(t)
                all_hists += numpy.histogram(t,bins)[0]
                count +=len(t)
                t = []
     
            #print all_hists 
            #t = 1/(k*k*((k/R)**3 - 12*(k/R) + 16))
            #t = 1/k**2
             
            
            
    #print all_hists 
    myhist = all_hists
    #plot_occupation(all_hists)
    bins = numpy.array(bins)
    cbins = bins[:-1]/2+bins[1:]/2
    N = numpy.sum(myhist/(cbins**2))
    myhist  /= cbins**2
    
    
                  
            
    
    #print all_dists/count,data
    #radius of giration for current size
    
    
    rad = cbins[numpy.argmax(myhist)]
    b = [cbins,numpy.cumsum(myhist)]
    
     
    #rad = sum([b[1][i]*b[0][i] for i in xrange(len(b[0]))])/sum(b[1])  
    #b[1] = [b[1][i] / (b[0][i] * b[0][i] * ((b[0][i]/R)**3 - 12*(b[0][i]/R)+16))  for i in xrange(len(b[0]))]  #adjusting density    
    b[0] = [i/rad for i in b[0]] #rescaling data and bins
    
    b[1] = [(i*rad/N) for i in b[1]]
    
    return b
    

    #autopointprint("screening","screening",1)





def one_endpoint(shift):    
    #b = give_endpoint_scaling("/home/magus/evo/trajectories/globule_creation/128000_RW_expanded/crumpled%d.dat",shift)
    b = give_endpoint_scaling("/home/magus/evo/topo242_million/expanded%d.dat",shift)    
    #b = give_endpoint_scaling("/home/magus/evo/topo25_128k_saw_folding/data/crumpled%d.dat",shift)
    
    b = list(b)
    return b

def run_endpoint():
    shifts = [1000,3000,10000,30000,100000,300000]
    cr = fmap(one_endpoint,shifts,n=4)
       
    plotting.pointplot(1,*cr)
    plotting.pointplot(3,*cr)

            

#cProfile.run("run_endpoint()")            
            
            
        


def give_area_scaling(data):
    min_distance = 6
    Nrepeat = 10
    N = len(data[0])
    
    
    datax, datay, dataz = numpy.array(data[0]),numpy.array(data[1]),numpy.array(data[2])
    def dist(i,j):
        return sqrt((datax[i] - datax[j]) ** 2 + (datay[i] - datay[j]) ** 2 + (dataz[i] - dataz[j]) ** 2)
    def cendist(i): return sqrt(datax[i]**2+datay[i]**2+dataz[i]**2)
    def closest_point(point):
        return ((datax-point[0])**2 + (datay - point[1])**2+(dataz-point[2])**2).argmin()
    def choose_point(dist):
        while True:
            a = numpy.random.random(3)*1.5*dist
            if numpy.sum(a**2) < dist**2: return a
    bins = range(0,50,2)
    cutoff_distances = [2,4,6]
    ret = [[0.] * (len(bins)-1) for j in xrange(len(cutoff_distances))]
    for j in xrange(len(cutoff_distances)):
        cutoff = cutoff_distances[j]
        for _ in xrange(Nrepeat):             
            point = 0
            while True:
                a = closest_point(choose_point(min_distance))
                if 400 < a < N - 400: 
                    point
                    point = a
                    break
            for k in xrange(len(bins)-1):
                for f in xrange(bins[k],bins[k+1]):
                    if dist(a-f,a) < cutoff:
                        ret[j][k] += 1
                    if dist(a+f,a) < cutoff:
                        ret[j][k] += 1
    ret = [numpy.array(i) for i in ret]
    for i in xrange(len(ret)):
        for j in xrange(len(ret[0])):
            ret[i][j] /= float(2*Nrepeat * (bins[j+1] - bins[j]))
    return ret


def give_area_distribution(data):
    min_distance = 5.5
    Nrepeat = 5
    N = len(data[0])
    
    
    datax, datay, dataz = numpy.array(data[0]),numpy.array(data[1]),numpy.array(data[2])
    def dist(i,j):
        return sqrt((datax[i] - datax[j]) ** 2 + (datay[i] - datay[j]) ** 2 + (dataz[i] - dataz[j]) ** 2)
    def cendist(i): return sqrt(datax[i]**2+datay[i]**2+dataz[i]**2)
    def closest_point(point):
        return ((datax-point[0])**2 + (datay - point[1])**2+(dataz-point[2])**2).argmin()
    def choose_point(dist):
        while True:
            a = numpy.random.random(3)*1.5*dist
            if numpy.sum(a**2) < dist**2: return a
    bins = range(0,1000,10)
    cutoff_distances = [2.5,4]
    ret = [[0.] * (len(bins)-1) for j in xrange(len(cutoff_distances))]
    for j in xrange(len(cutoff_distances)):
        cutoff = cutoff_distances[j]
        for _ in xrange(Nrepeat): 
            point = 0
            while True:
                a = closest_point(choose_point(min_distance))
                if 1000 < a < N - 1000: 
                    point
                    point = a
                    break
            norm = 0
            for k in xrange(N):
                if dist(k,a) < cutoff:
                    norm += 1

            for k in xrange(len(bins)-1):
                for f in xrange(bins[k],bins[k+1]):
                    if dist(a-f,a) < cutoff:
                        ret[j][k] += 1/float(norm)
                    if dist(a+f,a) < cutoff:
                        ret[j][k] += 1./norm
            
                    
    ret = [numpy.array(i) for i in ret]
    for i in xrange(len(ret)):
        for j in xrange(len(ret[0])):
            ret[i][j] /= Nrepeat
    print "dist"
    #return [integrate(i) for i in ret]
        

                  

                    
                
        
            
#rets = [give_area_distribution(load("data//DNA_conf//cr//conforation_N4000_k1_%i_t200.dat" % i)) for i in xrange(20,60)]
#
#rets2 = [give_area_distribution(load("data//DNA_conf//uneq//new_protocol//conformation_N4000_k1_R_0_%i_t160.dat" % i)) for i in xrange(20,60)]
#
#plot_occupation(*retsum(rets))
#plot_occupation(*(retsum(rets)+retsum(rets2)) )
         
class drift():
    def __init__(self,filename):
        data = load(filename)
        self.datax,self.datay,self.dataz =  data[0],data[1],data[2]
        self.N = len(data[0])
        self.covered = numpy.zeros(self.N,int)
        self.maxR =  sqrt((self.datax**2+self.datay**2+self.dataz**2).max())
        self.stats = []
        self.k=16000
        self.curStep = 0 
        self.time = numpy.zeros(self.N,int)
    def closest_point(self,point):
        return ((self.datax-point[0])**2 + (self.datay - point[1])**2+(self.dataz-point[2])**2).argmin()
    def nearest_search(self,i,cutoff=1):
        x,y,z = self.datax[i],self.datay[i],self.dataz[i]
        b = numutils.random_on_sphere(cutoff)
        x+=b[0]
        y+=b[1]
        z+=b[2]
        return self.closest_point((x,y,z))
    def stat(self):
        return float(numpy.sum(self.covered))/len(self.covered)
    def step(self,i):
        self.curStep += 1
        shift = 70
        N = self.N
        covered = self.covered
        a = numpy.random.randint(1,3)
        if a==1:
            if i+shift <N:
                start = i
                stop = i + shift                
                ret = stop
            else:
                start = i 
                stop = N - 1
                ret = N - 1
            
        else:
            if (i-shift >= 0):
                start = i - shift
                stop = i 
                ret = i - shift                
            else:
                start = 0
                stop = i 
                ret = 0
        cur = covered[start:stop]
        mytime = self.time[start:stop]
        mytime[cur <0.5] = self.curStep 
        covered[start:stop] = 1
        return ret
                
                
    def doit(self,M):
        k = self.k
        for _ in xrange(M):
            k = self.nearest_search(k)
            #print k,
            k = self.step(k)
            #print k,
            #print self.stat()
            #self.stats.append(self.stat())
        self.k = k
        return self.time + self.curStep * (self.time <1)
        return numpy.array(self.covered,float)
    def give_stat(self):
        return numpy.array(self.stats,float)
    
def give_coverage_graph():    
    s1 = []
    s2 = []            
    for i in xrange(10,4000):
        if i%10 == 0: print i
        a = drift("data//DNA_conf//uneq//new_protocol//conformation_N4000_k1_R_0_%i_t160.dat" % (i%400+10))
        b = drift("data//DNA_conf//cr//conforation_N4000_k1_%i_t200.dat" % (i%400+10))
        a.doit(400)    
        s1.append(a.give_stat())
        b.doit(400)
        s2.append(b.give_stat())
    
    for i in xrange(1,len(s1)):
    
        s1[0] += s1[i]
        s2[0] += s2[i] 
    
    s1[0] /= float(len(s1))
    s2[0] /= float(len(s2))
    plotting.plot_occupation(s1[0],s2[0])
    
        

#give_coverage_graph()
          
def give_one_drift(filename):
    b = []
    for _ in xrange(30): 
        a = drift(filename)

        b.append(a.doit(250))
        
    b = numpy.array(b)
    
    return numpy.reshape(numpy.min(b,axis = 0),(1,-1))
def run_drifts():
    labels = []
    t = 0
    for _ in xrange(1):
        t += 3000 
        labels.append("time = " + str(t))
    
    plt.xlabel("Position")
    plt.ylabel("Time")
    a1 = [numutils.expsmeer(j,20) for j in base.fmapav(give_one_drift,["/home/magus/evo/trajectories/globule_creation/32000_SAW/crumpled%d.dat" % ((i % 35)+1) for i in xrange(500)])]
    a2 = [numutils.expsmeer(j,20) for j in base.fmapav(give_one_drift,["/home/magus/evo/trajectories/globule_creation/32000_equilibrium/equilibrium%d.dat" % ((i % 35)+1) for i in xrange(500)])]
    plotting.plot_occupation([a1[0],a2[0]],["fractal","equilibrium"])

        
#some_code()            
         
#exitProgram()            
def give_proximity_scaling(data):
        N = len(data[0])
        bins = range(0, N + 50, 50)
        radii = [0. for i in xrange(len(bins))]
        def dist(n):
            return sqrt(data[0][n] ** 2 + data[1][n] ** 2 + data[2][n] ** 2)
        for i in xrange(N):
            radii[i / 50] += dist(i) / N
        bins, radii = [exp(i) for i in range(len(bins))], [exp(i) for i in radii]
        return (bins, radii)
def give_density_scaling(data):
    mybins = numpy.arange(0, 20, 0.8)
    bins = [(mybins[i] + mybins[i - 1]) / 2 for i in xrange(1, len(mybins))]
    myres = [0. for i in xrange(len(mybins) - 1)]
    for i in xrange(len(data[0])):
        a = data[0][i] ** 2 + data[1][i] ** 2 + data[2][i] ** 2
        a = sqrt(a)
        for j in xrange(1, len(mybins)):
            if mybins[j] > a and mybins[j - 1] < a:
                
                myres[j - 1] += 1
                break
    for j in xrange(1, len(mybins)):
        myres[j - 1] /= (mybins[j] ** 3 - mybins[j - 1] ** 3)
    
    return [[exp(i) for i in bins], [exp(i) for i in myres]]
def new_external_contacts(filename, bins = None):
    data = load(filename)
    
    N = len(data[0])
    if bins == None: bins = logbins(4,N-4,1.1)
    b = give_contacts(data,cutoff=1.7)
    b = numpy.array(b)
    b = b [abs(b[:,1] - b[:,0])> 2.5]
    #print len(b), numpy.sum(b[:,0] == 0), numpy.sum(b[:,1] == 299),numpy.sum(b[:,0] == 1), numpy.sum(b[:,1] == 298)
    #print numpy.transpose(b)
    Ncontacts = len(b)
    begs = b[:,0]
    ends = b[:,1]
    lenses = ends - begs
    probs = []
    for M in bins:
        mask = lenses < M 
        mybegs = begs[mask]
        myends = ends[mask]
        mylenses = lenses[mask]
        myN = len(mylenses)
        Nrepeats = numpy.zeros(myN,int) + M - mylenses
        #print mybegs
        #print myends
        #print Nrepeats
        eatenlow = numpy.maximum(0,M - myends-1) 
        eatenhigh = numpy.maximum(0,M+mybegs - N )
        realrepeats = Nrepeats - eatenlow - eatenhigh
        localN = realrepeats.sum()
        
        Nrepeat = numpy.zeros(Ncontacts,int) + 2*M
        eatenbetween = numpy.maximum(0,M-lenses)
        eatenlow = numpy.maximum(0,M-begs - 1)
        eatenhigh = numpy.maximum(0,M+ends - N)
        realrepeats = Nrepeat - eatenbetween - eatenlow - eatenhigh
        globalN = realrepeats.sum()
        print localN / float(N-M+1), globalN/float(N-M+1) 
        probs.append((2 * localN) / float(localN + globalN))
    
    return  numpy.array([bins,probs],float)
    
        #print realrepeats
def run_external_contacts():    
    #a = fmapav(new_external_contacts,["/home/magus/evo/trajectories/globule_creation/32000_SAW_expanded/crumpled%d.dat" % i for i in xrange(1,40)])
    #a = fmapav(new_external_contacts,["/home/magus/evo/trajectories/RWs/128000/rw%d.dat" % i for i in xrange(1,16)])
    a = fmapav(new_external_contacts,["/home/magus/evo/trajectories/globule_creation/32000_equilibrium/equilibrium%d.dat" % i for i in xrange(1,16)])
    
    #a = new_external_contacts("/home/magus/evo/trajectories/hilbert.dat")
    plt.xscale("log")
    plt.plot(*a)
    plt.show()
    exit()    
    
    

#N = 32000
#pointplot(1,[bins,eq],[bins,fr],[bins,[2*float(i)**-0.3333 for i in bins]],[bins,[(log(N)-log(i))/(log(N)-1) for i in bins]],[bins,[2*((N/float(i))**0.33333333 - 1 )/(N**0.3333333 - 1) for i in bins]])
#autopointprint("external_cont","external_cont",1
def give_dist_scaling_ext(data, bins0=None, cutoff=1.1,integrate = False,ring=False):
    "gives Pc scaling for data for rings only using external program"
    lines = give_contacts(data, cutoff)
    lines = lines[:,1] - lines[:,0]
    
    
    N = len(data[0])
 
    lines[lines>N/2] = N - lines[lines>N/2]
    
    lines = numpy.bincount(lines)
    lines.resize(N/2)
    #plot_occupation(lines)
    #plot_occupation(lines) 
    globconnections = 0      
    #print data,t,lines
    
    #print lines[:140]
    bins = [(bins0[i], bins0[i + 1]) for i in xrange(len(bins0) - 1)]
    connections = [1 for i in xrange(len(bins))]
    counts = [0 for i in xrange(len(bins))]    
    count = 0
    
    #print bins
    for i in xrange(len(bins)):
        count = 0
        
        oneBin = bins[i] 
        for j in xrange(oneBin[0], oneBin[1],1):
            
         
            
              
            
            #yo = numpy.sum(numpy.sum((data[:,0:N - j] - data[:,j:N])**2 ,0) < cutoff*cutoff)
            #yo = numpy.sum(numpy.sum((data[0:N - j,:] - data[j:N,:])**2 ,1) < cutoff*cutoff)
            count += lines[j]        
                        
            #for k in xrange(0,N-j):             
            #    if dist(k,k+j) < cutoff:
            #        count += 1
            counts[i] += N 
                        
        connections[i] = count
        globconnections += float(count)
        if connections[i] == 0: connections[i] = 0.001 
    print "average connections per unit:", float(globconnections) / N
    if integrate == False: 
        for i in xrange(len(connections)): connections[i] = connections[i] / float(counts[i])
    if integrate == True:
        print "integrate!"
        for i in xrange(1,len(connections)):
            connections[i] = connections[i] + connections[i-1]
        
        connections = [i/float(connections[-1])  for i in connections]         
    
    a = [sqrt(i[0] * i[1]) for i in bins]
    print connections
    b = [i for i in connections]
    #pointplot(1,[a,b])  
    
    #print datax
    #plot_occupation(numpy.array(datax),numpy.array(datay))
    return [a, b]
def give_dist_scaling_old(data, bins0=None, cutoff=1.1,integrate = False,ring=False):
    #t = int(numpy.sum(data[0][50:500]))
    #print t
    globconnections = 0         
    N = len(data[0])
    if bins0 == None:  bins0 = logbins(4,3998,1.25)      
#    print lines
    bins = [(bins0[i], bins0[i + 1]) for i in xrange(len(bins0) - 1)]
    connections = [0 for i in xrange(len(bins))]
    counts = [0 for i in xrange(len(bins))]    
    count = 0
    
    datax,datay,dataz = data[0],data[1],data[2]
    for i in xrange(len(bins)):
        count = 0
        
        oneBin = bins[i] 
        for j in xrange(oneBin[0], oneBin[1],1):
            
         
            
             
            
            #yo = numpy.sum(numpy.sum((data[:,0:N - j] - data[:,j:N])**2 ,0) < cutoff*cutoff)
            #yo = numpy.sum(numpy.sum((data[0:N - j,:] - data[j:N,:])**2 ,1) < cutoff*cutoff)
                         
            yo = numpy.sum((datax[0:N - j] - datax[j:N])**2 + (datay[0:N - j] - datay[j:N])**2+ (dataz[0:N - j] - dataz[j:N])**2  < cutoff*cutoff)
            count += yo
            #for k in xrange(0,N-j):             
            #    if dist(k,k+j) < cutoff:
            #        count += 1
            counts[i] += N-j            
        connections[i] = count
        globconnections += float(count)
        if connections[i] == 0: connections[i] = 0.1 
    print "average connections per unit:", float(globconnections) / N
    if integrate == False: 
        for i in xrange(len(connections)): connections[i] = connections[i] / float(counts[i])
    if integrate == True:
        print "integrate!"
        for i in xrange(1,len(connections)):
            connections[i] = connections[i] + connections[i-1]
        
        connections = [i/float(connections[-1])  for i in connections]         
    
    a = [sqrt(i[0] * i[1]) for i in bins]
    print connections
    b = [i for i in connections]
    #pointplot(1,[a,b])  
    
    #print datax
    #plot_occupation(numpy.array(datax),numpy.array(datay))
    return [a, b]
def give_dist_scaling(data, bins0=None, cutoff=1.1,integrate = False,ring=False,project=False,intContacts = False):
        #t = int(numpy.sum(data[0][50:500]))
    #print t
    globconnections = 0         
    N = len(data[0])          
    bins = [(bins0[i], bins0[i + 1]) for i in xrange(len(bins0) - 1)]
    connections = [0 for i in xrange(len(bins))]
    counts = [0 for i in xrange(len(bins))]    
    count = 0
    if intContacts == False:
        if (project == False): contacts = numpy.array(give_contacts(data,cutoff))
        else: contacts = numpy.array(dnaUtils.contacts_project(data,cutoff))
    else:
        contacts = dnaUtils.giveIntContacts(data)
    contacts = contacts[:,1] - contacts[:,0]
    for i in xrange(len(bins)):
        count = 0
        
        oneBin = bins[i] 
        for j in xrange(oneBin[0], oneBin[1],1):
            
            counts[i] += N-j
        count = numpy.sum((contacts>=oneBin[0])*(contacts<oneBin[1]))            
        connections[i] = count
        globconnections += float(count)
        if connections[i] == 0: connections[i] = 0.0000000000001 
    print "average connections per unit:", float(globconnections) / N
    if integrate == False: 
        for i in xrange(len(connections)): connections[i] = connections[i] / float(counts[i])
    if integrate == True:
        print "integrate!"
        for i in xrange(1,len(connections)):
            connections[i] = connections[i] + connections[i-1]
        
        connections = [i/float(connections[-1])  for i in connections]         
    
    a = [sqrt(i[0] * i[1]) for i in bins]
    print connections
    b = [i for i in connections]
    #pointplot(1,[a,b])  
    
    #print datax
    #plot_occupation(numpy.array(datax),numpy.array(datay))
    return [a, b]




#print give_dist_scaling_ext(load("/home/magus/evo/topo35_108_battery/run1/expanded1.dat",False)
def new_dist_scaling(data, bins0, cutoff=1.1):
    bins = bins0
    dists = numpy.sqrt(numpy.sum(data**2,0))
    touse = dists >15
    beg = numpy.nonzero((touse[:-1]==0) * (touse[1:] == 1))[0]+1
    end = numpy.nonzero((touse[:-1] == 1) * (touse[1:] == 0))[0]
    if beg[0] > end[0]: end = end[1:]
    areas =  zip(beg,end)
    areas = filter(lambda i:i[1]-i[0]>10,areas)
    connections = [0] * len(bins) 
    counts = [0 for i in xrange(len(bins))]    

     
    def subchain_counts(data,bins,cutoff):
        
        connections,counts       
        N = len(data[0])
        bins = [(bins0[i], bins0[i + 1]) for i in xrange(len(bins0) - 1)]
        datax,datay,dataz = data[0],data[1],data[2]
        for i in xrange(len(bins)):
            count = 0
            
            oneBin = bins[i] 
            for j in xrange(oneBin[0], oneBin[1],1):
                if j >=N:
                    break
        
         
            
             
            
            #yo = numpy.sum(numpy.sum((data[:,0:N - j] - data[:,j:N])**2 ,0) < cutoff*cutoff)
            #yo = numpy.sum(numpy.sum((data[0:N - j,:] - data[j:N,:])**2 ,1) < cutoff*cutoff)
                         
                yo = numpy.sum((datax[0:N - j] - datax[j:N])**2 + (datay[0:N - j] - datay[j:N])**2+ (dataz[0:N - j] - dataz[j:N])**2  < cutoff*cutoff)
                count += yo
                #for k in xrange(0,N-j):             
                #    if dist(k,k+j) < cutoff:
                #        count += 1
                counts[i] += N-j            
            connections[i] += count
    
            
    for i in areas: subchain_counts(data[:,i[0]:i[1]],bins0,cutoff)
    return connections,counts
def scale_new_distance():
    mybins = logbins(4,8000,1.3)
    bins0 = mybins
    bins = [(bins0[i], bins0[i + 1]) for i in xrange(len(bins0) - 1)]
    b = [sqrt(i[0] * i[1]) for i in bins]
    def myfun(x):
        return new_dist_scaling(load(x),mybins,1.7)
    
    filenames = ["/home/magus/evo/topo202_scale_fold_32000_55000/run%d/expanded%d.dat" % (i,j)  for j in xrange(550,700,20) for i in xrange(1,41)] 
    #["/home/magus/evo/topo202_scale_fold_32000_65000/run%d/expanded%d.dat" % (i,j)  for j in xrange(550,700,20) for i in xrange(1,41)]
    
    
    
    counts = zip(*fmap(myfun,filenames))
    
    connections = numpy.sum(counts[0],0)
    counts = numpy.sum(counts[1],0)
    #print counts,connections
    m = max([i for i  in xrange(len(connections)) if connections[i] > 0]) + 1
    counts = counts[:m]
    connections = connections[:m]
    b = b[:m]
    p = plotting.plot([b,connections/numpy.array(counts,float)],"point",label="external>13")
    cPickle.dump([[p],[p],[p]],open("data/DNA_conf/plots/32_external",'wb'))
    
    
    
#scale_new_distance()
#exitProgram()
def give_saw_scaling(data, bins0=None, cutoff=6,integrate = False):
    datax, datay, dataz = data[0], data[1], data[2]
    N = len(data[0])
    if bins0 == None:  bins0 = logbins(4,3998,1.25)  
    bins = [(bins0[i], bins0[i + 1]) for i in xrange(len(bins0) - 1)]
    connections = [1 for i in xrange(len(bins))]
    counts = [0 for i in xrange(len(bins))]    

    for i in xrange(len(bins)):
        oneBin = bins[i] 
        for j in xrange(oneBin[0], oneBin[1], 1):
            counts[i] += N-j
            
    
    for i in xrange(0,N):
        dists = (datax - datax[i])**2 + (datay - datay[i])**2 + (dataz - dataz[i])**2
        mons = numpy.abs((numpy.argsort(dists)[1:cutoff+1]) - i)
        #print mons         
        for k in mons:
            for j in xrange(len(bins)):
                oneBin = bins[j]
                if ((oneBin[1] > k) and (oneBin[0] <= k)):
                    connections[j] +=1;
                    #print k,j
    print connections                      
    for i in xrange(len(connections)): connections[i] = connections[i] / float(counts[i])
    if integrate == True:
        print "integrate!"          
        connections[0]*= (bins[0][1]-bins[0][0])
        for i in xrange(1,len(connections)): connections [i] = connections[i-1] + connections[i]*(bins[i][1]-bins[i][0])

    
    a = [sqrt(i[0] * i[1]) for i in bins]
    b = [i for i in connections]
    #pointplot(1,[a,b])  
    
    print connections
    #print datax
    #plot_occupation(numpy.array(datax),numpy.array(datay))
    return [a, b]


def give_fractality(data):
    N = len(data[0])
    data = numpy.array(data)
    gdatax = numpy.array(data[0])
    gdatay = numpy.array(data[1])
    gdataz = numpy.array(data[2])

    def radius_gyration(data, i, j):
        datax = gdatax[i:j]
        datay = gdatay[i:j]
        dataz = gdataz[i:j]
        avx, avy, avz = numpy.sum(datax), numpy.sum(datay), numpy.sum(dataz)
        N = len(datax)
        avx, avy, avz = avx / N, avy / N, avz / N        
        dists = numpy.sqrt(((datax - avx) ** 2) + ((datay - avy) ** 2) + ((dataz - avz) ** 2))
        return numpy.sum(dists) / N
    

    R1 = radius_gyration(data,0,N-1)
    n = int(N**(0.6666666))

    k = 0
    r2 = 0
    for i in xrange(10,N-10-n,n/10+1):
        k+=1
        r2 += radius_gyration(data,i,i+n)
    return (R1/(r2/k))


def give_distance(data, bins=None,ring = False ):
    "main working horse for end-to-end distance"
    N = len(data[0])
    if ring == True:
        data = numpy.concatenate([data,data],axis = 1)
    bins = [(bins[i], bins[i + 1]) for i in xrange(len(bins) - 1)]
    
    rads = [0. for i in xrange(len(bins))]  
    for i in xrange(len(bins)):
        oneBin = bins[i] 
        rad = 0.
        count = 0
        for j in xrange(oneBin[0],oneBin[1], (oneBin[1]-oneBin[0])/10 + 1):
            length = j 
            if ring == True: 
                rad += numpy.mean(numpy.sqrt(numpy.sum((data[:,:N]-data[:,length:length+N])**2,0)))
            else: 
                rad += numpy.mean(numpy.sqrt(numpy.sum((data[:,:-length]-data[:,length:])**2,0))) 
            count += 1
         
        rads[i] = rad/count
    bins = [sqrt(i[0] * i[1]) for i in bins] 
    return (bins, rads) 



def give_radius_scaling(data, bins=None,ring = False):
    "main working horse for radius of gyration"
    "uses dymanic programming algorithm"
    
    data = numpy.array(data,float)
    coms = numpy.cumsum(data,1)   #cumulative sum of locations to calculate COM
    coms2 = numpy.cumsum(data**2,1)  #cumulative sum of locations^2 to calculate RG    
    def radius_gyration(len2):
        data 
        if ring == True:
            comsadd = coms[:,:len2].copy()
            coms2add = coms2[:,:len2].copy()
            comsadd += coms[:,-1][:,None]
            coms2add += coms2[:,-1][:,None]
            comsw = numpy.concatenate([coms,comsadd],axis = 1)
            coms2w = numpy.concatenate([coms2,coms2add],axis = 1)  #for rings we extend need longer chain
        else:
            comsw = coms
            coms2w = coms2
            
        coms2d = (-coms2w[:,:-len2]+coms2w[:,len2:])/len2
        comsd = ((comsw[:,:-len2]-comsw[:,len2:])/len2)**2
        diffs = coms2d - comsd
        sums = numpy.sqrt(numpy.sum(diffs,0))
        return numpy.mean(sums)
        
    rads = [0. for i in xrange(len(bins))]
    for i in xrange(len(bins)):
        rads[i] = radius_gyration(bins[i])
    return (copy(bins), rads) 



def give_radius_scaling_eig(data, bins=None):
    #gives volume^^0.33 as  defined through  eigenvectors
    if bins==None: bins = [2*i for i in logbins(1,0.45*len(data[0]),1.3,40)]
    x,y,z = data[0],data[1],data[2]
    coords = [x,y,z]
    sums = [[i*j for j in coords] for i in coords]
    sums = numpy.array(sums) 
    N = len(data[0])
    coms = numpy.cumsum(data,1)
    sums = numpy.cumsum(sums,2)
    def tensor(a,b):        
        newsums = (sums[:,:,b]-sums[:,:,a])/float(b-a)
        newcoms = (coms[:,b] - coms[:,a])/float(b-a)
        tensor =  newsums - numpy.outer(newcoms,newcoms)
        return numpy.linalg.eigvalsh(tensor)
    ret = []
    for i in bins:
        av = 0.
        for j in xrange(1000):
            t = numpy.random.randint(5,N-5-i)
            res = tensor(t,t+i)
            av += sqrt(3) * (res[0] * res[1] * res[2] * 1. )**(1/6.) 
        ret.append(av/1000)
    retret = (copy(bins),ret)
    return  retret


def give_slices(base, tosave, slices, sliceParams, multipliers, mode = "chain", loadFunction = Cload, integrate = False,normalize=False):
    numpy.seterr(invalid='raise')
    sliceplots1 = []
    sliceplots2 = []
    sliceplots3 = []
    
       



    

    for cur_slice in slices: 
        files = []

        def slice2D(a, b,mult=[1]):             
            tm = []
            if type(b) == tuple:
                for i in xrange(b[0], b[1]+1):
                    tm.append((i, a))
            elif type(b) == int:
                for i in xrange(1, b + 1):
                    tm.append((i, a))
            elif type(b) == list: tm=[(i,a) for i in b]
            tm2 = sorted(list(set([(i[0],int(float(i[1])*m)) for i in tm for m in mult])))
            print tm2
            return tm2  
        
        def slice3D(a, b,c,mult=[1]):
            tm = []
            for i in xrange(b[0], b[1]+1):
                for t in xrange(c[0],c[1]+1):
                    tm.append((i, a,t))
            tm2 = sorted(list(set([(i[0],int(float(i[1])*m)) for i in tm for m in mult])))
            print tm2
            return tm2   
        
        #sluces actually are defined
        runs = slice2D(cur_slice, sliceParams ,multipliers)
        #runs = slice3D(cur_slice, (1,14),(1,10),multipliers)
        
        for i in runs:            
            #filename is replaced in slices
            try: files.append(base.replace("DATA1", str(i[0])).replace("DATA2", str(i[1])).replace("DATA3",str(i[2])))
            except: files.append(base.replace("DATA1", str(i[0])).replace("DATA2", str(i[1])))
            
    
        datas = []
        MM = len(datas)
        
        plots1 = [0 for i in xrange(MM)]
        plots2 = [0 for i in xrange(MM)]
        plots3 = [0 for i in xrange(MM)]

        def newload(i):
            #loads a file  
            try:                    
                data =  loadFunction(i,False)
                return data
            except: 
                print "file not found" , i
                return None
            
        
        #use this for determining the file size 
        datas = filter(lambda x:x!=None,fmap(newload,files[::len(files)/20 + 1],n=3))
        datlen = len(datas[0][0])
        if mode == "chain": bins2 =  logbins(4, datlen-100,1.25)        
        if mode == "parts": bins2 =  logbins(4, datlen-100,1.25)
        if (mode == "ring") or (mode == "intring"): 
            b1 = logbins(2, datlen/4-1,1.34)
            bins2 =  [2*i for i in b1]
            print bins2
        binsrg = logbins(4, datlen-100,1.25)
        
        
         
        def give_plots(i):
            data = newload(i)
            if data == None: return None                  
            i = data                  

            if (mode == "ring") or (mode == "intring"):
                b = give_radius_scaling(i,binsrg,ring = True)
            else: 
                b = give_radius_scaling(i,binsrg,ring = False)
            
            
            if (mode == "chain") : a = give_dist_scaling(i, bins2,1.7,integrate,project=False)
            if (mode == "ring"): a = give_dist_scaling(i, bins2,1.7,integrate,ring = True)
            if (mode == "intring"): a = give_dist_scaling(i, bins2,1.7,integrate,ring = True,project=False,intContacts = True) 
            if (mode  == 'parts'): a =new_dist_scaling(i, bins2,1.7)
            if (mode == "project"): a = give_dist_scaling(i, bins2,1.450,integrate,project=True)
            
            if (mode == "ring") or (mode == "intring"):
                c = give_distance(i, bins2,ring = True)
            else: 
                c = give_distance(i, bins2,ring = False)
                
            
                     
            flatten = normalize
            if (flatten == True):
            #some weired corrections... 
                print "using normal correction!!!"
                a[1] = [a[0][i]*a[1][i] for i in xrange(len(a[0]))]
                b = list(b)
                c = list(c)
                b[1] = [(b[1][i]/(b[0][i]**0.333333)) for i in xrange(len(b[0]))]
                c[1] = [(c[1][i]/(c[0][i]**0.333333)) for i in xrange(len(c[0]))]
            elif flatten=="super":
                print "using super correction!!!"
                b = list(b)
                c = list(c)
                for i in b[0]:
                    e = (log(datlen) - log(b[0][i]))/(log(datlen) - log(1))
                    #print e
                    b[1][i] = b[1][i]  * ((2-2*e)/(2-e))**(1/3.)
                    e = (log(datlen) - log(c[0][i]))/(log(datlen) - log(1))
                    c[1][i] = c[1][i] * ((2-2*e)/(2-e))**(1/3.)                                     
            return [a,b,c]
                
        
        parPlots = fmap(give_plots,files,n=7)
        parPlots = filter(lambda x:x!=None,parPlots)
        
        for i in xrange(len(parPlots)):
            p1 = parPlots[i][0]
            p2 = parPlots[i][1]
            p3 = parPlots[i][2]
            plots3.append(p3[1])
            plots2.append(p2[1]) 
            plots1.append(p1[1])
            b3 = p3[0]
            b2 = p2[0]
            b1 = p1[0]
            
        
        plot1 = plotting.plot([b1,plots1],"sigmaplot",label="t="+str(cur_slice))
        plot1.base = base
        plot1.multipliers = multipliers
        plot1.mode = mode
        plot1.runs = runs
        plot2 = plotting.plot([b2,plots2],"sigmaplot",label="t="+str(cur_slice))
        plot3 = plotting.plot([b3,plots3],"sigmaplot",label="t="+str(cur_slice))
        sliceplots1.append(plot1)
        sliceplots2.append(plot2)
        sliceplots3.append(plot3)        
    
        
    cPickle.dump([sliceplots1,sliceplots2,sliceplots3],open(tosave,'wb'),-1)
    print "Finished!!!"
    return sliceplots1,sliceplots2,sliceplots3


#give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run4eq/expandedDATA2.dat", 
#            tosave  = "data/DNA_conf/plots/paper_scalings/ring13eq", 
#             slices = [11000], sliceParams = (1,1), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)

#give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/runDATA1/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/ring13", 
#             slices = [300,1000,70000], sliceParams = (3,4), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)


#give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run5eq/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/chain32eq", 
#             slices = [7000], sliceParams = (1,1), multipliers = numpy.arange(0.5,1,0.0001), mode = "chain", loadFunction = Cload)

#give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run5/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/chain32", 
#             slices = [300,1000,30000], sliceParams = (1,1), multipliers = numpy.arange(0.5,1,0.0001), mode = "chain", loadFunction = Cload)
#  
#give_slices(base = "/home/magus/evo/topo371_256_equilibrium/runDATA1/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/ring_256_eq", 
#             slices = [2,3,7,13], sliceParams = (1,10), multipliers = numpy.arange(0.7,1,0.0001), mode = "intring", loadFunction = intload)

#-------- equilibrium Rings low density -----------

#give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run6_lowden/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/ring13_lowd", 
#             slices = [6200], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)


#give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run6_eq/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/ring13_lowd_eq", 
#             slices = [4600], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)

#---------------ETE distance for N=108000 ----------------
#give_slices(base = "/home/magus/evo/topo36_108_grow/runDATA1/expandedDATA2.dat",
#             tosave  = "data/DNA_conf/plots/paper_scalings/ring_108_eq", 
#             slices = [20000], sliceParams = (1,10), multipliers = numpy.arange(0.7,1,0.01), mode = "intring", loadFunction = intload)

#---------------------------Equilibrium rings, 4k ---------
give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run7_small/expandedDATA2.dat",
             tosave  = "/home/magus/workspace/testnucl2/data/DNA_conf/plots/paper_scalings/ring4", 
             slices = [4000], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)
exit()

give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run7_eq/expandedDATA2.dat",
             tosave  = "data/DNA_conf/plots/paper_scalings/ring4_eq", 
             slices = [7000], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)

give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run8_tiny/expandedDATA2.dat",
             tosave  = "data/DNA_conf/plots/paper_scalings/ring2", 
             slices = [14000], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)


give_slices(base = "/home/magus/evo/GO37_6k_diffusion/equilibration_new/run8_tiny_eq/expandedDATA2.dat",
             tosave  = "data/DNA_conf/plots/paper_scalings/ring2_eq", 
             slices = [9000], sliceParams = (3), multipliers = numpy.arange(0.5,1,0.0001), mode = "ring", loadFunction = Cload)

