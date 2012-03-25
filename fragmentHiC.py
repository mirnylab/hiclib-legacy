import base
base #Eclipse warning remover 
import os,cPickle
 
  
import numpy
from numpy import array as na  
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt 

import plotting 
from plotting import mat_img,removeAxes
import numutils 
from numutils import arrayInArray, arraySumByArray, sumByArray, correct, ultracorrect
import dnaUtils
import numexpr
import joblib 

import h5py
r_ = numpy.r_

def corr(x,y): return stats.spearmanr(x, y)[0]        

def saveData(data,filename,override = True):
    "function to save numpy arrays on HDD"
    if override == True:
        if os.path.exists(filename):
            os.remove(filename)
    joblib.dump(data,filename,compress = 3)
    
"TODO: Rewrite glue between my and Anton's code; make universal interface for other types of data" 


class track(object):
    def __init__(self,folder , genomeFolder , maximumMoleculeLength = 500,override = True):
        #These are fields that will be kept on a hard drive
        self.vectors = {"chrms1":"int8","chrms2":"int8", #chromosomes. If >chromosomeCount, then it's second chromosome arm! 
                        "mids1":"int64","mids2":"int64",  #midpoint of a fragment, determined as "(start+end)/2"
                        "fraglens1":"int32","fraglens2":"int32", #fragment lengthes                        
                        "distances":"int32", #distance between fragments. If -1, different chromosomes. If -2, different arms.                         
                        "fragids1":"int64","fragids2":"int64",  #IDs of fragments. fragIDmult * chromosome + location                          
                        "dists1":"int32","dists2":"int32",        #distance to rsite
                        "cuts1":"int32","cuts2":"int32",           #precise location of cut-site 
                        "strands1":"bool","strands2":"bool",
                        "DS":"bool","SS":"bool"}
        self.saveExtension = "hdf5"

        self.genome = dnaUtils.Genome(genomeFolder)
        self.chromosomeCount = len(self.genome.chromosomes)  #used for building heatmaps
        self.fragIDmult = self.genome.fragIDmult
        print "----> New dataset opened, genome %s,  %s chromosomes" % (self.genome.type, self.chromosomeCount)

        self.maximumMoleculeLength = maximumMoleculeLength  #maximum length of a molecule for SS reads

        self.folder = folder #folder to save the data. May be empty.         
         
        if os.path.isfile(self.folder) == True:
            self.exitProgram("Specified folder is a file")        
        
        if os.path.exists(os.path.join(self.folder,"%s.%s" % (self.vectors.keys()[0],self.saveExtension))):
            if override == False: 
                print "----->!!!Another dataset found in the folder. It will be used, or overriden if loaded." 
                print "     If you want to use it, be sure it is consistent."
                print "     Otherwise, loading datasets is adviced as it checks for consistency"
                print
            else:
                self.delete(angry = False)  
            
        if os.path.exists(self.folder) == False:
            os.mkdir(self.folder)
            print "Folder created: %s" % (self.folder,)
                         
        

    def setData(self,name,data,folder = None):
        "a method to save numpy arrays to HDD quickly"
        if folder == None: folder = self.folder        
        if name not in self.vectors.keys():
            raise
        dtype = numpy.dtype(self.vectors[name]) 
        if data.dtype != dtype:
            data = numpy.asarray(data,dtype=dtype)
        f = h5py.File(os.path.join(folder , '%s.%s' % (name, self.saveExtension))  ,'w')
        f.create_dataset("MyDataset",data = data,compression = "lzf")        
        f.close()
    
    def getData(self,name,folder = None):
        "a method to load numpy arrays from HDD quickly"
        if folder == None: folder = self.folder 
        if name not in self.vectors.keys():
            raise                
        try: 
            f = h5py.File(os.path.join(folder , '%s.%s' % (name,self.saveExtension)) ,'r')
        except IOError:             
            print "cannot open file, ", os.path.join(folder , '%s.%s' % (name,self.saveExtension))
            raise IOError("HDF5 file do not exist. Are you loading existing dataset?")
        data = numpy.array(f["MyDataset"])
        f.close()
        return data
    
    def __getattribute__(self,x):
        "a method that overrides set/get operation for self.vectors so that they're always on HDD"
        if x == "vectors": return object.__getattribute__(self,x)
        
        if x in self.vectors.keys():            
            a =  self.getData(x)            
            return a        
         
        else:
            return object.__getattribute__(self,x)

    def __setattr__(self,x,value):
        "a method that overrides set/get operation for self.vectors so that they're always on HDD"
        if x == "vectors": return object.__setattr__(self,x,value)
        
        if x in self.vectors.keys():        
            self.setData(x,value)            
            
        else:
            return object.__setattr__(self,x,value)
        
    
    def merge(self,folders):
        "combines data from multiple datasets"
        if self.folder in folders:
            self.exitProgram("----> Cannot merge folder into itself! Create a new folder")        
        for name in self.vectors.keys():
            res = []
            for folder in folders:
                res.append(self.getData(name,folder))
            res = numpy.concatenate(res)
            self.setData(name,res)

    
    def parseAnton(self,filename):
        "This method unpickles files saved by Anton"
        #data = shelve.open(filename)                        
        #self.parsedfilename = filename+".dat"
        
        data = dict(numpy.load(filename))      
        print "----> File loaded: ", filename, "total reads: ", len(data["chrms1"])        

        DE = (data['rfrags1']  == data['rfrags2'] ) * (data["chrms1"] == data["chrms2"])  #filtering out same fragment reads  
        mask = ((DE == False)  * (data['chrms1'] <= self.chromosomeCount) * (data["chrms2"] <= self.chromosomeCount)) #filtering chromosomes        
        dist = data["cuts1"] * (2 * data["dirs1"]-1) + data["cuts2"] * (2 * data["dirs2"] - 1) #distance between reads pointing to each other
        readsMolecules = (data["chrms1"] == data["chrms2"])*(data["dirs1"] != data["dirs2"]) *  (dist <=0) * (dist >= -self.maximumMoleculeLength)#filtering out DE 
        mask *= (readsMolecules  == False)
        print "     number of reads ---> (<=%d) <--- : " % (self.maximumMoleculeLength),readsMolecules.sum() 
        print "     resulting number of reads: ", mask.sum()

        for i in data.keys():
            data[i] = data[i][mask]
        mask = data['chrms1'] == 0  #moving SS reads to start at first side only 
        variables = set([i[:-1] for i in data.keys()])
        for i in variables:
            exec("data['%s1'][mask]  = data['%s2'][mask]" % (i,i))
            exec("data['%s2'][mask] = 0" % (i))
        
        up1 = data["uprsites1"]
        up2 = data["uprsites2"]
        down1 = data["downrsites1"]
        down2 = data["downrsites2"]
        rsites1 = data["rsites1"]
        rsites2 = data["rsites2"]
        pos1 = data["cuts1"]
        pos2 = data["cuts2"]
        self.setData("chrms1",data['chrms1'])        
        self.setData("chrms2",data['chrms2'])                
        self.setData("DS",self.getData("chrms2") != 0)
        self.setData("SS",self.getData("DS") == False) 
        self.setData("strands1", data["dirs1"])
        self.setData("strands2",data["dirs2"])
        self.setData("cuts1" , data["cuts1"])
        self.setData("cuts2" ,data["cuts2"])    
        self.setData("dists1", numpy.abs(rsites1 - pos1))  #distance to the downstream (where read points) restriction site 
        self.setData("dists2", numpy.abs(rsites2 - pos2))  #distance to the downstream (where read points) restriction site
        self.setData('mids1', (up1 + down1)/2) 
        self.setData('mids2', (up2 + down2)/2)
        self.setData("fraglens1",numpy.abs(up1 - down1))
        self.setData("fraglens2", numpy.abs(up2 - down2))
        self.fragids1 = self.mids1 + numpy.array(self.chrms1,dtype = "int64") * self.fragIDmult
        self.fragids2 = self.mids2 + numpy.array(self.chrms2,dtype = "int64") * self.fragIDmult 
        distances = numpy.abs(self.getData("mids1") - self.getData("mids2"))
        distances[self.getData("chrms1") != self.getData("chrms2")] = -1
        self.setData("distances",distances)
        print "finished parsing Anton file %s \n" % filename
 
            
    def saveFragments(self):
        "saves fragment data to make correct expected estimates after applying a heavy mask"
        self.ufragmentsOriginal = numpy.array(self.ufragments)
        self.ufragmentlenOriginal = numpy.array(self.ufragmentlen)

    def originalFragments(self):
        "loads original fragments"
        self.ufragments = numpy.array(self.ufragmentsOriginal)
        self.ufragmentlen = numpy.array(self.ufragmentlenOriginal)

    def buildHeatmap(self,chromosome = 36,chromosome2 = None,resolution = 1000000,show = True,weights = False):
        "builds two chromosome heatmap, possibly accepts weights"

        
        if weights == True: 
            try:
                self.weights
            except:
                print "!!!!!Forced calculation of weights"
                self.calculateWeights()
        if chromosome2 == None: 
            chromosome2 = chromosome                        
            mask = (self.chrms1 == chromosome) * (self.chrms2 == chromosome)
            p1 = self.mids1[mask]
            p2 = self.mids2[mask]
            p1 = na(p1,int)
            p2 = na(p2,int)
            b1 = numpy.arange(0,max(p1.max()+resolution,p2.max()+resolution),resolution)            
            hist = numpy.histogram2d(p1,p2,(b1))[0]
            hist = hist + numpy.transpose(hist)
            
        else:
            mask = (self.chrms1 == chromosome) * (self.chrms2 == chromosome2)
            p11 = self.mids1[mask]
            p21 = self.mids2[mask]            
            mask = (self.chrms1 == chromosome2) * (self.chrms2 == chromosome)
            p12 = self.mids2[mask]
            p22 = self.mids1[mask]            
            p1 = numpy.r_[p11,p12]
            p2 = numpy.r_[p21,p22]            
            if (len(p1) == 0) or (len(p2) == 0):                 
                return numpy.zeros((500,500),float) + 0.000001
            b1 = numpy.arange(0,(p1.max()+resolution),resolution)            
            b2 = numpy.arange(0,(p2.max()+resolution),resolution)
            l1 = len(b1)
            l2 = len(b2)
            hist = numpy.histogram2d(p1,p2,(b1,b2))[0]
        "now calculating weights"
        if chromosome2 == None: chromosome2 = chromosome
        try: self.weights
        except:weights = False
        if weights == True:            
            m1 = self.ufragments / self.fragIDmult == chromosome
            m2 = self.ufragments / self.fragIDmult == chromosome2
            p1 = self.ufragments[m1]  % self.fragIDmult 
            p2 = self.ufragments[m2]  % self.fragIDmult
            w1 = self.weights[m1]
            w2 = self.weights[m2]
            p1mod = p1 / resolution
            p2mod = p2 / resolution            
            myarray = numpy.array(range(numpy.max(p1mod)+1))
            vec1 =  arraySumByArray(p1mod,myarray,w1)
            myarray = numpy.array(range(numpy.max(p2mod)+1))
            vec2 = arraySumByArray(p2mod,myarray,w2)
            vec1 = vec1[:l1]
            vec2 = vec2[:l2]
            correction = na(vec1[:,None] * vec2[None,:],float)
            mask = numpy.logical_or(numpy.isnan(correction), correction == 0) 
            correction[mask] = 1 
            hist2 = hist / correction
            hist2 /= (numpy.mean(hist2[mask==False]) / numpy.mean(hist[mask==False]))
            if show == True: mat_img(numpy.log(numpy.array(hist2,float) ))
            else: return numpy.array(hist2 )
        if show == True: mat_img(numpy.log(numpy.array(hist,float) ))
        else: return numpy.array(hist)
    
    def buildAllHeatmap(self,resolution):
        "creates an all-by-all heatmap in accordance with mapping provided by 'genome' class" 
        self.genome.createMapping(resolution)
        dr = self.DS 
        label1 = self.genome.chromosomeStarts[self.chrms1[dr] - 1] + self.mids1[dr] / resolution
        label1 = numpy.array(label1, dtype = "uint32")
        label2 = self.genome.chromosomeStarts[self.chrms2[dr] - 1] + self.mids2[dr] / resolution
        label2 = numpy.array(label2, dtype = "uint32")       
        label = label1 * numpy.int64(self.genome.N) + label2
        del label1
        del label2
        if self.genome.N < 65000:
            label = numpy.array(label, dtype = "uint32")                    
         
        counts = numpy.bincount(label, minlength = self.genome.N**2)
        if len(counts) > self.genome.N**2:
            print "heatmap exceed length of the genome!!! Check genome"
            exit()
            
        counts.shape = (self.genome.N,self.genome.N)
        for i in xrange(len(counts)):
            counts[i,i:] += counts[i:,i]
            counts[i:,i] = counts[i,i:]        
        return counts 
    
    def buildSinglesCoverage(self,resolution):
        "creates an SS coverage vector heatmap in accordance with the output of the 'genome' class"
        self.genome.createMapping(resolution)
        ds = self.DS == False        
        label = self.genome.chromosomeStarts[self.chrms1[ds] - 1] + self.mids1[ds] / resolution
        counts = sumByArray(label, numpy.arange(self.genome.N))
        return counts
    
    def buildFragmetCoverage(self,resolution):
        "creates HindIII density vector (visible sites only) heatmap in accordance with the output of the 'genome' class"         
        self.genome.createMapping(resolution)
        try: self.ufragments
        except: self.rebuildFragments()
        chroms = self.ufragments / self.fragIDmult
        positions = self.ufragments % self.fragIDmult
        label = self.genome.chromosomeStarts[chroms - 1] + positions / resolution
        counts = sumByArray(label, numpy.arange(self.genome.N))
        return counts
        

    
    def splitFragmentsBystrands(self):
        "Splits fragments: those with strands = 1 gets location += 100. This is totally safe!  "
        "This might be fun if you want to analyze two sides of the fragment separately, but unnecessary otherwise"
        f1 = self.fragids1
        f1 += ((f1 % self.fragIDmult) % 2)
        f1 [self.strands1 == 1] += 1 
        f2 = self.fragids2
        f1 += ((f2 % self.fragIDmult) % 2)
        f2 [self.strands1 == 1] += 1
        self.rebuildFragments()
        

    
    def filterCentromeres(self):
        "separates chromosomal arms as different chromosomes"
        #mask = self.chrms1 == self.chrms2
        centromeres = self.genome.centromeres        
        centromeres = r_[1000000000,centromeres]
        mids1 = self.mids1 /  centromeres[self.chrms1]            
        tochange = (mids1 > 1)         
        self.chrms1[tochange] += self.chromosomeCount
        mids2 = self.mids2 /  centromeres[self.chrms2]    
        tochange = (mids2 > 1)         
        self.chrms2[tochange] += self.chromosomeCount
        locationsSingles = self.locationsSingles /  centromeres[self.chromosomesSingles]            
        tochange = (locationsSingles > 1)         
        self.chromosomesSingles[tochange] += self.chromosomeCount
        self.distances[numpy.abs(self.chrms1 - self.chrms2) == self.chromosomeCount] = -2
        self.fragids1 = self.mids1  + numpy.array(self.chrms1 * self.fragIDmult,"int64")
        self.fragids2 = self.mids2  + numpy.array(self.chrms2 * self.fragIDmult,"int64")
        self.rebuildFragments()
        
    def recoverCentromeres(self):
        "recovers whole chromosomes after chromosomal arm separation" 
        self.chrms1[self.chrms1 > self.chromosomeCount] -= self.chromosomeCount
        self.chrms2[self.chrms2 > self.chromosomeCount] -= self.chromosomeCount
        self.rebuildFragments()
      
    def fragmentFilter(self,fragments):
        "keeps only reads that originate from fragments in 'fragments' variable, for DS - on both sides"
        if fragments.dtype == numpy.bool:
            fragments = self.ufragments[fragments]        
        m1 = arrayInArray(self.fragids1,fragments)
        m2 = arrayInArray(self.fragids2,fragments) + self.SS
        mask = numpy.logical_and(m1,m2)
        self.maskFilter(mask)

    
    def maskFilter(self,mask):
        "keeps only reads designated by mask"
        print "          Number of reads changed  %d ---> %d" % (len(mask),mask.sum()),
        length = 0 
        for name in self.vectors:
            data = self.getData(name)
            ld = len(data)
            if length == 0: length = ld
            else:
                if ld != length: 
                    self.delete()             
            self.setData(name,data[mask])  
        self.rebuildFragments()
        
    def rebuildFragments(self):
        "recreates a set of fragments - runs when read have changed"
        try: 
            past = len(self.ufragments)        
        except:
            past = 0        
        self.ufragids1,self.ufragids1ind = numpy.unique(self.fragids1,return_index=True)
        self.ufragids2,self.ufragids2ind = numpy.unique(self.fragids2[self.DS],return_index=True)
                
        #Funding unique fragments and unique fragment IDs
        self.ufragment1len = self.fraglens1[self.ufragids1ind]
        self.ufragment2len = self.fraglens2[self.DS][self.ufragids2ind]
 
        uall = numpy.r_[self.ufragids1,self.ufragids2]
        ulen = numpy.r_[self.ufragment1len,self.ufragment2len]
         
        self.ufragments,ind = numpy.unique(uall, True)        
        self.ufragmentlen = ulen[ind]            
        print "          Fragments number changed -   %d --->  %d" % (past, len(self.ufragments))
        
    def calculateWeights(self):
        "calculates weights for reads based on fragment length correction similar to Tanay's, may be used for scalings or creating heatmaps"        
        fragmentLength = self.ufragmentlen
        pls = numpy.sort(fragmentLength)
        pls = numpy.r_[pls,pls[-1]+1]
        N = len(fragmentLength)    
        mysum = numpy.array(self.fragmentSum(),float)
        self.weights = numpy.ones(N,float)
        meanSum = numpy.mean(numpy.array(mysum,float))
        #watch = numpy.zeros(len(mysum),int)
        for i in numpy.arange(0,0.991,0.01):
            b1,b2 = pls[i*N],pls[(i+0.01)*N]
            p = (b1 <= fragmentLength)* (b2 > fragmentLength)
            #watch[p] += 1                        
            value = numpy.mean(mysum[p])      
            if p.sum() > 0:
                self.weights[p] =  value / meanSum
            else:
                print "no weights",i,b1,b2

            

    def filterExtreme(self,cutH = 0.005,cutL = 0):
        "removes fragments with most and/or least # counts"
        print "----->Extreme fragments filter: remove top %lf, bottom %lf fragments" % (cutH, cutL)
        s  = self.fragmentSum()
        ss = numpy.sort(s)
        print "     #Top fragments are: ",ss[-10:]                 
        N = len(ss)
        print "     # Cutoff for low # counts is (counts): ",ss[int(cutL * N)],"; cutoff for large # counts is: ",ss[int((1-cutH)*N)]         
        news = (s >= ss[int(cutL * N)]) * (s <= ss[int((1-cutH)*N)])        
        self.fragmentFilter(self.ufragments[news])
        print
        
    def filterLarge(self,cutlarge = 100000,cutsmall = 100 ):
        "removes very large and small fragments"
        print "----->Small/large fragments filter: keep strictly less than %d, strictly more than %d bp" % (cutlarge, cutsmall)                
        p = (self.ufragmentlen < (cutlarge ) ) * (self.ufragmentlen > cutsmall)                      
        self.fragmentFilter(self.ufragments[p])
        print 

    
    def filterRsiteStart(self,offset = 5):
        "removes reads that start within x bp near rsite"
        print "----->Semi-dangling end filter: remove guys who start %d bp near the rsite" % offset
#        d1 = self.dist1
#        l1 = self.fraglens1        
#        d2 = self.dist2
#        l2 = self.fraglens2
#        ds = self.DS
#        ss = self.SS
#        mask = numexpr.evaluate("(abs(d1 - l1) >= offset) and (((abs(d2 - l2) >= offset) and ds)  or  ss)")   #is buggy        
                
        mask = (numpy.abs(self.dists1 - self.fraglens1) >=offset) * ((numpy.abs(self.dists2 - self.fraglens2) >= offset )* self.DS + self.SS)
        self.maskFilter(mask)
        print
        
    def filterDuplicates(self):
        "removes duplicate molecules in DS reads"
        print "----->Filtering duplicates in DS reads: "
        l1 = numpy.array(self.cuts1)
        l2 = numpy.array(self.cuts2)
        ch1 = self.chrms1
        ch2 = self.chrms2
        n1 = numpy.int64(715827883)
        cc = numpy.int64(self.chromosomeCount)
        cc,n1 #Eclipse warning         
        uid = numexpr.evaluate("(l1 * cc + ch1) * n1 + (l2 * cc + ch2)")
        del l1
        del l2
        del ch1
        del ch2
        un = numpy.unique(uid,True)[1]
        del uid         
        stay = self.SS.copy()
        stay[un] = True
        ds = self.DS.sum()
        print "     Number of DS reads changed - %d ---> %d" % (ds,ds - len(self.DS) + stay.sum()) 
        del un
        self.maskFilter(stay)
        print
        

    def fragmentSum(self,fragments = None, strands = "both"):
        "returns sum of all counts for a set or subset of fragments"
        if fragments == None: fragments = self.ufragments                
        if strands == "both":  
            return sumByArray(self.fragids1,fragments) + sumByArray(self.fragids2[self.DS],fragments) 
        if strands == 1: return sumByArray(self.fragids1,fragments)
        if strands == 2:return sumByArray(self.fragids2[self.DS],fragments)
        
        
    def printStats(self):
        print "-----> Statistics for the folder %s!" % self.folder
        print "     Single sided reads: " ,self.SS.sum()
        print "     Double sided reads: " , self.DS.sum()
        ss1 = self.strands1[self.chrms1>0]
        ss2 = self.strands2[self.chrms2>0]
        sf = ss1.sum() + ss2.sum()
        sr = len(ss1) + len(ss2) - sf
        print "     reverse/forward bias",float(sr)/sf
        
    def delete(self,folder = None, angry = True):         
        if folder == None: folder = self.folder
        if angry == True: print "     Folder  %s removed due to inconsistency" % folder
        else: print "     Override: Folder %s cleaned up" % (folder,)
        for i in os.listdir(folder):
            os.remove(os.path.join(folder,i))
                    
        
    def save(self,folderName):
        if self.folder == folderName: self.exitProgram("Cannot save to the working folder")
        try: os.mkdir(folderName)
        except: print "folder exists", folderName
        for name in self.vectors.keys(): self.setData(name, self.getData(name),folderName)
        print "----> Data saved to folder %s" % (folderName,)
    
    def load(self,folderName):
        if folderName == self.folder: self.exitProgram("Cannot load from the working folder")
        length = 0 
        for name in self.vectors:
            data = self.getData(name,folderName)
            ld = len(data)
            if length == 0: 
                length = len(data)
            else:
                if ld != length: 
                    print("---->!!!!!Folder %s contains inconsistend data<----" % folderName)
                    print("     Both folders will be removed")
                    self.delete()
                    self.delete(folderName)
                    self.exitProgram("----> Data removed! Sorry...")
                    
            self.setData(name,data,self.folder) 
        print "---->Loaded data from folder %s, contains %d reads" % (folderName, length)
            

        
    def saveHeatmap(self,filename,resolution = 1000000):
        
        try: os.remove(filename)
        except: pass        
        heatmap = self.buildAllHeatmap(resolution)        
        singles = self.buildSinglesCoverage(resolution)        
        frags = self.buildFragmetCoverage(resolution)
        chromosomeStarts = numpy.array(self.genome.chromosomeStarts)
        tosave = {}
        tosave["resolution"] = resolution
        tosave["heatmap"] = heatmap
        tosave["singles"] = singles
        tosave["frags"] = frags
        tosave["genome"] = self.genome.type
        tosave["chromosomeSTarts"] = chromosomeStarts
        saveData(tosave,filename)
        print "----> Heatmap saved to '%s' at %d resolution" % (filename,resolution)
        
    def exitProgram(self,a):
        print a
        print "     ----> Bye! :) <----"
        exit()

class HiCStatistics(track):
    "a sub-class of a 'track' class used to do statistics on Hi-C reads" 

    def multiplyForTesting(self,N = 100):
        "used for heavy load testing only"
        for name in self.vectors:
            print "multipliing",name
            one = self.getData(name)
            if name in ["mids1","mids2"]:
                one += numpy.random.randint(0,5*N,len(one))            
            blowup = numpy.hstack(tuple([one]*N))
            self.setData(name,blowup)
            
    
    def buildLengthDependencePlot(self,label = "plot", strands = "both",color = None):
        "plots dependence of counts on fragment length. May do based on one strands only"
        "please run  plt.legend & plt.show() after calling this for all datasets you want to consider"     
        fragmentLength = self.ufragmentlen
        pls = numpy.sort(fragmentLength)
        N = len(fragmentLength)
        sums = []
        sizes = []            
        mysum = self.fragmentSum(None,strands)
    
        for i in numpy.arange(0,0.98,0.015):
            b1,b2 = pls[i*N],pls[(i+0.015)*N-1]
            p = (b1 < fragmentLength)* (b2 > fragmentLength)
            value = numpy.mean(mysum[p])            
            sums.append(value)
            sizes.append(numpy.mean(fragmentLength[p]))    
        if color == None: plt.plot(sizes,sums,'x-',markersize = 3,linewidth = 2,label = label)
        else:  plt.plot(sizes,sums,color,linewidth = 2,label = label)
    
        
    def plotScaling(self,fragids1 = None,fragids2 = None,   #IDs of fragments for which to plot scaling. 
                    #One can, for example, limit oneself to only fragments shorter than 1000 bp 
                    #Or calculate scaling only between different arms
                    color = None, label = "", #plot parameters
                    useWeights = False ,        #use weights associated with fragment length
                    excludeNeighbors = None, enzyme = None,   #number of neighboring fragments to exclude. Enzyme is needed for that!   
                    normalize = True,          #normalize the final plot to sum to one
                    withinArms = True,                #Treat chromosomal arms separately
                    mindist = 10000,  #Scaling was proved to be unreliable under 10000 bp for 6-cutter enzymes
                    #----------Calculating scaling within a set of regions only------
                    regions = None   #Array of tuples (chrom, start, end) for which scaling should be calculated
                    #Note that calculation might be extremely long (it might be proportional to # of regions for # > 100)
                                           
                    ):               #Sad smiley, because this method is very painful and complicated
        """plots scaling over, possibly uses subset of fragmetns, or weigts, possibly normalizes after plotting
        
        Plan of scaling calculation: 
        
        1. Subdivide all genome into regions. 
            a. Different chromosomes
            b. Different arms
            c. User defined squares/rectangles on a contact map
               -(chromosome, start,end) square around the diagonal 
               -(chr1, start1, end1, start2, end2) rectangle
        
        2. Use either all fragments, or only interactions between two groups of fragments 
            e.g. you can calculate how scaling for small fragments is different from that for large
            It can be possibly used for testing Hi-C protocol issues.
            One can see effect of weights by doing this
        
        3. (optional) Calculate correction associated with fragment length dependence 
        
        4. Subdivide all possible genomic separation into log-spaced bins 
        
        5. Calculate expected number of fragment pairs within each bin (possibly with weights from step 3). 
            If exclusion of neighbors is specificed, expected number of fragments knows about this                   
         
        """
        if excludeNeighbors <= 0: excludeNeighbors = None   #Not excluding neighbors  
        #use all fragments if they're not specified 
        if fragids1 == None: fragids1 = self.ufragments
        if fragids2 == None: fragids2 = self.ufragments
        
        #Calculate regions if not specified 
        if regions == None: 
            if withinArms == False: 
                regions = [(i,0,self.genome.chromosomes[i-1]) for i in xrange(1,self.genome.chromosomeCount+1)]
            else:
                regions = [(i,0,self.genome.centromeres[i-1]) for i in xrange(1,self.genome.chromosomeCount+1)] + [(i,self.genome.centromeres[i-1],self.genome.chromosomes[i-1]) for i in xrange(1,self.genome.chromosomeCount+1)]
                
        maxdist = max ( max([i[2] - i[1] for i in regions]),  #normal regions 
                        max([abs(i[2] - i[3]) for i in regions if len(i) > 3] + [0]),   #rectangular regions
                        max([abs(i[1] - i[4]) for i in regions if len(i) > 3] + [0]))

        regionID = numpy.zeros(len(self.chrms1),numpy.int16) - 1  #Region to which a read belongs             
        chr1 = self.chrms1
        chr2 = self.chrms2
        pos1 = self.mids1
        pos2 = self.mids2
        fragRegions1 = numpy.zeros(len(fragids1),int)
        fragRegions2 = numpy.zeros(len(fragids2),int)
        fragch1 =  fragids1 / self.fragIDmult
        fragch2 =  fragids2 / self.fragIDmult
        fragpos1 = fragids1 % self.fragIDmult
        fragpos2 = fragids2 % self.fragIDmult            
        
        for regionNum,region in enumerate(regions):
            if len(region) == 3: 
                chrom,start1,end1  = region
                mask = (chr1 == chrom) * (pos1 > start1) * (pos1 < end1) * (chr2 == chrom) * (pos2 > start1) * (pos2 < end1)                                     
                regionID[mask] = regionNum 
                mask1 = (fragch1 == chrom) * (fragpos1 > start1) * (fragpos1 < end1)
                mask2 = (fragch2 == chrom) * (fragpos2 > start1) * (fragpos2 < end1)   
                fragRegions1[mask1] = regionNum
                fragRegions2[mask2] = regionNum  

            if len(region) == 5: 
                chrom,start1,end1,start2,end2 = region
                mask1 = (chr1 == chrom) * (chr2 == chrom) * (pos1 > start1) * (pos1 < end1) * (pos2 > start2) * (pos2  < end2)
                mask2 = (chr1 == chrom) * (chr2 == chrom) * (pos1 > start2) * (pos1 < end2) * (pos2 > start1) * (pos2  < end1)
                mask = mask1 + mask2
                regionID[mask] = regionNum                    
                mask1 = (fragch1 == chrom) * (fragpos1 > start1) * (fragpos1 < end1)
                mask2 = (fragch2 == chrom) * (fragpos2 > start2) * (fragpos2 < end2)   
                fragRegions1[mask1] = regionNum
                fragRegions2[mask2] = regionNum
        del chr1,chr2,pos1,pos2                        

                                        
        lens = numpy.array(numutils.logbins(mindist,maxdist,1.25),float)+0.1   #bins of lengths        
        
        positions = []   #arrays to write the result 
        values = []
        
        if excludeNeighbors != None: 
            if enzyme == None: raise ValueError("Please specify enzyme if you're excluding Neighbors")
            fragmentDists = self.genome.getFragmentDistance(self.fragids1,self.fragids2, enzyme)
            mask3 = fragmentDists > excludeNeighbors   #keep only guys more than excludeNeighbors fragments apart 
        
        #Keeping reads for fragments in use
        p11 = arrayInArray(self.fragids1,fragids1)
        p12 = arrayInArray(self.fragids1,fragids2)
        p21 = arrayInArray(self.fragids2,fragids1)
        p22 = arrayInArray(self.fragids2,fragids2)
        mask = (p11* p22) + (p12 * p21)   #reads between fragids1 and fragids2 
        mask2 = (mask) * (regionID >= 0) * (self.strands1 == self.strands2)
        #Reads only between fragments
        #Reads from the same region
        #Reads like --> -->    or <-- <--, discarding --> <-- and <-- -->             
        if excludeNeighbors != None: mask2 = mask2 * mask3    #remove all neighbors 
                                     
        distances = numpy.sort(self.distances[mask2])            
                
        "calculating fragments lengths for exclusions to expected # of counts"        
        #sorted fragment IDs and lengthes
        args = numpy.argsort(self.ufragments)
        usort = self.ufragments[args]        
        
        if useWeights == True:   #calculating weights if needed  
            try:      self.weights
            except:   self.calculateWeights()                
            uweights = self.weights[args]   #weights for sorted fragment IDs 
            weights1 = uweights[numpy.searchsorted(usort,fragids1)]
            weights2 = uweights[numpy.searchsorted(usort,fragids2)]   #weghts for fragment IDs under  consideration        
        
        lenmins,lenmaxs  = lens[:-1], lens[1:]        
                 
        N = len(lenmins)   #number of bins 
        count = [0 for _ in xrange(N)]      #count of reads in each min 
        chr1 = fragids1 / self.fragIDmult    
        chr2 = fragids2 / self.fragIDmult
        pos1 = fragids1 % self.fragIDmult        
        pos2 = fragids2 % self.fragIDmult
        
        for regionNumber,region  in enumerate(regions):  
            
            mask = numpy.nonzero(fragRegions1  == regionNumber)[0]
            mask2 = numpy.nonzero(fragRegions2  == regionNumber)[0]  #filtering fragments that correspont to current region
            if (len(mask) == 0) or (len(mask2) == 0): continue                             
            bp1, bp2 = pos1[mask], pos2[mask2]         #positions of fragments on chromosome
             
            p2arg = numpy.argsort(bp2)   
            p2 = bp2[p2arg]           #sorted positions on the second fragment

            if excludeNeighbors != None:
                "calculating excluded fragments (neighbors) and their weights to subtract them later" 
                excFrag1, excFrag2 = self.genome.getPairsLessThanDistance(fragids1[mask] , fragids2[mask2], excludeNeighbors, enzyme)                                                
                excDists = numpy.abs(excFrag2 - excFrag1)  #distances between excluded fragment pairs
                if useWeights == True:        
                    correctionWeights = weights1[numutils.arraySearch(fragids1,excFrag1)]
                    correctionWeights = correctionWeights *  weights2[numutils.arraySearch(fragids2,excFrag2)]   #weights for excluded fragment pairs            
            
            if useWeights == True: 
                w1,w2 = weights1[mask], weights2[mask2]                                  
                sw2 = numpy.r_[0,numpy.cumsum(w2[p2arg])]    #cumsum for sorted weights on 2 strand             
            
            for lenmin,lenmax,index in zip(lenmins,lenmaxs,range(N)):
                "Now calculating actual number of fragment pairs for a length-bin, or weight of all these pairs"                
                mymin,mymax = bp1 - lenmax, bp1 - lenmin   #calculating boundaries for fragments on a second strand                  
                val1 = numpy.searchsorted(p2,mymin)  #Calculating indexes on the second strand 
                val2 = numpy.searchsorted(p2,mymax)
                if useWeights == False: curcount = numpy.sum(numpy.abs(val1 - val2))   #just # of fragments 
                else: curcount = numpy.sum(w1 * numpy.abs(sw2[val1] - sw2[val2]))      #(difference in cumsum of weights) * my weight 
                
                mymin,mymax = bp1 + lenmax, bp1 + lenmin   #repeating the same for 
                val1 = numpy.searchsorted(p2,mymin)
                val2 = numpy.searchsorted(p2,mymax)
                if useWeights == False: curcount += numpy.sum(numpy.abs(val1 - val2))
                else: curcount += numpy.sum(w1 * numpy.abs(sw2[val1] - sw2[val2]))
                
                if excludeNeighbors != None: #now modifying expected count because of excluded fragments 
                    if useWeights == False: 
                        ignore = ((excDists > lenmin) * (excDists < lenmax)).sum()
                    else:                        
                        ignore = (correctionWeights[((excDists > lenmin) * (excDists < lenmax))]).sum()
                    
                    if (ignore >= curcount) and (ignore != 0):
                        if ignore < curcount * 1.0001:
                            curcount = ignore = 0   #Check for error and numerical instabilities  
                        else:   print "error found" , "lenmin:",lenmin, "  curcount:",curcount, "  ignore:",  ignore
                    else: #Everything is all right                                                                      
                        curcount  -= ignore
                count[index] += curcount 
        maxcountsall = count        
         
        for i in xrange(len(lens) - 1):   #Dividing observed by expected
            beg,end  = lens[i], lens[i+1]            
            first,last = tuple(numpy.searchsorted(distances[:-1],[beg,end]))            
            mycounts = last - first
            maxcounts = maxcountsall[i]
                                     
            positions.append(sqrt(float(beg)*float(end)))
            values.append(mycounts/float(maxcounts))
        positions = numpy.array(positions)
        values = numpy.array(values)
         
        if normalize == True: values /= numpy.sum(1. * (positions * values) [numpy.logical_not(numpy.isnan(positions * values))])        
        if color!= None:
            plt.plot(positions,values,color , label = label+", %d reads" % len(distances))
        else: 
            plt.plot(positions,values, label = label+", %d reads" % len(distances))
        return (positions, values) 
        
    
    def plotRsiteStartDistribution(self,useSSReadsOnly = False,offset = 5,length = 200):
        if useSSReadsOnly == True:
            mask = self.SS
        else:
            mask = self.DS
        dists1 = self.fraglens1 - numpy.array(self.dists1,dtype = "int32") 
        dists2 = self.fraglens2 - numpy.array(self.dists2,dtype = "int32") 
        m = min(dists1.min(),dists2.min())
        if offset < -m: 
            offset = -m
            print "minimum negative distance is %d, larger than offset; offset set to %d" % (m,-m)
        dists1 += offset
        dists2 += offset          
        myrange = numpy.arange(-offset,length-offset)
        
        plt.subplot(141)
        plt.title("strands1, side 1")
        plt.plot(myrange,numpy.bincount(5+dists1[mask][self.strands1[mask] == True ])[:length])
        plt.subplot(142)
        plt.title("strands1, side 2")
        plt.plot(myrange,numpy.bincount(dists2[mask][self.strands1[mask] == True ])[:length])
        
        plt.subplot(143)
        plt.title("strands2, side 1")
        plt.plot(myrange,numpy.bincount(dists1[mask][self.strands1[mask] == False ])[:length])
        plt.subplot(144)
        plt.title("strands2, side 2")
        plt.plot(myrange,numpy.bincount(dists2[mask][self.strands1[mask] == False ])[:length])
        
        plt.show()







def plotFigure2c():
    TR = track()
    TR.load("GM-all.refined")
    hm = TR.buildHeatmap(1, 1, 1000000, False,False)
    TR.calculateWeights()    
    TR.weights = numpy.ones(len(TR.weights),float)   #if you want to correct just by fragment density, not by length dependence 
    hm2 = TR.buildHeatmap(1, 1, 1000000, False,weights = True )
    hm2[numpy.isnan(hm2)] = 0
    mask = numpy.sum(hm,axis = 0) > 0
    """p1-6 are 6 lines to be plotted, below is plotting only"""
    p1 = numpy.sum(hm, axis = 0)[mask]
    p3 = numpy.sum(correct(hm),axis = 0)[mask]
    p5 = numpy.sum(ultracorrect(hm,40),axis = 0)[mask]
    p4 = numpy.sum(correct(hm2),axis = 0)[mask]
    p2 = numpy.sum(hm2, axis = 0)[mask]
    p6 = numpy.sum(ultracorrect(hm2,40),axis = 0)[mask]
    matplotlib.rcParams['font.sans-serif']='Arial'        
    dashstyle = (3,3)    
    plt.figure(figsize = (4,4))
    
    ax = plt.subplot(2,1,1)
    plt.xlim((0,80))
    plt.ylim((0,2))
    plt.ylabel("Total coverage", fontsize = 8)
        
    line21 = plt.plot(p1 / p1.mean(),"-",linewidth = 1,color = "#e5a826")[0]
    line22 = plt.plot(p3/p3.mean(),"--",linewidth = 1,color = "#e5a826")[0]
    line22.set_dashes(dashstyle)
    line23 = plt.plot(p5/p5.mean(), linewidth = 1, color = "grey")[0]
        
    for xlabel_i in ax.get_xticklabels(): xlabel_i.set_fontsize(8)
    for xlabel_i in ax.get_yticklabels(): xlabel_i.set_fontsize(8)                    
    legend = plt.legend([line21, line22, line23],
                        ["Raw data","Single correction","Iterative correction"],prop={"size":6},loc = 1,handlelength=2)
    legend.draw_frame(False)
    removeAxes(shift = 0,ax = ax)
    
    for i in ax.spines.values(): i.set_color('none')
    ax.axhline(linewidth=1, color='black')
    ax.axvline(linewidth=1, color='black')
    

    ax2 = plt.subplot(2,1,2,sharex = ax )
    plt.xlim((0,80))
    plt.ylim((0,2))
    plt.xlabel("Position on chom 1 (MB)", fontsize = 8)
    plt.ylabel("Total coverage", fontsize = 8)

    line1 = plt.plot(p4/p4.mean(),"--",color = "#9b3811",linewidth = 1)[0]    
    line1.set_dashes(dashstyle)
    line2 = plt.plot(p2 / p2.mean(),"-",color = "#9b3811",linewidth = 1)[0]    
    line3 = plt.plot(p6/p6.mean(), linewidth = 1, color = "grey")[0]

    for xlabel_i in ax2.get_xticklabels(): xlabel_i.set_fontsize(8)
    for xlabel_i in ax2.get_yticklabels(): xlabel_i.set_fontsize(8)
    
    legend = plt.legend([line2,line1,line3],
                        ["HindIII corrected","Single correction","Iterative correction"],prop={"size":6},loc = 1,handlelength=2)
    legend.draw_frame(False)
    removeAxes(shift = 0,ax = ax2)
    plotting.niceShow() 
    
#plotFigure2c()


def doSupplementaryCoveragePlot():
    TR = track()
    TR.load("GM-all.refined")
    s1 = TR.fragmentSum(strands = 1)
    TR.saveFragments()
    TR.maskFilter(TR.dists1 > 500)
    TR.originalFragments()
    s2 = TR.fragmentSum(strands = 1)
    resolution = 1000000
    def coverage(s1,s2,TR):
        genome = dnaUtils.Genome("HG18")
        genome.createMapping(resolution)            
        label = genome.chromosomeStarts[TR.ufragments / TR.fragIDmult - 1] + (TR.ufragments % TR.fragIDmult ) / resolution
        counts = numpy.bincount(label, weights = s1)
        counts2 = numpy.bincount(label,weights = s2)
        data = cPickle.load(open("GC1M",'rb'))
        eigenvector = numpy.zeros(genome.chromosomeEnds[-1],float)
        inds = numpy.argsort(counts)
        mask = inds[int(0.02 * len(inds)):]                
        for chrom in range(1,24):
            eigenvector[genome.chromosomeStarts[chrom-1]:genome.chromosomeStarts[chrom-1] + len(data[chrom-1])] = data[chrom-1]
        eigenvector[eigenvector < 35] = 35        
        plt.scatter(counts[mask],counts2[mask], c = eigenvector[mask],s = 6,linewidth = 0)
        print stats.spearmanr(counts[mask],counts2[mask])
        plt.xlabel("Coverage from all reads")
        plt.xticks([0,5000,10000,15000])
        plt.ylabel("Coverage from RBs")
        b = plt.colorbar()
        b.ax.set_xlabel("GC content")
    plt.subplot(121)
    plt.title("HinIII")
    coverage(s1,s2,TR)
    
    TR = track()
    TR.load("GM-NcoI.refined")
    s1 = TR.fragmentSum(strands = 1)
    TR.saveFragments()
    TR.maskFilter(TR.dists1 > 500)
    TR.originalFragments()
    s2 = TR.fragmentSum(strands = 1)
    resolution = 1000000
    plt.subplot(122)
    plt.title("NcoI")
    coverage(s1,s2,TR)
    plt.show() 

 





  
