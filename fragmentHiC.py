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
    if override == True:
        if os.path.exists(filename):
            os.remove(filename)
    joblib.dump(data,filename,compress = 3)


class track(object):
    def __init__(self,folder , genomeFolder , maximumMoleculeLength = 500,override = True):
        #These are fields that will be kept on a hard drive
        self.vectors = {"chromosomes1":"int8","chromosomes2":"int8", #chromosomes. If >chromosomeCount, then it's second chromosome arm! 
                        "locations1":"int64","locations2":"int64",  #midpoint of a fragment, determined as "(start+end)/2"
                        "lengthes1":"uint32","lengthes2":"uint32", #fragment lengthes                        
                        "distances":"uint32", #distance between fragments. If -1, different chromosomes. If -2, different arms.                         
                        "fragments1":"int64","fragments2":"int64",  #IDs of fragments. 100 * location + chromosome.                         
                        "dists1":"uint32","dists2":"uint32",
                        "cuts1":"uint32","cuts2":"uint32",
                        "strand1":"bool","strand2":"bool",
                        "DS":"bool","SS":"bool"}
        self.saveExtension = "hdf5"

        self.genome = dnaUtils.Genome(genomeFolder)
        self.chromosomeCount = len(self.genome.chromosomes)  #used for building heatmaps
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
            data = numpy.array(data,dtype=dtype)
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
            raise "HDF5 file do not exist. Are you loading existing dataset?"
            
        
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
        "merges multiple files in one"
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
        self.setData("chromosomes1",data['chrms1'])        
        self.setData("chromosomes2",data['chrms2'])                
        self.setData("DS",self.getData("chromosomes2") != 0)
        self.setData("SS",self.getData("DS") == False) 
        self.setData("strand1", data["dirs1"])
        self.setData("strand2",data["dirs2"])
        self.setData("cuts1" , data["cuts1"])
        self.setData("cuts2" ,data["cuts2"])    
        self.setData("dists1", numpy.abs(rsites1 - pos1))  #distance to the downstream (where read points) restriction site 
        self.setData("dists2", numpy.abs(rsites2 - pos2))  #distance to the downstream (where read points) restriction site
        self.setData('locations1', (up1 + down1)/2) 
        self.setData('locations2', (up2 + down2)/2)
        self.setData("lengthes1",numpy.abs(up1 - down1))
        self.setData("lengthes2", numpy.abs(up2 - down2))
        self.setData("fragments1", self.getData("locations1") * 100 + self.getData("chromosomes1"))
        self.setData("fragments2", self.getData("locations2") * 100 + self.getData("chromosomes2"))                
        
        distances = numpy.abs(self.getData("locations1") - self.getData("locations2"))
        distances[self.getData("chromosomes1") != self.getData("chromosomes2")] = -1
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
            mask = (self.chromosomes1 == chromosome) * (self.chromosomes2 == chromosome)
            p1 = self.locations1[mask]
            p2 = self.locations2[mask]
            p1 = na(p1,int)
            p2 = na(p2,int)
            b1 = numpy.arange(0,max(p1.max()+resolution,p2.max()+resolution),resolution)            
            hist = numpy.histogram2d(p1,p2,(b1))[0]
            hist = hist + numpy.transpose(hist)
            
        else:
            mask = (self.chromosomes1 == chromosome) * (self.chromosomes2 == chromosome2)
            p11 = self.locations1[mask]
            p21 = self.locations2[mask]            
            mask = (self.chromosomes1 == chromosome2) * (self.chromosomes2 == chromosome)
            p12 = self.locations2[mask]
            p22 = self.locations1[mask]            
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
            m1 = self.ufragments % 100 == chromosome
            m2 = self.ufragments % 100 == chromosome2
            p1 = self.ufragments[m1] /100 
            p2 = self.ufragments[m2] / 100 
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
        label1 = self.genome.chromosomeStarts[self.chromosomes1[dr] - 1] + self.locations1[dr] / resolution
        label1 = numpy.array(label1, dtype = "uint32")
        label2 = self.genome.chromosomeStarts[self.chromosomes2[dr] - 1] + self.locations2[dr] / resolution
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
        label = self.genome.chromosomeStarts[self.chromosomes1[ds] - 1] + self.locations1[ds] / resolution
        counts = sumByArray(label, numpy.arange(self.genome.N))
        return counts
    
    def buildFragmetCoverage(self,resolution):
        "creates HindIII density vector (visible sites only) heatmap in accordance with the output of the 'genome' class"         
        self.genome.createMapping(resolution)
        try: self.ufragments
        except: self.rebuildFragments()
        chroms = self.ufragments % 100
        positions = self.ufragments / 100 
        label = self.genome.chromosomeStarts[chroms - 1] + positions / resolution
        counts = sumByArray(label, numpy.arange(self.genome.N))
        return counts
        

    
    def splitFragmentsByStrand(self):
        "Splits fragments: those with strand = 1 gets location += 100. This is totally safe!  "
        "This might be fun if you want to analyze two sides of the fragment separately, but unnecessary otherwise"
        f1 = self.fragments1
        f1 += 100*((f1/100) % 2)
        f1 [self.strand1 == 1] += 100 
        f2 = self.fragments2
        f2 += 100*((f2/100) % 2)
        f2 [self.strand1 == 1] += 100
        self.rebuildFragments()
        

    
    def filterCentromeres(self):
        "separates chromosomal arms as different chromosomes"
        #mask = self.chromosomes1 == self.chromosomes2
        centromeres = self.genome.centromeres        
        centromeres = r_[1000000000,centromeres]
        locations1 = self.locations1 /  centromeres[self.chromosomes1]            
        tochange = (locations1 > 1)         
        self.chromosomes1[tochange] += self.chromosomeCount
        locations2 = self.locations2 /  centromeres[self.chromosomes2]    
        tochange = (locations2 > 1)         
        self.chromosomes2[tochange] += self.chromosomeCount
        locationsSingles = self.locationsSingles /  centromeres[self.chromosomesSingles]            
        tochange = (locationsSingles > 1)         
        self.chromosomesSingles[tochange] += self.chromosomeCount
        self.distances[numpy.abs(self.chromosomes1 - self.chromosomes2) == self.chromosomeCount] = -2
        self.fragments1 = self.locations1 * 100 + self.chromosomes1
        self.fragments2 = self.locations2 * 100 + self.chromosomes2
        self.rebuildFragments()
        
    def recoverCentromeres(self):
        "recovers whole chromosomes after chromosomal arm separation" 
        self.chromosomes1[self.chromosomes1 > self.chromosomeCount] -= self.chromosomeCount
        self.chromosomes2[self.chromosomes2 > self.chromosomeCount] -= self.chromosomeCount
        self.rebuildFragments()
      
    def fragmentFilter(self,fragments):
        "keeps only reads that originate from fragments in 'fragments' variable, for DS - on both sides"
        if fragments.dtype == numpy.bool:
            fragments = self.ufragments[fragments]        
        m1 = arrayInArray(self.fragments1,fragments)
        m2 = arrayInArray(self.fragments2,fragments) + self.SS
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
        self.ufragments1,self.ufragments1ind = numpy.unique(self.fragments1,return_index=True)
        self.ufragments2,self.ufragments2ind = numpy.unique(self.fragments2,return_index=True)
                
        #Funding unique fragments and unique fragment IDs
        self.ufragment1len = self.lengthes1[self.ufragments1ind]
        self.ufragment2len = self.lengthes2[self.ufragments2ind]
 
        uall = numpy.r_[self.ufragments1,self.ufragments2]
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
#        l1 = self.lengthes1        
#        d2 = self.dist2
#        l2 = self.lengthes2
#        ds = self.DS
#        ss = self.SS
#        mask = numexpr.evaluate("(abs(d1 - l1) >= offset) and (((abs(d2 - l2) >= offset) and ds)  or  ss)")   #is buggy        
                
        mask = (numpy.abs(self.dists1 - self.lengthes1) >=offset) * ((numpy.abs(self.dists2 - self.lengthes2) >= offset )* self.DS + self.SS)
        self.maskFilter(mask)
        print
        
    def filterDuplicates(self):
        "removes duplicate molecules in DS reads"
        print "----->Filtering duplicates in DS reads: "
        l1 = numpy.array(self.cuts1)
        l2 = numpy.array(self.cuts2)
        ch1 = self.chromosomes1
        ch2 = self.chromosomes2
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
        

    def fragmentSum(self,fragments = None, strand = "both"):
        "returns sum of all counts for a set or subset of fragments"
        if fragments == None: fragments = self.ufragments                
        if strand == "both":  
            return sumByArray(self.fragments1,fragments) + sumByArray(self.fragments2[self.DS],fragments) 
        if strand == 1: return sumByArray(self.fragments1,fragments)
        if strand == 2:return sumByArray(self.fragments2[self.DS],fragments)
        
        
    def printStats(self):
        print "-----> Statistics for the folder %s!" % self.folder
        print "     Single sided reads: " ,self.SS.sum()
        print "     Double sided reads: " , self.DS.sum()
        ss1 = self.strand1[self.chromosomes1>0]
        ss2 = self.strand2[self.chromosomes2>0]
        sf = ss1.sum() + ss2.sum()
        sr = len(ss1) + len(ss2) - sf
        print "     reverse/forward bias",float(sr)/sf
        
    def delete(self,folder = None, angry = True):         
        if folder == None: folder = self.folder
        if angry == True: print "     Folder  %s removed due to inconsistency" % folder
        else: print "     Override: Folder %s cleaned up" % (folder,)        
        for name in self.vectors:
            try: 
                filename = os.path.join(folder,name+"."+self.saveExtension)
                os.remove(filename)
                #print "removed inconsistent file %s" % filename
            except:
                filename = os.path.join(folder,name+"."+self.saveExtension)
                #print "inconsistent file not found: %s" % filename
        
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
            if name in ["locations1","locations2"]:
                one += numpy.random.randint(0,5*N,len(one))            
            blowup = numpy.hstack(tuple([one]*N))
            self.setData(name,blowup)
            
    
    def buildLengthDependencePlot(self,label = "plot", strand = "both",color = None):
        "plots dependence of counts on fragment length. May do based on one strand only"
        "please run  plt.legend & plt.show() after calling this for all datasets you want to consider"     
        fragmentLength = self.ufragmentlen
        pls = numpy.sort(fragmentLength)
        N = len(fragmentLength)
        sums = []
        sizes = []            
        mysum = self.fragmentSum(None,strand)
    
        for i in numpy.arange(0,0.98,0.015):
            b1,b2 = pls[i*N],pls[(i+0.015)*N-1]
            p = (b1 < fragmentLength)* (b2 > fragmentLength)
            value = numpy.mean(mysum[p])            
            sums.append(value)
            sizes.append(numpy.mean(fragmentLength[p]))    
        if color == None: plt.plot(sizes,sums,'x-',markersize = 3,linewidth = 2,label = label)
        else:  plt.plot(sizes,sums,color,linewidth = 2,label = label)
    
        
    def plotScaling(self,fragments1 = None,fragments2 = None,color = None,label = "",weights = False ,normalize = True):
        "plots scaling over, possibly uses subset of fragmetns, or weigts, possibly normalizes after plotting"
        if fragments1 == None: fragments1 = self.ufragments
        if fragments2 == None: fragments2 = self.ufragments
        lens = numpy.array(numutils.logbins(10000,180000000,1.25),float)+0.1        
        positions = []
        values = []
        p11 = arrayInArray(self.fragments1,fragments1)
        p12 = arrayInArray(self.fragments1,fragments2)
        p21 = arrayInArray(self.fragments2,fragments1)
        p22 = arrayInArray(self.fragments2,fragments2)
        mask = numpy.logical_or(p11* p22,p12 * p21)
        mask2 = numpy.logical_and(mask,self.distances > 0)
        dind = numpy.argsort(self.distances[mask2])
        distances = self.distances[mask2][dind]        
        "calculating fragments lengths for exclusions to expected # of counts"
        args = numpy.argsort(self.ufragments)
        usort = self.ufragments[args]
        ulen = self.ufragmentlen[args]
        if weights == True:
            try:
                self.weights
            except:
                self.calculateWeights()
            uweights = self.weights[args]
            weights1 = uweights[numpy.searchsorted(usort,fragments1)]
            weights2 = uweights[numpy.searchsorted(usort,fragments2)]
        len1 = ulen[numpy.searchsorted(usort,fragments1)]
        len2 = ulen[numpy.searchsorted(usort,fragments2)]
        if (weights == True):
            maxcountsall = self.calculateExpectedCount(fragments1, fragments2, lens[:-1], lens[1:],len1,len2,weights1,weights2)
        else:
            maxcountsall = self.calculateExpectedCount(fragments1, fragments2, lens[:-1], lens[1:],len1,len2)
        for i in xrange(len(lens) - 1):
            beg = lens[i]
            end = lens[i+1]
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
        

    def calculateExpectedCount(self,fragments1,fragments2,  #fragments to exclude
                               lenmins,lenmaxs,  #beginnings and ends of bins
                               alllen1, alllen2,   #fragment lengthes
                               weights1=None,weights2=None
                               ):
        "Used for a fancy scaling calculation; is completely messed up"
        N = len(lenmins)
        count = [0]*N
        chr1 = fragments1 % 100 
        chr2 = fragments2 % 100 
        pos1 = fragments1 / 100        
        pos2 = fragments2 / 100
        chromosomes = set(chr1)
        #summedWeights = numpy.cumsum(weights2) 
        for chrom in chromosomes:    
            mask = (chr1 == chrom)
            mask2 = (chr2 == chrom)
            bp1 = pos1[mask]
            bp2 = pos2[mask2]
            p2arg = numpy.argsort(bp2)
            p2 = bp2[p2arg]
            if weights1 != None: 
                w1 = weights1[mask]
                w2 = weights2[mask2]
                w2sort = w2[p2arg]
                sw2 = numpy.r_[numpy.cumsum(w2sort),0]
            p1 = bp1
            
            for lenmin,lenmax,index in zip(lenmins,lenmaxs,range(N)):
                mymin = p1 - lenmax
                mymax = p1 - lenmin
                val1 = numpy.searchsorted(p2[:-1],mymin)
                val2 = numpy.searchsorted(p2[:-1],mymax)
                cval1 = numpy.maximum(val1-1,0)
                cval2 = numpy.maximum(val2-1,0)
                
                if weights1 == None: curcount = numpy.sum(numpy.abs(val1 - val2))
                else: curcount = numpy.sum(w1 * numpy.abs(sw2[cval1] - sw2[cval2]))
                mymin = p1 + lenmax
                mymax = p1 + lenmin
                val1 = numpy.searchsorted(p2[:-1],mymin)
                val2 = numpy.searchsorted(p2[:-1],mymax)
                cval1 = numpy.maximum(val1-1,-1)
                cval2 = numpy.maximum(val2-1,-1)

                if weights1 == None: curcount += numpy.sum(numpy.abs(val1 - val2))
                else: curcount += numpy.sum(w1 * numpy.abs(sw2[cval1] - sw2[cval2]))
                count[index] += curcount
        return count
    
    def plotRsiteStartDistribution(self,useSSReadsOnly = False,offset = 5,length = 200):
        if useSSReadsOnly == True:
            mask = self.SS
        else:
            mask = self.DS
        dists1 = self.lengthes1 - numpy.array(self.dists1,dtype = "int32") 
        dists2 = self.lengthes2 - numpy.array(self.dists2,dtype = "int32") 
        m = min(dists1.min(),dists2.min())
        if offset < -m: 
            offset = -m
            print "minimum negative distance is %d, larger than offset; offset set to %d" % (m,-m)
        dists1 += offset
        dists2 += offset          
        myrange = numpy.arange(-offset,length-offset)
        
        plt.subplot(141)
        plt.title("strand1, side 1")
        plt.plot(myrange,numpy.bincount(5+dists1[mask][self.strand1[mask] == True ])[:length])
        plt.subplot(142)
        plt.title("strand1, side 2")
        plt.plot(myrange,numpy.bincount(dists2[mask][self.strand1[mask] == True ])[:length])
        
        plt.subplot(143)
        plt.title("strand2, side 1")
        plt.plot(myrange,numpy.bincount(dists1[mask][self.strand1[mask] == False ])[:length])
        plt.subplot(144)
        plt.title("strand2, side 2")
        plt.plot(myrange,numpy.bincount(dists2[mask][self.strand1[mask] == False ])[:length])
        
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
    s1 = TR.fragmentSum(strand = 1)
    TR.saveFragments()
    TR.maskFilter(TR.dists1 > 500)
    TR.originalFragments()
    s2 = TR.fragmentSum(strand = 1)
    resolution = 1000000
    def coverage(s1,s2,TR):
        genome = dnaUtils.Genome("HG18")
        genome.createMapping(resolution)            
        label = genome.chromosomeStarts[TR.ufragments % 100 - 1] + (TR.ufragments / 100 ) / resolution
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
    s1 = TR.fragmentSum(strand = 1)
    TR.saveFragments()
    TR.maskFilter(TR.dists1 > 500)
    TR.originalFragments()
    s2 = TR.fragmentSum(strand = 1)
    resolution = 1000000
    plt.subplot(122)
    plt.title("NcoI")
    coverage(s1,s2,TR)
    plt.show() 

 





  
