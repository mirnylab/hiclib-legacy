"""
This is a module class for fragment-level Hi-C data analysis.
The base class "HiCdataset" can load, save and merge Hi-C datasets, perform certain filters, and save binned heatmaps.   

Additional class HiCStatistics contains methods to analyze HiC data on a fragment level. 
This includes read statistics, scalings, etc.


"""

import warnings
from mirnylab import systemutils
systemutils.setExceptionHook() 
import os,cPickle
from mirnylab.genome import Genome 
import mirnylab.hic.mapping 
import numpy
from numpy import array as na  
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt 
from mirnylab.h5dict import h5dict 
from mirnylab import plotting 
from mirnylab.plotting import mat_img,removeAxes

from mirnylab import numutils  
from mirnylab.numutils import arrayInArray,  sumByArray, correct, ultracorrect

r_ = numpy.r_

def corr(x,y): return stats.spearmanr(x, y)[0]        


class HiCdataset(object):
    """Base class to operate on HiC dataset. 
    
    This class stores all information about HiC reads on a hard drive. 
    Whenever a variable corresponding to any record is used, it is loaded/saved from/to the HDD. 
    
    If you apply any filters to a dataset, it will actually modify the content of the current working file.  
    Thus, to preserve the data, loading datasets is advised. """
    
    def __init__(self, filename , genome , maximumMoleculeLength = 500 , override = True , autoFlush = True):
        """        
        __init__ method 
        
        Initializes empty dataset by default. 
        If "override" is False, works with existing dataset.   
        
        Parameters
        ----------
        filename : string 
            A filename to store HiC dataset in an HDF5 file.  
        genome : folder with genome, or Genome object 
            A folder with fastq files of the genome and gap table from Genome browser.
            Alternatively, mirnylab.genome.Genome object.               
        maximumMoleculeLength : int, optional 
            Maximum length of molecules in the HiC library, used as a cutoff for dangling ends filter
        override : bool, optional
            Use specified dataset, do not remove it's contents. 
            By default, if filename exists, it is deleted upon initialization.
        autoFlush : bool, optional
            Set to True to disable autoflush - possibly speeds up read/write operations. 
            Don't forget to run flush then!             
        
        """                
        #---------->>> Important::: do not define any variables before vectors!!! <<<-------- 
        #These are fields that will be kept on a hard drive 
        self.vectors = {"chrms1":"int8","chrms2":"int8", #chromosomes. If >chromosomeCount, then it's second chromosome arm! 
                        "mids1":"int32","mids2":"int32",  #midpoint of a fragment, determined as "(start+end)/2"
                        "fraglens1":"int32","fraglens2":"int32", #fragment lengthes                        
                        "distances":"int32", #distance between fragments. If -1, different chromosomes. If -2, different arms.                         
                        "fragids1":"int64","fragids2":"int64",  #IDs of fragments. fragIDmult * chromosome + location                          
                        "dists1":"int32","dists2":"int32",        #distance to rsite
                        "cuts1":"int32","cuts2":"int32",           #precise location of cut-site 
                        "strands1":"bool","strands2":"bool",
                        "DS":"bool","SS":"bool"}
        
        self.autoFlush = autoFlush
                       
        if type(genome) == str: 
            self.genome = Genome(genomePath = genome, readChrms = ["#","X"])
        else:
            self.genome = genome             
        assert isinstance(self.genome, Genome)        
        
        self.chromosomeCount = self.genome.chrmCount  #used for building heatmaps
        self.fragIDmult = self.genome.fragIDmult
        print "----> New dataset opened, genome %s,  %s chromosomes" % (self.genome.folderName, self.chromosomeCount)

        self.maximumMoleculeLength = maximumMoleculeLength  #maximum length of a molecule for SS reads        
        self.filename = filename #File to save the data                             
                 
        if os.path.exists(self.filename):
            if override == False: 
                print "----->!!!File already exists! It will be opened in the 'append' mode."  
                print
            else:
                os.remove(self.filename)
        
        self.h5dict = h5dict(self.filename,autoflush = self.autoFlush )

    def _setData(self,name,data):
        "an internal method to save numpy arrays to HDD quickly"
        if name not in self.vectors.keys():
            raise ValueError("Attept to save data not specified in self.vectors")
        dtype = numpy.dtype(self.vectors[name])         
        data = numpy.asarray(data,dtype=dtype)        
        self.h5dict[name] = data

    
    def _getData(self,name):
        "an internal method to load numpy arrays from HDD quickly"         
        if name not in self.vectors.keys():
            raise ValueError("Attept to load data not specified in self.vectors")
        return self.h5dict[name]                
    
    def __getattribute__(self,x):
        "a method that overrides set/get operation for self.vectors so that they're always on HDD"
        if x == "vectors": return object.__getattribute__(self,x)
        
        if x in self.vectors.keys():            
            a =  self._getData(x)            
            return a                 
        else:
            return object.__getattribute__(self,x)

    def __setattr__(self,x,value):
        "a method that overrides set/get operation for self.vectors so that they're always on HDD"
        if x == "vectors": return object.__setattr__(self,x,value)
        
        if x in self.vectors.keys():        
            self._setData(x,value)                        
        else:
            return object.__setattr__(self,x,value)
        
    def flush(self):
        "Flushes h5dict if used in autoFlush = False mode"
        self.h5dict.flush()
        
    
    def merge(self,filenames):
        """combines data from multiple datasets
        
        Parameters
        ----------
            folders : list of strings
                List of folders to merge to current working folder                
        """
        if self.filename in filenames:
            self.exitProgram("----> Cannot merge folder into itself! Create a new folder")
        h5dicts = [h5dict(i,mode = 'r') for i in filenames]        
        for name in self.vectors.keys():
            res = []
            for mydict in h5dicts:
                res.append(mydict[name])
            res = numpy.concatenate(res)
            self._setData(name,res)

    
    def parseInputData(self,dictLike,zeroBaseChrom = True ,enzymeToFillRsites = None):        
        """Inputs data from a dictionary-like object, containing coordinates of the reads. 
        Performs filtering of the reads.   
        
        A good example of a dict-like object is a numpy.savez 
        
        ..warning::  Restriction fragments MUST be specified exactly as in the Genome class.
        
        ..warning:: Strand information is needed for proper scaling calculations, but will be imitated if not provided
        
        ..note::  
                    
        
        Parameters
        ----------
        dictLike["chrms1,2"] : Chromosomes of 2 sides of the read 
        dictLike["cuts1,2"] : Exact position of cuts
        dictLike["strands1,2"], essential : Direction of the read
        dictLike["rsites1,2"], optional : Position of rsite to which the read is pointing                         
        dictLike["uprsites1,2"] , optional : rsite upstream (larger genomic coordinate) of the cut position
        dictLike["downrsites1,2"] , optional  : rsite downstream (smaller genomic coordinate) of the cut position
        
        zeroBaseChrom : bool , optional
            Use zero-base chromosome counting if True, one-base if False
        enzymeToFillRsites = None or str, optional if rsites are specified
            Enzyme name to use with Bio.restriction
                 
        """
        
        rsite_related = ["rsites1","rsites2","uprsites1","uprsites2","downrsites1","downrsites2"]
        
        if False not in [i in dictLike.keys() for i in rsite_related]:
            noRsites = False            
        else:
            noRsites = True
        
        "Filling in chromosomes and positions - mandatory objects"
        a = dictLike["chrms1"]
        self.trackLen = len(a)         
        
        if zeroBaseChrom == True:
            self.chrms1 = a        
            self.chrms2 = dictLike["chrms2"]
        else:
            self.chrms1 = a - 1
            self.chrms2 = dictLike["chrms2"] - 1
            
        self.cuts1 = dictLike["cuts1"]
        self.cuts2 = dictLike["cuts2"]    
        

        if not (("strands1" in dictLike.keys()) and ("strands2" in dictLike.keys())):
            print "No strand information provided, assigning random strands."
            t = numpy.random.randint(0,2,self.trackLen)
            self.strands1 = t
            self.strands2 = 1-t
            noStrand = True 
        else:
            self.strands1 = dictLike["strands1"]
            self.strands2 = dictLike["strands2"]            
            noStrand = False   #strand information filled in 
        
        if noRsites == True:   #We have to fill rsites ousrlves. Let's see what enzyme to use! 
            if (enzymeToFillRsites == None) and (self.genome.hasEnzyme() == False) :
                raise ValueError("Please specify enzyme if your data has no rsites")
            if enzymeToFillRsites != None:
                if self.genome.hasEnzyme() == True:
                    if enzymeToFillRsites != self.genome.enzymeName:
                        warnings.warn("genome.enzymeName different from supplied enzyme")                                         
                self.genome.setEnzyme(enzymeToFillRsites)
            #enzymeToFillRsites has preference over self.genome's enzyme
                
            print "Filling rsites"            
            rsitedict = h5dict()  #creating dict to pass to anton's code 
            rsitedict["chrms1"] = self.chrms1
            rsitedict["chrms2"] = self.chrms2
            rsitedict["cuts1"] = self.cuts1
            rsitedict["cuts2"] = self.cuts2            
            rsitedict["strands1"] = self.strands1
            rsitedict["strands2"] = self.strands2            
            mirnylab.hic.mapping.fill_rsites(lib = rsitedict, genome_db = self.genome)
        else:
            rsitedict = dictLike #rsite information is in our dictionary        
        
        self.DS = (self.chrms1 >= 0) * (self.chrms2 >=0)   #if we have reads from both chromosomes, we're a DS read
        self.SS = (self.DS == False)
        
        self.dists1 = numpy.abs(rsitedict["rsites1"] - self.cuts1)
        self.dists2 = numpy.abs(rsitedict["rsites2"] - self.cuts2)
        
        self.mids1 = (rsitedict["uprsites1"] + rsitedict["downrsites1"])/2
        self.mids2 = (rsitedict["uprsites2"] + rsitedict["downrsites2"])/2
        
        self.fraglens1 = numpy.abs((rsitedict["uprsites1"] - rsitedict["downrsites1"]))
        self.fraglens2 = numpy.abs((rsitedict["uprsites2"] - rsitedict["downrsites2"]))
        
        del rsitedict   #deletes hdf5 file, so it's important 
        
        self.fragids1 = self.mids1 + numpy.array(self.chrms1,dtype = "int64") * self.fragIDmult
        self.fragids2 = self.mids2 + numpy.array(self.chrms2,dtype = "int64") * self.fragIDmult 
        
        distances = numpy.abs(self.mids1 - self.mids2)
        distances[self.chrms1 != self.chrms2] = -1
        self.distances = distances   #distances between restriction fragments
        
        mask = self.chrms1 == -1  #moving SS reads to the first side
        mask #Eclipse warning removal 
        variables = set([i[:-1] for i in self.vectors.keys() if i[-1] == "1"])  #set of variables to change
        for i in variables:  
            exec("a = self.%s1" % i)
            exec("b = self.%s2" % i)
            exec("a[mask] = b[mask]")
            exec("b[mask] = -1")
            exec("self.%s1 = a" % i)
            exec("self.%s2 = b" % i)
                
        mask = (self.fragids1 != self.fragids2)   #Discard dangling ends and self-circles                
        mask *=  ((self.chrms1 < self.chromosomeCount) * (self.chrms2 < self.chromosomeCount))  #Discard unused chromosomes                
        mask *= ((self.chrms1 >=0) + (self.chrms2 >=0))   #Has to have at least one side mapped                
        if noStrand == True: 
            dist = numpy.abs(self.cuts1 - self.cuts2)  #Can't tell if reads point to each other. 
        else: 
            dist = - self.cuts1 * (2 * self.strands1 -1) - self.cuts2 * (2 * self.strands2 - 1)  #distance between sites facing each other
                                
                                
        readsMolecules = (self.chrms1 == self.chrms2)*(self.strands1 != self.strands2) *  (dist >=0) * (dist <= self.maximumMoleculeLength)#filtering out DE         
        mask *= (readsMolecules  == False)                                
        
        
                        
        print len(mask),mask.sum()
        self.maskFilter(mask)
        
            
    def saveFragments(self):
        "saves fragment data to make correct expected estimates after applying a heavy mask"
        self.ufragmentsOriginal = numpy.array(self.ufragments)
        self.ufragmentlenOriginal = numpy.array(self.ufragmentlen)

    def originalFragments(self):
        "loads original fragments"
        self.ufragments = numpy.array(self.ufragmentsOriginal)
        self.ufragmentlen = numpy.array(self.ufragmentlenOriginal)

    def calculateWeights(self):
        """Calculates weights for reads based on fragment length correction similar to Tanay's;
         may be used for scalings or creating heatmaps"""        
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



    def buildHeatmap(self,chromosome = 14,chromosome2 = None,resolution = 1000000,show = False,useRsiteDensity = True):
        """Builds heatmaps between any two chromosomes at a given resolution
        
        Parameters
        ----------
        chromosome1 : int
            First chromosome
        chromosome2 : int, optional 
            Second chromosome, if the map is trans. 
        resolution : int
            Resolution of a heatmap. Default is 1M  
        show : bool , optional
            Show heatmap, or just output it?
        useRsiteDensity : bool, optional
            Correct map by density of rsites.              
        
        """        
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

        if chromosome2 == None: chromosome2 = chromosome
        if useRsiteDensity == True:            
            m1 = self.ufragments / self.fragIDmult == chromosome
            m2 = self.ufragments / self.fragIDmult == chromosome2
            p1 = self.ufragments[m1]  % self.fragIDmult 
            p2 = self.ufragments[m2]  % self.fragIDmult
            p1mod = p1 / resolution
            p2mod = p2 / resolution            
            myarray = numpy.array(range(numpy.max(p1mod)+1))
            vec1 =  sumByArray(p1mod,myarray)
            myarray = numpy.array(range(numpy.max(p2mod)+1))
            vec2 = sumByArray(p2mod,myarray)
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
        """Creates an all-by-all heatmap in accordance with mapping provided by 'genome' class
        
        Parameters
        ----------
        Resolution : int
            Resolution of a heatmap 
        """ 
        self.genome.setResolution(resolution)
        dr = self.DS 
        
        label1 = self.genome.chrmStartsBinCont[self.chrms1[dr] ] + self.mids1[dr] / resolution
        label1 = numpy.array(label1, dtype = "uint32")
        label2 = self.genome.chrmStartsBinCont[self.chrms2[dr] ] + self.mids2[dr] / resolution
        label2 = numpy.array(label2, dtype = "uint32")       
        label = label1 * numpy.int64(self.genome.numBins) + label2
        del label1
        del label2
        if self.genome.numBins < 65000:
            label = numpy.array(label, dtype = "uint32")                    
         
        counts = numpy.bincount(label, minlength = self.genome.numBins**2)
        if len(counts) > self.genome.numBins**2:
            print "heatmap exceed length of the genome!!! Check genome"
            exit()
            
        counts.shape = (self.genome.numBins,self.genome.numBins)
        for i in xrange(len(counts)):
            counts[i,i:] += counts[i:,i]
            counts[i:,i] = counts[i,i:]        
        return counts 
    
    def buildSinglesCoverage(self,resolution):
        "creates an SS coverage vector heatmap in accordance with the output of the 'genome' class"
        self.genome.setResolution(resolution)
        ds = self.DS == False        
        label = self.genome.chrmStartsBinCont[self.chrms1[ds] ] + self.mids1[ds] / resolution
        counts = sumByArray(label, numpy.arange(self.genome.numBins))
        return counts
    
    def buildFragmetCoverage(self,resolution):
        "creates restriction site density vector (visible sites only) in accordance with the 'genome' class"         
        self.genome.setResolution(resolution)
        try: self.ufragments
        except: self.rebuildFragments()
        chroms = self.ufragments / self.fragIDmult
        positions = self.ufragments % self.fragIDmult
        label = self.genome.chrmStartsBinCont[chroms - 1] + positions / resolution
        counts = sumByArray(label, numpy.arange(self.genome.numBins))
        return counts
        

    
        
    
    def fragmentFilter(self,fragments):
        """keeps only reads that originate from fragments in 'fragments' variable, for DS - on both sides
        Parameters 
        ----------
        fragments : np.array of fragment IDs or bools
            List of fragments to keep, or their indexes in self.ufragments        
        """
        if fragments.dtype == numpy.bool:
            fragments = self.ufragments[fragments]        
        m1 = arrayInArray(self.fragids1,fragments)
        m2 = arrayInArray(self.fragids2,fragments) + self.SS
        mask = numpy.logical_and(m1,m2)
        self.maskFilter(mask)

    
    def maskFilter(self,mask):
        """keeps only reads designated by mask
        Parameters
        ----------
        mask : array of bools
            Indexes of reads to keep 
        """
        
        print "          Number of reads changed  %d ---> %d" % (len(mask),mask.sum()),
        length = 0 
        for name in self.vectors:
            data = self._getData(name)
            ld = len(data)
            if length == 0: 
                length = ld
                
            else:
                if ld != length: 
                    self.delete()             
            self._setData(name,data[mask])
        self.N = mask.sum()   
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
        
            

    def filterExtreme(self,cutH = 0.005,cutL = 0):
        """removes fragments with most and/or least # counts
        
        Parameters
        ----------
        cutH : float, 0<=cutH < 1, optional
            Fraction of the most-counts fragments to be removed
        cutL: float, 0<=cutL<1, optional 
            Fraction of the least-counts fragments to be removed
        """
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
        """removes very large and small fragments
        
        Parameters
        ----------
        cutlarge : int 
            remove fragments larger than it
        cutsmall : int
            remove fragments smaller than it
        """
        print "----->Small/large fragments filter: keep strictly less than %d, strictly more than %d bp" % (cutlarge, cutsmall)                
        p = (self.ufragmentlen < (cutlarge ) ) * (self.ufragmentlen > cutsmall)                      
        self.fragmentFilter(self.ufragments[p])
        print 

    
    def filterRsiteStart(self,offset = 5):
        """removes reads that start within x bp near rsite
        
        Parameters
        ----------
        offset : int 
            Number of bp to exclude next to rsite, not including offset
        """
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
        
        dups = numpy.zeros((self.N,2),dtype = "int64",order = "C")
        dups[:,0] = numpy.array(self.cuts1 , dtype = "int64") + self.chrms1 * self.fragIDmult 
        dups[:,1] = numpy.array(self.cuts1 , dtype = "int64") + self.chrms1 * self.fragIDmult   
        dups.shape = (self.N * 2)
        strings = dups.view("|S16")
        assert len(strings) == self.N
        uids = numpy.unique(strings,return_index = True)[1]
        del strings, dups 
        stay = self.SS.copy()
        stay[uids] = True
        ds = self.DS.sum()
        print "     Number of DS reads changed - %d ---> %d" % (ds,ds - len(self.DS) + stay.sum()) 
        del uids
        self.maskFilter(stay)
        print
        

    def fragmentSum(self,fragments = None, strands = "both"):
        """returns sum of all counts for a set or subset of fragments
        Parameters
        ---------- 
        fragments : list of fragment IDs, optional
            Use only this fragments. By default all fragments are used
        strands : 1,2 or "both" (default) 
            Use only first or second side of the read (first has SS, second - doesn't) 
        """
        if fragments == None: fragments = self.ufragments                
        if strands == "both":  
            return sumByArray(self.fragids1,fragments) + sumByArray(self.fragids2[self.DS],fragments) 
        if strands == 1: return sumByArray(self.fragids1,fragments)
        if strands == 2:return sumByArray(self.fragids2[self.DS],fragments)
        
        
    def printStats(self):
        print "-----> Statistics for the file  %s!" % self.filename
        print "     Single sided reads: " ,self.SS.sum()
        print "     Double sided reads: " , self.DS.sum()
        ss1 = self.strands1[self.chrms1>=0]
        ss2 = self.strands2[self.chrms2>=0]
        sf = ss1.sum() + ss2.sum()
        sr = len(ss1) + len(ss2) - sf
        print "     reverse/forward bias",float(sr)/sf
                    
        
    def save(self,filename):
        "Saves dataset to filename, does not change the working file."
        if self.filename == filename: self.exitProgram("Cannot save to the working file")
        newh5dict = h5dict(filename,mode = 'w')
        for name in self.vectors.keys(): newh5dict[name] = self.h5dict[name]
        print "----> Data saved to file %s" % (filename,)
    
    def load(self,filename):
        "Loads dataset from file to working file; check for inconsistency"
        otherh5dict = h5dict(filename,'r')
        length = 0 
        for name in self.vectors:
            data = otherh5dict[name]
            ld = len(data)
            if length == 0: 
                length = ld
            else:
                if ld != length: 
                    print("---->!!!!!File %s contains inconsistend data<----" % filename)
                    self.exitProgram("----> Sorry...")
                    
            self._setData(name,data) 
        print "---->Loaded data from file %s, contains %d reads" % (filename, length)
            

        
    def saveHeatmap(self,filename,resolution = 1000000):
        """
        Saves heatmap to filename at given resolution. 
        
        .. note:: More than one file may be used to store the data, as joblib creates extra files
        to store huge numpy arrays. All files will start with filename, but some may end with _01.npy.z  
        
        Parameters
        ----------
        filename : str
        resolution : int              
        """
        
        try: os.remove(filename)
        except: pass        
        heatmap = self.buildAllHeatmap(resolution)        
        singles = self.buildSinglesCoverage(resolution)        
        frags = self.buildFragmetCoverage(resolution)
        chromosomeStarts = numpy.array(self.genome.chrmStartsBinCont)
        tosave = h5dict(path = filename,mode = "w")
        tosave["resolution"] = resolution
        tosave["heatmap"] = heatmap
        tosave["singles"] = singles
        tosave["frags"] = frags
        tosave["genomeBinNum"] = self.genome.numBins
        tosave["chromosomeSTarts"] = chromosomeStarts
        

        
        print "----> Heatmap saved to '%s' at %d resolution" % (filename,resolution)
        
    def exitProgram(self,a):
        print a
        print "     ----> Bye! :) <----"
        exit()

class HiCStatistics(HiCdataset):
    "a sub-class of a 'HiCdataset' class used to do statistics on Hi-C reads" 

    def multiplyForTesting(self,N = 100):
        "used for heavy load testing only"
        for name in self.vectors:
            print "multipliing",name
            one = self._getData(name)
            if name in ["mids1","mids2"]:
                one += numpy.random.randint(0,5*N,len(one))            
            blowup = numpy.hstack(tuple([one]*N))
            self._setData(name,blowup)
            
    
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
                    useWeights = False ,        #use weights associated with fragment length                    
                    excludeNeighbors = None, enzyme = None,   #number of neighboring fragments to exclude. Enzyme is needed for that!   
                    normalize = True,          #normalize the final plot to sum to one
                    withinArms = True,                #Treat chromosomal arms separately
                    mindist = 10000,  #Scaling was proved to be unreliable under 10000 bp for 6-cutter enzymes
                    maxdist = None,
                    #----------Calculating scaling within a set of regions only------
                    regions = None,   #Array of tuples (chrom, start, end) for which scaling should be calculated
                    #Note that calculation might be extremely long (it might be proportional to # of regions for # > 100)
                    appendReadCount = True,   #Append read count to the plot label 
                    **kwargs #kwargs to be passed to plotting                       
                    ):               #Sad smiley, because this method is very painful and complicated
        """plots scaling over, possibly uses subset of fragmetns, or weigts, possibly normalizes after plotting
        
        Plan of scaling calculation: 
        
        1. Subdivide all genome into regions. \n         
            a. Different chromosomes \n            
            b. Different arms \n            
            c. User defined squares/rectangles on a contact map \n            
               -(chromosome, start,end) square around the diagonal \n                
               -(chr1, start1, end1, start2, end2) rectangle \n
                       
        2. Use either all fragments, or only interactions between two groups of fragments \n         
            e.g. you can calculate how scaling for small fragments is different from that for large \n            
            It can be possibly used for testing Hi-C protocol issues. \n             
            One can see effect of weights by doing this \n        
            
        3. (optional) Calculate correction associated with fragment length dependence 
                         
        4. Subdivide all possible genomic separation into log-spaced bins
                 
        5. Calculate expected number of fragment pairs within each bin (possibly with weights from step 3).
         
            If exclusion of neighbors is specificed, expected number of fragments knows about this
            
        Parameters
        ----------
        fragids1, fragids2 : numpy.array of fragment IDs, optional 
            Scaling is calculated only for interactions between fragids1 and fragids2
            If omitted, all fragments are used
            If boolean array is supplied, it serves as a mask for fragments. 
        useWeights : bool, optional 
            Use weights calculated from fragment length
        excludeNeighbors : int or None, optional
            If None, all fragment pairs are considered. 
            If integer, only fragment pairs separated by this rsites are considered.
        enzyme : string ("HindIII","NcoI")
            If excludeNeighbors is used, you have to specify restriction enzyme
        normalize : bool, optional
            Do an overall normalization of the answer, by default True. 
        withinArms : bool, optional
            Set to false to use whole chromosomes instead of arms
        mindist, maxdist : int, optional
            Use lengthes from mindist to maxdist
        regions : list of (chrom, start,end), optional
            Restrict scaling calculation to only certain squares of the map
        appendReadCount : bool, optional 
            Append read count to the plot label
        kwargs : dict, optional
            Dictionary of kw args to be passed to plt.plot
            
        Returns
        -------
        (bins,probabilities) - values to plot on the scaling plot        
         
        """
        if excludeNeighbors <= 0: excludeNeighbors = None   #Not excluding neighbors  
        #use all fragments if they're not specified 
        if fragids1 == None: fragids1 = self.ufragments
        if fragids2 == None: fragids2 = self.ufragments
        if fragids1.dtype == numpy.bool: 
            fragids1 = self.ufragments[fragids1]
        if fragids2.dtype == numpy.bool:
            fragids2 = self.ufragments[fragids2]
        
        #Calculate regions if not specified 
        if regions == None: 
            if withinArms == False: 
                regions = [(i,0,self.genome.chrmLens[i]) for i in xrange(self.genome.chrmCount)]
            else:
                regions = [(i,0,self.genome.cntrMids[i]) for i in xrange(self.genome.chrmCount)] + [(i,self.genome.cntrMids[i],self.genome.chrmLens[i]) for i in xrange(self.genome.chrmCount)]
                
        
        if maxdist == None: maxdist = max ( max([i[2] - i[1] for i in regions]),  #normal regions 
                        max([abs(i[2] - i[3]) for i in regions if len(i) > 3] + [0]),   #rectangular regions
                        max([abs(i[1] - i[4]) for i in regions if len(i) > 3] + [0]))

        regionID = numpy.zeros(len(self.chrms1),numpy.int16) - 1  #Region to which a read belongs             
        chr1 = self.chrms1
        chr2 = self.chrms2
        pos1 = self.mids1
        pos2 = self.mids2
        fragRegions1 = numpy.zeros(len(fragids1),int) - 1
        fragRegions2 = numpy.zeros(len(fragids2),int) - 1 
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
        print regions                         

                                        
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
            for lenmin,lenmax,index in zip(lenmins,lenmaxs,range(len(lenmins))):
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
                #print index,count[index]
        maxcountsall = count        
        rawValues = []
        for i in xrange(len(lens) - 1):   #Dividing observed by expected
            beg,end  = lens[i], lens[i+1]            
            first,last = tuple(numpy.searchsorted(distances,[beg,end]))            
            mycounts = last - first
            maxcounts = maxcountsall[i]                                     
            positions.append(sqrt(float(beg)*float(end)))
            values.append(mycounts/float(maxcounts))
            rawValues.append(mycounts)
        print rawValues
        print count
            
        positions = numpy.array(positions)
        values = numpy.array(values)
         
        if normalize == True: values /= numpy.sum(1. * (positions * values) [numpy.logical_not(numpy.isnan(positions * values))])        
        
        if appendReadCount == True: 
            if "label" in kwargs.keys(): 
                kwargs["label"]  = kwargs["label"] + ", %d reads" % len(distances)
        plt.plot(positions,values, **kwargs)
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






class experimentalFeatures(HiCdataset):
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



def plotFigure2c():
    TR = HiCdataset()
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
    TR = HiCdataset()
    TR.load("GM-all.refined")
    s1 = TR.fragmentSum(strands = 1)
    TR.saveFragments()
    TR.maskFilter(TR.dists1 > 500)
    TR.originalFragments()
    s2 = TR.fragmentSum(strands = 1)
    resolution = 1000000
    def coverage(s1,s2,TR):
        genome = Genome()
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
    
    TR = HiCdataset()
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

 





  