"""
Binned data - analysis of HiC, binned to resolution.

Concepts
--------

class Binned Data allows low-level manipulation of  multiple HiC datasets, 
binned to the same resolution from the same genome. 

When working with multiple datasets, all the filters will be synchronized,
so only bins present in all datasets will be considered for the analysis. 
Removal of bins from one dataset will remove them from the others. 
E.g. removing 1% of bins with lowest # of count might remove more than 1% of total bins, 
when working with 2 or more datasets.  

Class has significant knowledge about filters that have been applied. 
If an essential filter was not applied, it will throw an exception; 
if adviced filter is not applied, it will throw a warning. 
Most of the methods have an optional "force" argument that will ignore dependensies.
I
  

We provide example scripts that show ideal protocols for certain types of the analysis, 
but they don't cover the realm of all possible manipulations that can be performed with this class.

Input data
----------

method :py:func:`SimpleLoad <binnedData.simpleLoad>` may be used to load the data. 
It automatically checks for possible genome length mismatch. 

This method works best with h5dict files, created by fragmentHiC.
In this case you just need to supply the filename.  

It can also accept any dictionary-like object with the following keys, 
where all but "heatmap" is optional.  

* ["heatmap"] : all-by-all heatmap

* ["singles"] : vector of SS reads, optional 

* ["frags"] : number of rsites per bin, optional

* ["resolution"] : resolution 


Variables
---------

self.dataDict - dictionary with heatmaps; keys are provided when loading the data.

self.singlesDict - dictionary with SS read vectors. Keys are the same. 

self.trackDict - dictionary with genomic tracks, such as GC content. 
Custom tracks may be also added to this dictionary.

--------------------------------------------------------------- 
"""

import os 
from mirnylab import systemutils,numutils

import mirnylab.plotting
from mirnylab.plotting import  removeBorder 
from mirnylab.numutils import PCA, EIG,correct, ultracorrectSymmetricWithVector
from mirnylab.genome import Genome 
import  numpy
from math import exp
from mirnylab.h5dict import h5dict  
from scipy import weave
import warnings 
from scipy.stats.stats import spearmanr
import matplotlib.pyplot as plt 


    
class binnedData(object):
    """Base class to work with binned data, the most documented and robust part of the code. 
    Further classes for other analysis are inherited from this class. 
    """
    
    def __init__(self, resolution,genome):
        """
        
        self.__init__ - initializes an empty dataset. 
        
        This method sets up a Genome object and resolution. 
        
        Genome object specifies genome version and inclusion/exclusion of sex chromosomes.
        
        Parameters
        ----------
        resolution : int 
            Resolution of all datasets
        genome : genome Folder or Genome object
        
        """        
        if type(genome) == str: 
            self.genome = Genome(genomePath = genome, readChrms = ["#","X"])
        else:
            self.genome = genome 
            
        assert isinstance(self.genome, Genome)

        if resolution != None: self.resolution = resolution
        self.chromosomes = self.genome.chrmLens
        self.resolution = resolution                    
        self.genome.setResolution(self.resolution)                        
        self._initChromosomes()        
        self.dataDict = {}
        self.biasDict = {}
        self.trackDict = {}
        self.singlesDict = {}
        self.fragsDict = {}
        self.PCDict = {}
        self.EigDict = {}
        self.dicts = [self.trackDict, self.biasDict, self.singlesDict, self.fragsDict]
        self.eigDicts = [self.PCDict, self.EigDict]
        
        
        self.appliedOperations = {}
        
    def _initChromosomes(self):
        "internal: loads mappings from the genome class based on resolution"        
        self.chromosomeStarts = self.genome.chrmStartsBinCont
        self.centromerePositions = self.genome.cntrMidsBinCont
        self.chromosomeEnds = self.genome.chrmEndsBinCont
        self.trackLength = self.genome.numBins
                
        self.chromosomeCount = self.genome.chrmCount
        self.chromosomeIndex = self.genome.chrmIdxBinCont
        self.positionIndex = self.genome.posBinCont        
        self.armIndex = self.chromosomeIndex * 2 + numpy.array(self.positionIndex > self.genome.cntrMids[self.chromosomeIndex],int)

    def _giveMask(self):
        "Returns index of all bins with non-zero read counts"
        self.mask = numpy.ones(len(self.dataDict.values()[0]),numpy.bool)
        for data in self.dataDict.values():
            datasum = numpy.sum(data,axis = 0)
            datamask = datasum > 0
            self.mask *= datamask
        return self.mask
    
    def _giveMask2D(self):
        "Returns outer product of _giveMask with itself, i.e. bins with possibly non-zero counts"
        self._giveMask()
        self.mask2D = self.mask[:,None] * self.mask[None,:]
        return self.mask2D   
                
    
    def simpleLoad(self,in_data,name):
        """Loads data from h5dict file or dict-like object
        
        Parameters
        ----------
        
        in_data : str or dict-like
            h5dict filename or dictionary-like object with input data            
        name : str
            Key of the dataset in self.dataDict
        
        """        
        if type(in_data) == str:
            if os.path.exists(in_data) == False: raise IOError("HDF5 dict do not exist, %s" % in_data) 
            alldata = h5dict(in_data,mode = "r")
        else:
            alldata = in_data
                             
        self.dataDict[name] = alldata["heatmap"]        
        try: self.singlesDict[name] = alldata["singles"]
        except: print "No SS reads found"
        try: self.fragsDict[name] = alldata["frags"]
        except:pass 
        
        if self.resolution != alldata["resolution"]:
            print "resolution mismatch!!!"
            print "--------------> Bye <-------------"
            raise StandardError("Resolution mismatch! ") 
        
        if self.genome.numBins != len(alldata["heatmap"]):              
            print "Genome length mismatch!!!"
            print "source genome",len(alldata["heatmap"])
            print "our genome",self.genome.numBins
            raise StandardError("Genome size mismatch! ")

              
        
        
    def loadGC(self):        
        "loads GC content at given resolution"
        self.trackDict["GC"] = numpy.concatenate(self.genome.GCBin)
  
    
    def removeDiagonal(self,m=1):
        """Removes all bins on a diagonal, and bins that are up to m away from the diagonal, including m. 
        By default, removes all bins touching the diagonal.
        
        Parameters
        ----------
        m : int, optional 
            Number of bins to remove
        """
        for i in self.dataDict.keys():
            data = self.dataDict[i] * 1.             
            N = len(data)
            N   #Eclipse warning remover 
            code = r"""
            #line 841 "binary_search.py"
            using namespace std; 
            for (int i = 0; i < N; i++)    
            {    
                for (int j = max(i-m,0); j<min(i+m+1,N); j++)
                {
                    data[i*N + j] = 0;
                }
            } 
            """
            support = """
            #include <math.h>
            """
            weave.inline(code, ['m','data',"N"], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
            self.dataDict[i] = data
        self.appliedOperations["RemovedDiagonal"] = True        
        
        
    
    def removeStandalone(self,offset = 3):
        """removes standalone bins (groups of less-than-offset bins)
        
        Parameters
        ----------
        offset : int 
            Maximum length of group of bins to be removed
        """                
        diffs = numpy.diff(numpy.array(numpy.r_[False, self._giveMask(), False],int))
        begins = numpy.nonzero(diffs == 1)[0] 
        ends = numpy.nonzero(diffs == -1)[0]
        beginsmask = (ends - begins) <= offset
        newbegins = begins[beginsmask]
        newends = ends[beginsmask]        
        print "removing %d standalone bins"% numpy.sum(newends - newbegins)
        mask = self._giveMask()
        for i in xrange(len(newbegins)): mask[newbegins[i]:newends[i]] = False
        mask2D = mask[:,None] * mask[None,:]        
        for i in self.dataDict.values(): i[mask2D == False] = 0
        self.appliedOperations["RemovedStandalone"] = True
        
    def removePoorRegions(self,names = None, cutoff = 2):
        """removes cutoff persent of bins with least counts
        
        Parameters
        ----------
        names : list of str
            List of datasets to perform the filter. All by default. 
        cutoff : int, 0<cutoff<100
            Percent of lowest-counts bins to be removed         
        """
        statmask = numpy.zeros(len(self.dataDict.values()[0]),numpy.bool)
        mask = numpy.ones(len(self.dataDict.values()[0]),numpy.bool)
        if names == None: names =self.dataDict.keys()  
        for i in names:
            data = self.dataDict[i]
            datasum = numpy.sum(data,axis = 0)            
            datamask = datasum > 0
            mask *= datamask  
            try: countsum = numpy.sum(data,axis = 0) + self.singlesDict[i]
            except: countsum = numpy.sum(data,axis = 0) 
            newmask = countsum >= numpy.percentile(countsum[datamask],cutoff)
            mask *= newmask  
            statmask [(newmask == False) * (datamask == True)] = True
        print "removed %d poor bins", statmask.sum() 
        mask2D = mask[:,None] * mask[None,:]
        for i in self.dataDict.values(): i[mask2D == False] = 0
        self.appliedOperations["RemovedPoor"] = True
              
    def truncTrans(self,high = 0.0005):
        """Truncates trans contacts to remove blowouts
        
        Parameters
        ----------
        high : float, 0<high<1, optional 
            Fraction of top trans interactions to be removed
        """
        for i in self.dataDict.keys():
            data = self.dataDict[i]
            transmask = self.chromosomeIndex[:,None] != self.chromosomeIndex[None,:]
            lim = numpy.percentile(data[transmask],100*(1 - high))
            tdata = data[transmask]
            tdata[tdata > lim] = lim            
            self.dataDict[i][transmask] = tdata
        self.appliedOperations["TruncedTrans"] = True                        
        
    def removeCis(self):
        "sets to zero all cis contacts"
        
        mask = self.chromosomeIndex[:,None] == self.chromosomeIndex[None,:]         
        for i in self.dataDict.keys():                
            self.dataDict[i][mask] = 0   
        self.appliedOperations["RemovedCis"] = True           
        print("All cis counts set to zero")
            
    def fakeCisOnce(self,mask = "CisCounts", silent = False):
        """Used to fake cis counts.
        If extra mask is supplied, it is used instead of cis counts. 
        
        Parameters
        ----------
        mask : NxN boolean array or "CisCounts" 
            Mask of elements to be faked. 
            If set to "CisCounts", cis counts will be faked
        silent : bool
            Do not print anything
             
        
        """
        if silent == False: print("All cis counts are substituted with matching trans count")
        for i in self.dataDict.keys():             
            data = numpy.asarray(self.dataDict[i],order = "C",dtype = float)
            if mask == "CisCounts": mask =  numpy.array(self.chromosomeIndex[:,None] == self.chromosomeIndex[None,:],int)
            else: assert mask.shape == self.dataDict.values()[0].shape  #check that mask has correct shape                              
            s = numpy.abs(numpy.sum(data,axis = 0)) <= 1e-20 
            mask[:,s]= 2
            mask[s,:] = 2              
            N = len(data)
            N
            code = r"""
            #line 310 "binnedData.py"
            using namespace std; 
            for (int i = 0; i < N; i++)    
            {    
                for (int j = i; j<N; j++)
                {
                    if (mask[i* N + j] == 1)                    
                    {
                    while (true) 
                        {
                        int r = rand() % 2;                         
                        if (r == 0)
                            {                            
                            int s = rand() % N;
                            if (mask[i * N + s] == 0)
                                {                                
                                data[i * N + j] = data[i * N + s];
                                data[j * N + i] = data[i * N + s];
                                break;
                                 
                                }
                            }
                        else
                            {
                            int s = rand() % N;
                            if (mask[s * N + j] == 0)
                                {
                                data[i * N + j] = data[i * N + s];
                                data[j * N + i] = data[i * N + s]; 
                                break;
                                }                            
                            }                        
                        }
                    }
                }
            } 
            """
            support = """
            #include <math.h>
            """
            weave.inline(code, ['mask','data',"N"], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
            self.dataDict[i] = data
            self.appliedOperations["RemovedCis"] = True 
            self.appliedOperations["FakedCis"] = True 
            #mat_img(data)

    def fakeCis(self):
        """This method fakes cis contacts in an interative way
        It is done to achieve faking cis contacts that is independent of normalization of the data. 
        """
        self.removeCis()
        self.iterativeCorrectWithoutSS(M=5)
        self.fakeCisOnce(silent = True)
        self.iterativeCorrectWithoutSS(M=5)
        self.fakeCisOnce(silent = True)
        self.iterativeCorrectWithoutSS(M=10) 
        print("All cis counts are substituted with faked counts")
        print("Data is iteratively corrected as a part of faking cis counts")


    def correct(self,names=None):
        """performs single correction without SS
        
        Parameters
        ----------
        names : list of str or None
            Keys of datasets to be corrected. If none, all are corrected. 
        """
        self.iterativeCorrectWithoutSS(names, M=1)         
        if ("RemovedDiagonal") not in self.appliedOperations.keys():        
            warnings.warn("Did you forget to remove diagonal?")
        
        self.appliedOperations["Corrected"] = True
    def iterativeCorrectWithoutSS(self, names = None,M=50):
        """performs iterative correction without SS
        
        Parameters
        ----------
        names : list of str or None, optional
            Keys of datasets to be corrected. By default, all are corrected.
        M : int, optional
            Number of iterations to perform. 
        """
        if names == None: names = self.dataDict.keys()
        
        for i in names:
            self.dataDict[i] = ultracorrectSymmetricWithVector(self.dataDict[i],M=M)[0]         
        if ("RemovedDiagonal") not in self.appliedOperations.keys():
            warnings.warn("Did you forget to remove diagonal?")
            
        self.appliedOperations["Corrected"] = True

    def iterativeCorrectWithSS(self,names = None,M = 55,force = False):
        """performs iterative correction with SS
        
        Parameters
        ----------
        names : list of str or None, optional
            Keys of datasets to be corrected. By default, all are corrected.
        M : int, optional
            Number of iterations to perform. 
        force : bool, optional 
            Force current operation 
        """
        if ("Corrected" in self.appliedOperations.keys()) and (force == False):
            raise StandardError("Cannot correct after previous correction was applied")
        if ("RemovedCis" in self.appliedOperations.keys()) and (force == False):
            raise StandardError("Cannot correct with SS if there are no cis reads")        
        
        if names == None: names = self.dataDict.keys()
        for i in names:
            data = self.dataDict[i]
            vec = self.singlesDict[i]
            ndata,nvec = ultracorrectSymmetricWithVector(data, vec,M=M)                         
            self.dataDict[i] = ndata
            self.singlesDict[i] = nvec
            self.biasDict[i] = (vec / nvec)
        if ("RemovedDiagonal") not in self.appliedOperations.keys():
            warnings.warn("Did you forget to remove diagonal?")
        self.appliedOperations["Corrected"] = True
        
    def removeChromosome(self,chromNum):
        "removes certain chromosome from all tracks and heatmaps, setting all values to zero "
        beg = self.genome.chrmStartsBinCont[chromNum]
        end = self.genome.chrmEndsBinCont[chromNum]
        for i in self.dataDict.values():
            i[beg:end] = 0
            i[:,beg:end] = 0 
        
        for mydict in self.dicts:
            for value in mydict.values():
                value[beg:end] = 0 
        
        for mydict in self.eigDicts:
            for value in mydict.values():
                value[beg:end] = 0 
                


    def removeZeros(self,zerosMask = None):
        """removes bins with zero counts
        keeps chromosome starts, ends, etc. consistent"""
        if zerosMask != None:
            s = zerosMask
        else:          
            s = numpy.sum(self._giveMask2D(),axis = 0) > 0
            for i in self.dataDict.values():
                s *= (numpy.sum(i,axis = 0) > 0)
        indices = numpy.zeros(len(s),int)
        count = 0 
        for i in xrange(len(indices)):
            if s[i] == True:
                indices[i] = count
                count +=1
            else: 
                indices[i] = count
        indices = numpy.r_[indices,indices[-1] + 1]  
        for i in self.dataDict.keys():
            a = self.dataDict[i]
            b = a[:,s]
            c = b[s,:]
            self.dataDict[i] = c            
        
        for mydict in self.dicts:
            for key in mydict.keys():
                mydict[key] = mydict[key][s]

        
        for mydict in self.eigDicts:
            for key in mydict.keys():
                
                mydict[key] = mydict[key][:,s]  
        
        self.chromosomeIndex = self.chromosomeIndex[s]
        self.armIndex = self.armIndex[s]
        self.chromosomeEnds = indices[self.chromosomeEnds]
        self.chromosomeStarts = indices[self.chromosomeStarts]
        self.centromerePositions = indices[self.centromerePositions]
        self.removeZerosMask = s
        if self.appliedOperations.get("RemovedZeros",False) == True:
            warnings.warn("You're removing zeros twice. You can't restore zeros now!")
        self.appliedOperations["RemovedZeros"] = True  
        return s 

    def restoreZeros(self, value = numpy.NAN):
        """Restores zeros that were removed by removeZeros command. 
        
        .. warning:: You can restore zeros only if you used removeZeros once.
        
        Parameters
        ----------
        value : number-like, optional. 
            Value to fill in missing regions. By default, NAN. 
        """
        if not hasattr(self,"removeZerosMask"): raise StandardError("Zeros have not been removed!")        
        
        s = self.removeZerosMask
        N = len(s)

        for i in self.dataDict.keys():
            a = self.dataDict[i]
            self.dataDict[i] = numpy.zeros((N,N),dtype = a.dtype) * value 
            tmp = numpy.zeros((N,len(a)),dtype = a.dtype) * value 
            tmp[s,:] = a
            self.dataDict[i][:,s] = tmp                     
        for mydict in self.dicts:
            for key in mydict.keys():
                a = mydict[key]
                mydict[key] = numpy.zeros(N,dtype = a.dtype) * value 
                mydict[key][s] = a
                        
        for mydict in self.eigDicts:
            for key in mydict.keys():
                a = mydict[key]
                mydict[key] = numpy.zeros((len(a),N),dtype = a.dtype) * value 
                mydict[key][:,s] = a
                                        
        self._initChromosomes()     
        self.appliedOperations["RemovedZeros"] = False           
        
    
            
    def doPCA(self,force = False):
        """performs PCA on the data
        creates dictionary self.PCADict with results
        Last column of PC matrix is first PC, second to last - second, etc. 
        
        Returns
        -------
        Dictionary of principal component matrices for different datasets
        """
        neededKeys = ["RemovedZeros","Corrected","FakedCis"]
        advicedKeys = ["TruncedTrans","RemovedPoor"]
        if (False in [i in self.appliedOperations for i in neededKeys]) and (force == False):
            print "needed operations:",neededKeys
            print "applied operations:", self.appliedOperations
            print "use 'force = True' to override this message" 
            raise StandardError("Critical filter not applied")
        if (False in [i in self.appliedOperations for i in advicedKeys]) and (force == False):
            print "Adviced operations:",advicedKeys
            print "Applied operations:", self.appliedOperations
            warnings.warn("Not all adviced filters applied")                        
        
        for i in self.dataDict.keys():
            self.PCDict[i] = PCA(self.dataDict[i])
        return self.PCDict
    
            
    def doEig(self,force = False):
        """performs eigenvector expansion on the data
        creates dictionary self.EigDict with results
        Last row of the eigenvector matrix is the largest eigenvector, etc. 
        
        Returns
        -------
        Dictionary of eigenvector matrices for different datasets
        """
        neededKeys = ["RemovedZeros","Corrected","FakedCis"]
        advicedKeys = ["TruncedTrans","RemovedPoor"]
        if (False in [i in self.appliedOperations for i in neededKeys]) and (force == False):            
            print "needed operations:",neededKeys
            print "applied operations:", self.appliedOperations
            print "use 'force = True' to override this message" 
            raise StandardError("Critical filter not applied")
        if (False in [i in self.appliedOperations for i in advicedKeys]) and (force == False):
            print "Not all adviced filters applied"
            print "Adviced operations:",advicedKeys
            print "Applied operations:", self.appliedOperations
            warnings.warn("Not all adviced filters applied")                        

        
        for i in self.dataDict.keys():
            self.EigDict[i] = EIG(self.dataDict[i])             
        return self.EigDict
    
    
    def cisToTrans(self,mode = "All", filename = "GM-all"):
        """
        Calculates cis-to-trans ratio.
        "All" - treating SS as trans reads
        "Dummy" - fake SS reads proportional to cis reads with the same total sum
        "Matrix" - use heatmap only        
        """
        data = self.dataDict[filename]
        cismap = self.chromosomeIndex[:,None] == self.chromosomeIndex[None,:]        
        cissums = numpy.sum(cismap * data, axis = 0) 
        allsums = numpy.sum(data,axis = 0)    
        if mode == "All":
            cissums += self.singlesDict[filename]
            allsums += self.singlesDict[filename]
        elif mode == "Dummy":
            sm = numpy.mean(self.singlesDict[filename])
            fakesm = cissums * sm / numpy.mean(cissums)
            cissums += fakesm
            allsums += fakesm
        elif mode == "Matrix":
            pass
        else:
            raise 
        return cissums / allsums
    


class binnedDataAnalysis(binnedData):
    """
    Class containing experimental features and data analysis scripts
    """
    
    def plotScaling(self,name,label = "BLA", color = None):
        "plots scaling of a heatmap,treating arms separately"
        data = self.dataDict[name]
        bins = numutils.logbins(2,self.genome.maxChrmArm/self.resolution,1.17)
        s = numpy.sum(data,axis = 0) > 0 
        mask = s[:,None] * s[None,:]
        chroms = []
        masks = []
        for i in xrange(self.chromosomeCount):
            beg = self.chromosomeStarts[i]
            end = self.centromerePositions[i]
            chroms.append(data[beg:end,beg:end])
            masks.append(mask[beg:end,beg:end])
            beg = self.centromerePositions[i]
            end = self.chromosomeEnds[i]
            chroms.append(data[beg:end,beg:end])
            masks.append(mask[beg:end,beg:end])         
        observed = []
        expected = []
        for i in xrange(len(bins) - 1):
            low = bins[i]
            high = bins[i+1]
            obs = 0
            exp = 0
            for j in xrange(len(chroms)):
                if low > len(chroms[j]): continue 
                high2 = min(high,len(chroms[j]))
                for k in xrange(low,high2):
                    obs += numpy.sum(numpy.diag(chroms[j],k))
                    exp += numpy.sum(numpy.diag(masks[j],k))                    
            observed.append(obs)
            expected.append(exp)
        observed = numpy.array(observed,float)
        expected = numpy.array(expected,float)
        values = observed/expected
        bins = numpy.array(bins,float)
        bins2 = 0.5 * (bins[:-1] + bins[1:])
        norm = numpy.sum(values * (bins[1:] - bins[:-1]))
        args = [self.resolution * bins2 / 1000000,values/(1. * norm)]
        if color != None: args.append(color)
        plt.plot(*args, label = label,linewidth = 2)
        


        
    def averageTransMap(self,name , mycmap = "hot_r",vmin = None,vmax = None):
        "plots and returns average inter-chromosomal inter-arm map"
        data = self.dataDict[name]
        #data = trunk(data,low = 0,high = 0.0001)         
        avarms = numpy.zeros((80,80))
        avmasks = numpy.zeros((80,80))
        discardCutoff = 20 
        
        for i in xrange(self.chromosomeCount):
            print i 
            for j in xrange(self.chromosomeCount):
                for k in [-1,1]:
                    for l in [-1,1]:
                        if i == j: continue 
                                                
                        cenbeg1 = self.chromosomeStarts[i] + self.genome.cntrStarts[i] / self.resolution
                        cenbeg2 = self.chromosomeStarts[j] + self.genome.cntrStarts[j] / self.resolution
                        cenend1 = self.chromosomeStarts[i] + self.genome.cntrEnds[i] / self.resolution
                        cenend2 = self.chromosomeStarts[j] + self.genome.cntrEnds[j] / self.resolution
                        
                        beg1 = self.chromosomeStarts[i]
                        beg2 = self.chromosomeStarts[j]
                        end1 = self.chromosomeEnds[i]
                        end2 = self.chromosomeEnds[j]
                        if k == 1:
                            bx = cenbeg1
                            ex = beg1-1
                            dx = -1                            
                        else:
                            bx = cenend1
                            ex = end1
                            dx = 1
                        if l ==1:
                            by = cenbeg2
                            ey = beg2-1
                            dy = -1
                        else:
                            by = cenend2
                            ey = end2
                            dy = 1
                            
                        
                        if abs(bx - ex) < discardCutoff: continue
                        if bx < 0: bx = None                          
                        if abs(by - ey) < discardCutoff: continue
                        if by < 0: by = None 
                                                 
                                                    
                        arms = data[bx:ex:dx,by:ey:dy]                        
                        assert max(arms.shape) <= self.genome.maxChrmArm / self.genome.resolution + 2
                                                
                        mx = numpy.sum(arms, axis = 0)
                        my = numpy.sum(arms, axis = 1)
                        maskx = mx == 0 
                        masky = my == 0
                        mask = (maskx[None,:] + masky[:,None]) == False                    
                        maskf = numpy.array(mask,float)
                        mlenx = (numpy.sum(mask, axis = 0) > 0 ).sum() 
                        mleny = (numpy.sum(mask, axis = 1) > 0 ).sum()
                        
                        if min(mlenx, mleny) < discardCutoff: continue
                    
                        add = numutils.zoomOut(arms,avarms.shape)
                        assert numpy.abs((arms.sum() - add.sum()) / arms.sum()) < 0.01
                                                                        
                        addmask = numutils.zoomOut(maskf,avarms.shape)
                        avarms += add 
                        avmasks += addmask
                        #mat_img(addmask)  
                          
        avarms /= numpy.mean(avarms)
        data = avarms / avmasks
        data /= numpy.mean(data) 
        plt.imshow(numpy.log(numutils.trunk(data)),cmap = "jet",interpolation = "nearest",vmin = vmin, vmax = vmax)
        removeBorder()
        return numpy.log(numutils.trunk(data))
            
    def perArmCorrelation(self,data1,data2,doByArms = []):
        """does inter-chromosomal spearman correlation 
        of two vectors for each chromosomes separately.
        Averages over chromosomes with weight of chromosomal length
        For chromosomes in "doByArms" treats arms as separatre chromosomes
        returns average Spearman r correlation
        """
        cr = 0 
        ln = 0
        for i in xrange(self.chromosomeCount):
            if i in doByArms:
                beg = self.chromosomeStarts[i]
                end = self.centromerePositions[i]
                if end > beg: 
                    cr += (abs(spearmanr(data1[beg:end],data2[beg:end])[0])) * (end - beg)            
                    ln += (end-beg)
                    print spearmanr(data1[beg:end],data2[beg:end])[0]
                beg = self.centromerePositions[i]
                end = self.chromosomeEnds[i]
                if end > beg: 
                    cr += (abs(spearmanr(data1[beg:end],data2[beg:end])[0])) * (end - beg)            
                    ln += (end-beg)
                    print spearmanr(data1[beg:end],data2[beg:end])[0]
            else: 
                beg = self.chromosomeStarts[i]
                end = self.chromosomeEnds[i]
                if end > beg: 
                    cr += (abs(spearmanr(data1[beg:end],data2[beg:end])[0])) * (end - beg)            
                    ln += (end-beg)                            
        return cr/ln 
 
 
    def divideOutAveragesPerChromosome(self):
        "divides each interchromosomal map by it's mean value"
        mask2D = self._giveMask2D()
        for chrom1 in xrange(self.chromosomeCount):
            for chrom2 in xrange(self.chromosomeCount):
                for i in self.dataDict.keys():
                    value = self.dataDict[i]
                    submatrix = value[self.chromosomeStarts[chrom1]:self.chromosomeEnds[chrom1],
                                      self.chromosomeStarts[chrom2]:self.chromosomeEnds[chrom2]]
                    masksum = numpy.sum(mask2D[self.chromosomeStarts[chrom1]:self.chromosomeEnds[chrom1],
                                      self.chromosomeStarts[chrom2]:self.chromosomeEnds[chrom2]])
                    valuesum = numpy.sum(submatrix)
                    mean = valuesum / masksum
                    submatrix /= mean
 

    def interchromosomalValues(self, filename = "GM-all",returnAll = False):
        "returns average inter-chromosome-interaction values, ordered always the same way"
        values = self.chromosomeIndex[:,None] + self.chromosomeCount * self.chromosomeIndex[None,:]
        values[self.chromosomeIndex[:,None] == self.chromosomeIndex[None,:]] = self.chromosomeCount * self.chromosomeCount - 1 
        #mat_img(values)
        uv = numpy.sort(numpy.unique(values))[1:-1] 
        probs = numpy.bincount(values.ravel(),weights = self.dataDict[filename].ravel())
        counts = numpy.bincount(values.ravel())
        if returnAll == False: return probs[uv] / counts[uv]
        else: 
            probs[self.chromosomeCount * self.chromosomeCount - 1] = 0 
            values = probs/counts
            values[counts==0] = 0             
            #mat_img(values.reshape((22,22)))
            return values.reshape((self.chromosomeCount, self.chromosomeCount))


    
class experimentalBinnedData(binnedData):
    "Contains some poorly-implemented new features"        
    def emulateCis(self):
        """if you want to have fun creating syntetic data, this emulates cis contacts. 
        adjust cis/trans ratio in the C code"""
        transmap = self.chromosomeIndex[:,None] == self.chromosomeIndex[None,:]
        len(transmap)
        for i in self.dataDict.keys():
            data = self.dataDict[i] * 1. 
            #mask = na(self.chromosomeIndex[:,None] == self.chromosomeIndex[None,:],int)
            N = len(data)
            N
            code = r"""
            #line 841 "binary_search.py"
            using namespace std; 
            for (int i = 0; i < N; i++)    
            {                 
                for (int j = 0; j<N; j++)
                {
                    if (transmap[N * i + j] == 1)
                    {
                        data[N * i + j] = data[N * i +j] * 300 /(abs(i-j) + 0.5);                     
                    }
                }
            } 
            """
            support = """
            #include <math.h>
            """
            weave.inline(code, ['transmap','data',"N"], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
            self.dataDict[i] = data
        self.removedCis = False 
        self.fakedCis = False
         


    def fakeMissing(self,iterative = True):
        """fakes megabases that have no reads. For cis reads fakes with cis reads at the same distance. For trans fakes with random trans read at the same diagonal. 
        """
        for i in self.dataDict.keys():
            data = self.dataDict[i] * 1.
            sm = numpy.sum(data,axis = 0) > 0  
            mask = sm[:,None] * sm[None,:]
            transmask = numpy.array(self.chromosomeIndex[:,None] == self.chromosomeIndex[None,:],int)
            #mat_img(transmask)
            N = len(data)
            N,transmask ,mask #to remove warning
            code = r"""
            #line 841 "binary_search.py"
            using namespace std; 
            for (int i = 0; i < N; i++)    
            {    
                for (int j = i; j<N; j++)
                {
                    if ((MASK2(i,j) == 0) )                    
                    {
                    for (int ss = 0; ss < 401; ss++) 
                        {
                        int k = 0;                            
                        int s = rand() % (N - (j-i));                        
                        if ((mask[s * N + s + j - i] == 1) && ((transmask[s * N + s + j - i] == transmask[i * N + j]) || (ss > 200)) )
                                {                                
                                data[i * N + j] = data[s * N + s + j - i];
                                data[j * N + i] = data[s * N + s + j - i];
                                break;
                                }
                        if (ss == 400) {printf("Cannot fake one point... skipping %d %d \n",i,j);}
                            
                            
                        }
                    }
                }
            } 
            """
            support = """
            #include <math.h>
            """
            for s in xrange(5):
                s #to remove warning
                weave.inline(code, ['transmask','mask','data',"N"], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
                #mat_img(numpy.log(data+1),trunk = True)
                data = correct(data)              
            self.dataDict[i] = data
             
            
            #mat_img(self.dataDict[i]>0)
            

    

    def iterativeCorrectByTrans(self,names = None):        
        """performs iterative correction by trans data only, corrects cis also
        
        Parameters
        ----------
        names : list of str or None, optional
            Keys of datasets to be corrected. By default, all are corrected.
        """

        if names == None: names = self.dataDict.keys()        
        self.transmap = self.chromosomeIndex[:,None] != self.chromosomeIndex[None,:]
        #mat_img(self.transmap)
        for i in names:
            data = self.dataDict[i]
            self.dataDict[i],self.biasDict[i] = numutils.ultracorrectSymmetricByMask(data,self.transmap,50)
            try: self.singlesDict[i] /= self.biasDict[i]
            except: print "bla"
    
    def loadWigFile(self,filename,label,wigFileType = "Auto"):
        filename = os.path.abspath(filename) 
        if wigFileType == "Auto":
            ext = os.path.splitext(filename)[1]
            if ext == "":
                raise StandardError("Wig file has no extension. Please specify it's type")
            elif ext.lower()  == ".wig":
                wigFileType = "wig"
            elif ext.lower() == ".bigwig":
                wigFileType = "bigwig"
            else: raise StandardError("Unknown extension of wig file: %s" % ext)
                         
        if wigFileType.lower() == "wig": 
            data = self.genome.parseFixedStepWigAtKbResolution(filename)
        elif wigFileType.lower() == "bigwig": 
            data = self.genome.parseBigWigFile(filename,resolution = 1000,divideByValidCounts = False)
        else:
            raise StandardError("Wrong type of wig file : %s" % wigFileType) 
        
        if self.genome.resolution % 1000 != 0: raise StandardError("Cannot parse wig file at non-kb resolution")
                
        vector = numpy.zeros(self.genome.numBins,float)
        for chrom,value in enumerate(data):
            value = numpy.array(value)             
            value.resize(self.genome.chrmLensBin[chrom] * (self.resolution/1000))            
            value.shape = (-1,self.genome.resolution / 1000 )
            if value.mean() == 0:
                raise StandardError("Chromosome %s contains zero data in wig file %s" % (self.genome.idx2label[chrom],filename))
            mask = value == 0
            value = numpy.log(value) 
            value[mask] = 0 
            av = numpy.sum(value,axis = 1) / numpy.sum(mask == False,axis = 1)
            av[numpy.isfinite(av) == False] = numpy.NAN
            
            vector[self.genome.chrmStartsBinCont[chrom]:self.genome.chrmEndsBinCont[chrom]] = av
        vector[numpy.isnan(vector)] = numpy.median(vector) 
        self.trackDict[label] = vector
            
             
            

            
            
            
            
         
