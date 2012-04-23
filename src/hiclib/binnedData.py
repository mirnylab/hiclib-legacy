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

All information about the genome, including GC content and restriction sites, can be obtained from the Genome class. 

Genomic tracks can be loaded using an automated parser that accepts bigWig files and fixed step wiggle files. 
See documentation for :py:func:`experimentalBinnedData.loadWigFile` that describes exactly how the data is averaged and parsed.    


Variables
---------

self.dataDict - dictionary with heatmaps; keys are provided when loading the data.

self.singlesDict - dictionary with SS read vectors. Keys are the same.

self.fragsDict - dictionary with fragment density data  

self.trackDict - dictionary with genomic tracks, such as GC content. Custom tracks should be added here. 

self.biasDict - dictionary with biases as calculated by iterative correction (incomplete)  

self.PCDict - dictionary with principal components of each datasets. Keys as in dataDict

self.EigEict - dictionary with eigenvectors for each dataset. Keys as in datadict. 


Hierarchy of filters
--------------------

hierarchy

--------------------------------------------------------------- 
"""

import os 
from mirnylab import numutils
import warnings
from mirnylab.plotting import  removeBorder , mat_img
from mirnylab.numutils import PCA, EIG,correct, ultracorrectSymmetricWithVector,\
    isInteger
from mirnylab.genome import Genome 
import  numpy
from math import exp
from mirnylab.h5dict import h5dict  
from scipy import weave 
from scipy.stats.stats import spearmanr
import matplotlib.pyplot as plt 
from mirnylab.systemutils import setExceptionHook


    
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
        self.eigEigenvalueDict = {}
        self.PCAEigenvalueDict = {}
        self.dicts = [self.trackDict, self.biasDict, self.singlesDict, self.fragsDict]
        self.eigDicts = [self.PCDict, self.EigDict]
        self._loadGC() 
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

    def _loadGC(self):        
        "loads GC content at given resolution"
        self.trackDict["GC"] = numpy.concatenate(self.genome.GCBin)
        
    def _checkItertiveCorrectionError(self):
        "internal method for checking if iterative correction might be bad to apply"
        for value in self.dataDict.values():           
             
            if isInteger(value).all() == True:
                s = numpy.sum(value,axis = 0)
                sums = numpy.sort(s[s!=0])
                if sums[0] < 100:
                    error = int(100. / numpy.sqrt(sums[0]))
                    message1 =  "Lowest 5 sums of an array rows are: " + str(sums[:5])                    
                    warnings.warn("\n%s\nIterative correction will lead to about %d %% relative error for certain columns" % (message1,error))
                    
                    if sums[0] < 5: 
                        raise StandardError("Iterative correction is very dangerous. Use force=true to override.") 
                          
            else:
                s = numpy.sum(value>0,axis = 0)
                sums = numpy.sort(s[s!=0])
                if sums[0] < min(100,len(value)/2):
                    error = int(100. / numpy.sqrt(sums[0]))
                    print "Got floating-point array for correction. Rows with 5 least entrees are:",sums[:5]
                    warnings.warn("\nIterative correction might lead to about %d %% relative error for certain columns" % error)
                    if sums[0] < 4: 
                        raise StandardError("Iterative correction is very dangerous. Use force=true to override.")
                    
    def _checkAppliedOperations(self,neededKeys=[],advicedKeys=[],excludedKeys=[]):
        "Internal method to check if all needed operations were applied"
        
        if (True in [i in self.appliedOperations for i in excludedKeys]):
            print "Operations that are not allowed:", excludedKeys
            print "applied operations: ", self.appliedOperations
            print "use 'force = True' to override this message" 
            raise StandardError("Prohibited filter was applied")
        
        if (False in [i in self.appliedOperations for i in neededKeys]): 
            print "needed operations:",neededKeys
            print "applied operations:", self.appliedOperations
            print "use 'force = True' to override this message" 
            raise StandardError("Critical filter not applied")
        
        if (False in [i in self.appliedOperations for i in advicedKeys]):
            print "Adviced operations:",advicedKeys
            print "Applied operations:", self.appliedOperations
            warnings.warn("\nNot all adviced filters applied")                        
        
    
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
            countsum = numpy.sum(data,axis = 0)  
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
            lim = numpy.percentile(data[transmask],100.*(1 - high))
            print "dataset %s truncated at %lf" % (i,lim) 
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
        """Used to fake cis counts or any other region with random trans counts. 
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
        for key in self.dataDict.keys():                         
            data = numpy.asarray(self.dataDict[key],order = "C",dtype = float)
            if mask == "CisCounts": 
                _mask =  numpy.array(self.chromosomeIndex[:,None] == self.chromosomeIndex[None,:],int,order = "C")
            else:
                assert mask.shape == self.dataDict.values()[0].shape  #check that mask has correct shape 
                _mask = numpy.array(mask,dtype = int, order = "C")
                _mask[self.chromosomeIndex[:,None] == self.chromosomeIndex[None,:]] = 2   #do not fake with cis counts                                             
            s = numpy.abs(numpy.sum(data,axis = 0)) <= 1e-10 
            _mask[:,s]= 2
            _mask[s,:] = 2              
            N = len(data)
            N
            code = r"""
            #line 310 "binnedData.py"
            using namespace std; 
            for (int i = 0; i < N; i++)    
            {    
                for (int j = i; j<N; j++)
                {
                    if (_mask[i* N + j] == 1)                    
                    {
                    while (true) 
                        {
                        int r = rand() % 2;                         
                        if (r == 0)
                            {                            
                            int s = rand() % N;
                            if (_mask[i * N + s] == 0)
                                {                                
                                data[i * N + j] = data[i * N + s];
                                data[j * N + i] = data[i * N + s];
                                break;
                                 
                                }
                            }
                        else
                            {
                            int s = rand() % N;
                            if (_mask[j * N + s] == 0)
                                {
                                data[i * N + j] = data[j * N + s];
                                data[j * N + i] = data[j * N + s]; 
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
            weave.inline(code, ['_mask','data',"N"], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
            self.dataDict[key] = data
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
        
    def fakeTranslocations(self,translocationRegions):
        """
        This method fakes reads corresponding to a translocation.
        
        Parameters
        ----------
        
        translocationRegions: list of tuples
            List of tuples (chr1,start1,end1,chr2,start2,end2), masking a high-count region around visible translocation.
            If chromosome end is None, it is treated as length of chromosome. So, use (chr1,0,None,chr1,0,None) to remove map entirely.          
        """
        self._checkAppliedOperations(excludedKeys = "RemovedZeros")
        mask = numpy.zeros((self.genome.numBins, self.genome.numBins),int)
        resolution = self.genome.resolution           
        
        for i in translocationRegions:
            st1 = self.genome.chrmStartsBinCont[i[0]]
            st2 = self.genome.chrmStartsBinCont[i[3]]
            beg1 = st1 +i[1]/resolution
            if i[2] != None: 
                end1 = st1 + i[2]/resolution+1
            else:
                end1 = self.genome.chrmEndsBinCont[i[0]]
            beg2 = st2 + i[4]/resolution
            if i[5] != None: 
                end2 = st2 + i[5]/resolution+1
            else:
                end2 = self.genome.chrmEndsBinCont[i[3]]

            mask[beg1:end1,beg2:end2] = 1
            mask[beg2:end2,beg1:end1] = 1
        self.fakeCisOnce(mask)
         
            
            


    def correct(self,names=None):
        """performs single correction without SS
        
        Parameters
        ----------
        names : list of str or None
            Keys of datasets to be corrected. If none, all are corrected. 
        """
        self.iterativeCorrectWithoutSS(names, M=1)                                          
                
    def iterativeCorrectWithoutSS(self, names = None,M=50,force = False):
        """performs iterative correction without SS
        
        Parameters
        ----------
        names : list of str or None, optional
            Keys of datasets to be corrected. By default, all are corrected.
        M : int, optional
            Number of iterations to perform. 
        """
        
        
        
        
        if force == False:
            self._checkItertiveCorrectionError()             
            self._checkAppliedOperations(advicedKeys = ["RemovedDiagonal","RemovedPoor"])        
                    
        if names == None: names = self.dataDict.keys()
        for i in names:
            self.dataDict[i] = ultracorrectSymmetricWithVector(self.dataDict[i],M=M)[0]         
          
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
        
        if force == False: 
            self._checkAppliedOperations(advicedKeys = ["RemovedDiagonal","RemovedPoor"], 
                                        excludedKeys = ["Corrected","RemovedCis"])
            self._checkItertiveCorrectionError()
                
        if names == None: names = self.dataDict.keys()
        for i in names:
            data = self.dataDict[i]
            vec = self.singlesDict[i]
            ndata,nvec = ultracorrectSymmetricWithVector(data, vec,M=M)                         
            self.dataDict[i] = ndata
            self.singlesDict[i] = nvec
            vec[nvec==0] = 1
            nvec[nvec==0] = 1
            self.biasDict[i] = (vec / nvec)

        
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
            warnings.warn("\nYou're removing zeros twice. You can't restore zeros now!")
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
        if force == False: self._checkAppliedOperations(neededKeys, advicedKeys)
        
        for i in self.dataDict.keys():
            currentPCA,eigenvalues = PCA(self.dataDict[i])
            self.PCAEigenvalueDict[i] = eigenvalues
            for j in xrange(len(currentPCA)):
                if spearmanr(currentPCA[j],self.trackDict["GC"])[0] < 0:
                    currentPCA[j] = -currentPCA[j]
            self.PCDict[i] = currentPCA
        return self.PCDict 
    
            
    def doEig(self,numPCs = 3,force = False):
        """performs eigenvector expansion on the data
        creates dictionary self.EigDict with results
        Last row of the eigenvector matrix is the largest eigenvector, etc. 
        
        Returns
        -------
        Dictionary of eigenvector matrices for different datasets
        """
        neededKeys = ["RemovedZeros","Corrected","FakedCis"]
        advicedKeys = ["TruncedTrans","RemovedPoor"]
        if force == False: self._checkAppliedOperations(neededKeys, advicedKeys)

        
        for i in self.dataDict.keys():             
            currentEIG,eigenvalues = EIG(self.dataDict[i],numPCs = numPCs)
            self.eigEigenvalueDict[i] = eigenvalues 
            for j in xrange(len(currentEIG)):
                if spearmanr(currentEIG[j],self.trackDict["GC"])[0] < 0:
                    currentEIG[j] = -currentEIG[j]
            self.EigDict[i] = currentEIG                                      
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
        avarms = numpy.zeros((80,80))
        avmasks = numpy.zeros((80,80))
        discardCutoff = 10
                
        
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
                        if ex <0: ex = None                          
                        if abs(by - ey) < discardCutoff: continue
                        if by < 0: by = None
                        if ey < 0: ey = None 
                                                 
                                                    
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
                        assert numpy.abs((arms.sum() - add.sum()) / arms.sum()) < 0.02
                        
                         
                                                                        
                        addmask = numutils.zoomOut(maskf,avarms.shape)
                        avarms += add 
                        avmasks += addmask
                          
        
        avarms /= numpy.mean(avarms)
        data = avarms / avmasks
        data /= numpy.mean(data) 
        plt.imshow(numpy.log(numutils.trunc(data)),cmap = "jet",interpolation = "nearest",vmin = vmin, vmax = vmax)
        removeBorder()
        return numpy.log(numutils.trunc(data))
            
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
    
    def loadWigFile(self,filename,label,control = None,wigFileType = "Auto"):
        """Currently only fixed-step wig files and bigWig files are supported!!!
        Import from fixedStep wig files is very fast, however is not the most reliable.
        
        TODO: rewrite parser to support any resolution.  
        
        for VariableStep files use wigToBigWig utility to convert them to bigWig format first. 
        To use it you will also need to have fetchChromSizes script. 
        
        Then you just run 
        $bash fetchChromSizes hg18 > hg18.chrom.sizes
        $./wigToBigWig myWig.wig hg18.chrom.sizes myWig.bigWig
        
        And you enjoy your favourite bigWig. 
        
        BigWig import is implemented using bx-python module.
        It is normally very fast; however, it has a bug at low resolutions. 
        I have an ugly workaround for it (chopping the quiery into many smaller pieces), but I hope 
        that they actually fix this bug. 
        
        Anyway, check their repo on BitBucket, maybe they've fixed my issue # 39 :) 
        https://bitbucket.org/james_taylor/bx-python/overview
        """
        
        filename = os.path.abspath(filename)
        
        def loadFile(name,wigFileType = wigFileType):
            """Choosing right method to load wigfile"""
            
            if os.path.exists(name) == False:
                raise IOError("\n Wig file not found : %s " % (os.path.abspath(name)))
                        
            if wigFileType == "Auto":
                ext = os.path.splitext(name)[1]
                if ext == "":
                    raise StandardError("Wig file has no extension. Please specify it's type")
                elif ext.lower()  == ".wig":
                    if open(filename).readline()[:2] != "fi":
                        raise StandardError("Cannot read non fixed-step wig files! Please use wigToBigWig utility. See docstring of this method.")
                    wigFileType = "wig"
                elif ext.lower() == ".bigwig":
                    wigFileType = "bigwig"
                else: raise StandardError("Unknown extension of wig file: %s" % ext)
                             
            if wigFileType.lower() == "wig": 
                data = self.genome.parseFixedStepWigAtKbResolution(name,resolution = 5000)
            elif wigFileType.lower() == "bigwig": 
                data = self.genome.parseBigWigFile(name,resolution = 5000,divideByValidCounts = False)
            else:
                raise StandardError("Wrong type of wig file : %s" % wigFileType)
            return data
        
        "Loading the data files"
        data = loadFile(filename)
        if control != None: 
            controlData = loadFile(control) 
        
        if self.genome.resolution % 5000 != 0: raise StandardError("Cannot parse wig file at resolution that is not a multiply of 5 kb")
                 
        vector = numpy.zeros(self.genome.numBins,float)  #resulting array
        for chrom,value in enumerate(data):            
            value = numpy.array(value)            
            if control != None: 
                chromControl = numpy.asarray(controlData[chrom])
                vmask = value != 0 
                cmask = chromControl != 0
                keepmask = vmask * cmask
                vmasksum,cmasksum = vmask.sum(),cmask.sum()
                if max(vmasksum,cmasksum) / (1. * min(vmasksum,cmasksum)) > 1.3:
                    warnings.warn("\nBig deviation: number of non-zero data points: %s, control points:%s." % (vmasksum,cmasksum))    
                value[-keepmask] = 0 
                value[keepmask] = value[keepmask] / chromControl[keepmask]                    
                                          
            value.resize(self.genome.chrmLensBin[chrom] * (self.resolution/5000))            
            value.shape = (-1,self.genome.resolution / 5000 )
            if value.mean() == 0:
                raise StandardError("Chromosome %s contains zero data in wig file %s" % (self.genome.idx2label[chrom],filename))
            mask = value == 0
            value[-mask] = numpy.log(value[-mask]) 
              
            valuesum = numpy.sum(value,axis = 1)
            masksum = numpy.sum(mask == False,axis = 1)
            valuesum[masksum==0] = 0 
            vmask = valuesum != 0 
            valuesum[vmask] /= masksum[vmask]
            valuesum[-vmask] = numpy.median(valuesum[vmask])  #setting all unknown points to the median of known points             
                        
            vector[self.genome.chrmStartsBinCont[chrom]:self.genome.chrmEndsBinCont[chrom]] = valuesum        
        if len(vector) != self.genome.numBins: 
            raise ValueError("Length mismatch. Length of vector: %d, length of genome:%d" % 
                             (len(vector), self.genome.numBins)) 
        self.trackDict[label] = vector
            
             
    def loadErezEigenvector1MB(self,erezFolder):
        "Loads Erez chromatin domain eigenvector for HindIII"
        if self.resolution != 1000000: raise StandardError("Erez eigenvector is only at 1MB resolution")
        if self.genome.folderName != "hg18": raise StandardError("Erez eigenvector is for hg18 only!")
        folder = os.path.join(erezFolder,"GM-combined.ctgDATA1.ctgDATA1.1000000bp.hm.eigenvector.tab")
        folder2 = os.path.join(erezFolder,"GM-combined.ctgDATA1.ctgDATA1.1000000bp.hm.eigenvector2.tab")
        eigenvector = numpy.zeros(self.genome.numBins,float)
        for chrom in range(1,24):
            filename = folder.replace("DATA1",str(chrom))
            if chrom in [4,5]: 
                filename = folder2.replace("DATA1",str(chrom))
            mydata = numpy.array([[float(j) for j in i.split()] for i in open(filename).readlines()])
            eigenvector[self.genome.chrmStartsBinCont[chrom-1]+ numpy.array(mydata[:,1],int)] = mydata[:,2]
        self.trackDict["Erez"] = eigenvector            

            
    def loadTanayDomains(self):             
        "domains, extracted from Tanay paper image"
        if self.genome.folderName != "hg18":
            raise StandardError("Tanay domains work only with hg18")
        data = """0 - 17, 1 - 13.5, 2 - 6.5, 0 - 2, 2 - 2; x - 6.5, 0 - 6, 1 - 13.5, 0 - 1.5, 1 - 14.5 
    1 - 8.5, 0 - 2.5, 1 - 14, 2 - 6; 0 - 1.5, 2 - 11.5, 1 - 35
    1 - 14, 0-6, 2 - 11; 2 - 4.5, 1 - 5, 0 - 4, 1 -20.5, 0 - 2
    0 - 3, 2 - 14; 2 - 5, 1 - 42
    2 - 16; 2 - 7, 0 - 3, 1 - 18.5, 0 - 1, 1 - 13, 0 - 2.5 
    0 - 2, 1 - 6.5, 0 - 7.5, 2 - 4; 2 - 6, 1 - 31
    0 - 2, 1 - 11, 2 - 7; 2 - 7.5, 1 - 5, 0 - 3, 1 - 19 
    2 - 9.5, 0 - 1, 2 - 5; 2 - 4, 1 - 27.5, 0 - 2.5
    2 - 11.5, 0 - 2.5, x - 2.5; x - 5, 2 - 8, 0 - 3.5, 1 - 9, 0 - 6
    2 - 13.5; 2 - 9, 0 - 3, 1 - 6, 0 - 3.5, 1 - 10.5
    0 - 3.5, 2 - 15; 2 - 1, 0 - 7.5, 1 - 13, 0 - 1.5, 1 - 4
    0 - 4, 2 - 8; 2 - 2, 0 - 5, 2 - 2.5, 1 - 13, 0 - 6.5, 1 - 3.5 
    x - 5.5; 2 - 8.5, 0 - 1, 2 - 7, 1 - 16
    x - 5.5; 2 - 14.5, 0 - 6, 2 - 3, 1 - 2.5, 2 - 1, 0 - 3
    x - 5.5; 2 - 6, 0 - 3.5, 2 - 1.5, 0 - 11.5, 2 - 5.5
    0 - 11, 2 - 1; x - 2.5, 2 - 6.5, 0 - 3, 2 - 2, 0 - 3.5 
    0 - 4, 2 - 1.5, 0 - 1.5; 0 - 19
    2 - 5; 2 - 20
    0 - 9.5, x - 1.5; x - 1, 2 - 2, 0 - 8.5
    0 - 2, 2 - 7; 0 - 8, 2 - 2, 0 - 1
    x - 0.5; 2 - 8.5, 0 - 3
    x - 4; 0 -12 
    x - 1.5, 1 - 13, 2 - 5.5; 2 - 2, 1 - 29"""
        chroms = [i.split(";") for i in data.split("\n")]
        result = []     
        for chrom in chroms:
            result.append([])
            cur = result[-1] 
            for arm in chrom:
                
                for enrty in arm.split(","):
                    spentry = enrty.split("-")
                    if "x" in spentry[0]: value = -1 
                    else: value = int(spentry[0])
                    cur += ([value] * int(2 * float(spentry[1])))
                cur += [-1] * 2
        #lenses = [len(i) for i in result]
            
        domains = numpy.zeros(self.genome.numBins,int)
        for i in xrange(self.genome.chrmCount):
            for j in xrange((self.genome.chrmLens[i] / self.resolution)):
                domains[self.genome.chrmStartsBinCont[i] + j] = result[i][(j * len(result[i])/ ((self.genome.chrmLens[i] /self.resolution))) ]
        self.trackDict['TanayDomains'] = domains 
            
            
            
         
