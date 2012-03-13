import base,dnaUtils,numutils
base  # Eclipse warning damper
from plotting import mat_img, removeBorder 
from numutils import PCA, EIG,correct, ultracorrectSymmetricWithVector

import  cPickle, numpy,joblib,scipy.stats  
from math import exp 
from scipy import weave

from scipy.stats.stats import spearmanr
import matplotlib.pyplot as plt 



    
class binnedData(object):
    "base class to work with binned data"
    
    def __init__(self, resolution,genomeFolder):
        self.genome = dnaUtils.Genome(genomeFolder = genomeFolder)
        if resolution != None: self.resolution = resolution
        self.chromosomes = self.genome.chromosomes
        self.resolution = resolution
        self.chromosomeExtensionLength = 0   #length of the gap between chromosomes                
        self.genome.createMapping(self.resolution)                        
        self.initChromosomes()        
        self.dataDict = {}
        self.biasDict = {}
        self.trackDict = {}
        self.singlesDict = {}
        self.fragsDict = {}
        
    def initChromosomes(self):
        "loads mappings from the genome class based on resolution"
        self.genome.createMapping(self.resolution,chromosomeExtensionLength = self.chromosomeExtensionLength)
        self.chromosomeStarts = self.genome.chromosomeStarts
        self.centromerePositions = self.genome.centromerePositions
        self.chromosomeEnds = self.genome.realChromosomeEnds
        self.trackLength = self.genome.chromosomeEnds[-1]
        self.N = self.genome.N
        self.chromosomeCount = self.genome.chromosomeCount
        self.chromosomeIndex = self.genome.chromosomeIndex
        self.positionIndex = self.genome.positionIndex        
        self.armIndex = self.chromosomeIndex * 2 + numpy.array(self.positionIndex > self.genome.centromeres[self.chromosomeIndex],int)
                
    
    def simpleLoad(self,filename,name,checkGenomes = True):
        "load standard type data"
        alldata = dict(joblib.load(filename))        
        self.dataDict[name] = alldata["heatmap"]        
        self.singlesDict[name] = alldata["singles"]
        self.fragsDict[name] = alldata["frags"]
        if self.genome.type != alldata["genome"]: 
            if checkGenomes == True: 
                print "Genome name mismatch!!!"
                print "source genome",alldata["genome"]
                print "our genome",self.genome.type
                self.exit()
        try: self.resolution
        except: self.resolution = alldata["resolution"]
        if self.resolution != alldata["resolution"]:
            print "resolution mismatch!!!"
            print "--------------> Bye <-------------"
            raw_input("Press any key to continue... ") 
            exit()  
        
        
    def loadGC(self):        
        "loads GC content at given resolution"
        data = self.genome.getBinnedGCContent(self.resolution)
        eigenvector = numpy.zeros(self.trackLength,float)        
        for chrom in range(1,self.chromosomeCount + 1):
            eigenvector[self.chromosomeStarts[chrom-1]:self.chromosomeStarts[chrom-1] + len(data[chrom-1])] = data[chrom-1]
        self.trackDict["GC"] = eigenvector
  
    
    def removeDiagonal(self,m):
        "does what it says from all datasets"
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
            #mat_img(data) 
            
        
        
    
    def removeStandalone(self,offset = 3):
        "removes standalone bins (groups of less-than-offset-long bins)"
        self.N = len(self.dataDict.values()[0])
        self.mask = numpy.ones(len(self.dataDict.values()[0]),numpy.bool)
        for data in self.dataDict.values():
            datasum = numpy.sum(data,axis = 0)
            datamask = datasum > 0
            self.mask *= datamask  
        diffs = numpy.diff(numpy.array(numpy.r_[False, self.mask, False],int))
        begins = numpy.nonzero(diffs == 1)[0] 
        ends = numpy.nonzero(diffs == -1)[0]
        beginsmask = (ends - begins) <= offset
        newbegins = begins[beginsmask]
        newends = ends[beginsmask]        
        print "removing %d standalone megabases"% numpy.sum(newends - newbegins)
        for i in xrange(len(newbegins)): self.mask[newbegins[i]:newends[i]] = False
        self.mask2D = self.mask[:,None] * self.mask[None,:]
        for i in self.dataDict.values(): i[self.mask2D == False] = 0
        
    def removePoorRegions(self,names = None, cutoff = 2):
        "removes cutoff persent of bins with least counts"
        statmask = numpy.zeros(len(self.dataDict.values()[0]),numpy.bool)
        self.mask = numpy.ones(len(self.dataDict.values()[0]),numpy.bool)
        if names == None: names =self.dataDict.keys()  
        for i in names:
            data = self.dataDict[i]
            datasum = numpy.sum(data,axis = 0)            
            datamask = datasum > 0
            self.mask *= datamask  
            try: countsum = numpy.sum(data,axis = 0) + self.singlesDict[i]
            except: countsum = numpy.sum(data,axis = 0) 
            newmask = countsum >= numpy.percentile(countsum[datamask],cutoff)
            self.mask *= newmask  
            statmask [(newmask == False) * (datamask == True)] = True
        print "removed %d poor megabases", statmask.sum() 
        self.mask2D = self.mask[:,None] * self.mask[None,:]
        for i in self.dataDict.values(): i[self.mask2D == False] = 0
              
    def TrunkTrans(self,high = 0.0005):
        "trunkates trans contacts to remove blowouts"
        for i in self.dataDict.keys():
            data = self.dataDict[i]
            transmask = self.chromosomeIndex[:,None] != self.chromosomeIndex[None,:]
            lim = numpy.percentile(data[transmask],100*(1 - high))
            tdata = data[transmask]
            tdata[tdata > lim] = lim            
            self.dataDict[i][transmask] = tdata                        
        
    def removeCis(self):
        "sets to zero all cis contacts"
        mask = self.chromosomeIndex[:,None] == self.chromosomeIndex[None,:]
        self.mask2D[mask] = False 
        for i in self.dataDict.keys():                
            self.dataDict[i][mask] = 0   
        self.removedCis = True           
        
            
    def fakeCisOnce(self):
        for i in self.dataDict.keys():
            data = self.dataDict[i] * 1. 
            data = numpy.array(data,order = "C")
            mask = numpy.array(self.chromosomeIndex[:,None] == self.chromosomeIndex[None,:],int)
            s = numpy.abs(numpy.sum(data,axis = 0)) <= 1e-20 
            mask[:,s]= 2
            mask[s,:] = 2              
            N = len(data)
            N
            code = r"""
            #line 841 "binary_search.py"
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
            self.removedCis = True 
            self.fakedCis = True 
            #mat_img(data)

    def fakeCis(self):
        "fakes cis contacts in a fancy way"
        self.removeCis()
        self.ultracorrect(M=5)
        self.fakeCisOnce()
        self.ultracorrect( M = 5)
        self.fakeCisOnce()
        self.ultracorrect(M = 10) 

    def emulateCis(self):
        """if you want to have fun creating syntetic data, this emulates cis contacts. 
        adjust cis/trans ratio in the code"""
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


    def fakeMissing(self, stay = False):
        "fakes missing megabases in the fancy way"
        for i in self.dataDict.keys():
            data = self.dataDict[i] * 1.
            sm = numpy.sum(data,axis = 0) > 0  
            mask = sm[:,None] * sm[None,:]
            transmask = numpy.array(self.chromosomeIndex[:,None] == self.chromosomeIndex[None,:],int)
            #mat_img(transmask)
            N = len(data)
            N,transmask  #to remove warning
            code = r"""
            #line 841 "binary_search.py"
            using namespace std; 
            for (int i = 0; i < N; i++)    
            {    
                for (int j = i; j<N; j++)
                {
                    if ((mask[i* N + j] == 0) )                    
                    {
                    while (true) 
                        {
                        int k = 0;                            
                        int s = rand() % (N - (j-i-1));
                        if (j - i > 100)   k = rand() % 2;
                        if ((mask[s * N + s + k + j - i] == 1) && (transmask[s * N + s + j - i] == transmask[i * N + j])) 
                                {                                
                                data[i * N + j] = data[s * N + s + j - i];
                                data[j * N + i] = data[s * N + s + j - i];
                                break;
                                }
                        }
                    }
                }
            } 
            """
            support = """
            #include <math.h>
            """
            for s in xrange(10):
                s #to remove warning
                weave.inline(code, ['transmask','mask','data',"N"], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
                #mat_img(numpy.log(data+1),trunk = True)
                data = correct(data)
            if stay == False: data[mask == 0] = 0  
            self.dataDict[i] = data
            #mat_img(self.dataDict[i]>0) 


    def correct(self,names=None):
        "performs single correction without SS"
        self.ultracorrect(names,M=1)
        
    def ultracorrect(self, names = None,M=50):
        "performs iterative correction wihtout SS"
        if names == None: names = self.dataDict.keys()
        
        for i in names:
            self.dataDict[i] = ultracorrectSymmetricWithVector(self.dataDict[i],M=M)[0]

    def ultracorrectAll(self,names = None,M = 55):
        "performs iterative correction with SS"
        if names == None: names = self.dataDict.keys()
        for i in names:
            data = self.dataDict[i]
            vec = self.singlesDict[i]
            ndata,nvec = ultracorrectSymmetricWithVector(data, vec,M=M)                         
            self.dataDict[i] = ndata
            self.singlesDict[i] = nvec
            self.biasDict[i] = (vec / nvec)
            
    def ultracorrectByTrans(self,names = None):
        "performs iterative correction by trans data only, corrects cis also"
        if names == None: names = self.dataDict.keys()        
        self.transmap = self.chromosomeIndex[:,None] != self.chromosomeIndex[None,:]
        #mat_img(self.transmap)
        for i in names:
            data = self.dataDict[i]
            self.dataDict[i],self.biasDict[i] = numutils.ultracorrectSymmetricByMask(data,self.transmap,50)
            try: self.singlesDict[i] /= self.biasDict[i]
            except: print "bla" 

    def removeZeros(self):
        """removes bins with zero counts
        keeps chromosome starts, ends, etc. consistent"""
        
        s = numpy.sum(self.mask2D,axis = 0) > 0
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
        dicts = [self.trackDict, self.biasDict, self.singlesDict, self.fragsDict]
        for mydict in dicts:
            for key in mydict.keys():
                mydict[key] = mydict[key][s]
        self.chromosomeIndex = self.chromosomeIndex[s]
        self.armIndex = self.armIndex[s]
        self.chromosomeEnds = indices[self.chromosomeEnds]
        self.chromosomeStarts = indices[self.chromosomeStarts]
        self.centromerePositions = indices[self.centromerePositions]
        self.mask2D = self.mask2D[:,s]
        self.mask2D = self.mask2D[s,:]
        return s 


        
    
            
    def doPCA(self):
        """performs PCA on the data
        creates dictionary self.PCA with results"""
        if (self.removedCis == False) or (self.fakedCis == False): 
            print "Cis contacts have not been removed and/or faked."
            print 'Are you sure you want to continue???'
            raw_input("press any button to continue... <-----")            
        self.PCA = {}
        for i in self.dataDict.keys():
            self.PCA[i] = PCA(self.dataDict[i])
        return self.PCA
    def doEig(self):
        """performs eigenvector expansion on the data
        creates dictionary self.EIG with results"""
        if (self.removedCis == False) or (self.fakedCis == False): 
            print "Cis contacts have not been removed and/or faked."
            print 'Are you sure you want to continue???'
            raw_input("press any button to continue... <-----")

        self.EIG = {}
        for i in self.dataDict.keys():
            self.EIG[i] = EIG(self.dataDict[i])
        return self.EIG
    
    
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
    
    def plotScaling(self,name,label = "BLA", color = None):
        "plots scaling of a heatmap,treating arms separately"
        data = self.dataDict[name]
        bins = numutils.logbins(2,self.genome.maximumChromosomeArm/self.resolution,1.17)
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
        print observed 
        print expected 
        values = observed/expected
        bins = numpy.array(bins,float)
        bins2 = 0.5 * (bins[:-1] + bins[1:])
        norm = numpy.sum(values * (bins[1:] - bins[:-1]))
        args = [self.resolution * bins2 / 1000000,values/(1. * norm)]
        if color != None: args.append(color)
        plt.plot(*args, label = label,linewidth = 2)
        


        
    def averageTransMap(self,name = "GM-all", mycmap = "hot_r",chromMax = 22,vmin = None,vmax = None):
        "plots and returns average inter-chromosomal inter-arm map"
        data = self.dataDict[name]
        #data = trunk(data,low = 0,high = 0.0001)         
        avarms = numpy.zeros((80,80))
        avmasks = numpy.zeros((80,80))
        for i in xrange(self.chromosomeCount):
            print i 
            for j in xrange(self.chromosomeCount):
                for k in [-1,1]:
                    for l in [-1,1]:
                        if i == j: continue 
                        cenbeg1 = self.chromosomeStarts[i] + self.genome.centromereStarts[i] / self.resolution
                        cenbeg2 = self.chromosomeStarts[j] + self.genome.centromereStarts[j] / self.resolution
                        cenend1 = self.chromosomeStarts[i] + self.genome.centromereEnds[i] / self.resolution
                        cenend2 = self.chromosomeStarts[j] + self.genome.centromereEnds[j] / self.resolution
                        
                        beg1 = self.chromosomeStarts[i]
                        beg2 = self.chromosomeStarts[j]
                        end1 = self.chromosomeEnds[i]
                        end2 = self.chromosomeEnds[j]
                        if k == 1:
                            bx = cenbeg1-1
                            ex = beg1+1
                            dx = -1
                        else:
                            bx = cenend1+1
                            ex = end1-1
                            dx = 1
                        if l ==1:
                            by = cenbeg2-1
                            ey = beg2+1
                            dy = -1
                        else:
                            by = cenend2+1
                            ey = end2-1
                            dy = 1                        
                        arms = data[bx:ex:dx,by:ey:dy]
                                                
                        mx = numpy.sum(arms, axis = 0)
                        my = numpy.sum(arms, axis = 1)
                        maskx = mx == 0 
                        masky = my == 0
                        mask = (maskx[None,:] + masky[:,None]) == False                    
                        maskf = numpy.array(mask,float)
                        mlenx = (numpy.sum(mask, axis = 0) > 0 ).sum() 
                        mleny = (numpy.sum(mask, axis = 1) > 0 ).sum()
                        if min(mlenx, mleny) < 20: continue
                        #arms = trunk(arms,low = 0.01,high = 0.025)
                        add = numutils.smartZoomOut(arms,avarms.shape)
                        addmask = numutils.smartZoomOut(maskf,avarms.shape)
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
        for chrom1 in xrange(self.chromosomeCount):
            for chrom2 in xrange(self.chromosomeCount):
                for i in self.dataDict.keys():
                    value = self.dataDict[i]
                    submatrix = value[self.chromosomeStarts[chrom1]:self.chromosomeEnds[chrom1],
                                      self.chromosomeStarts[chrom2]:self.chromosomeEnds[chrom2]]
                    masksum = numpy.sum(self.mask2D[self.chromosomeStarts[chrom1]:self.chromosomeEnds[chrom1],
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


    
        
