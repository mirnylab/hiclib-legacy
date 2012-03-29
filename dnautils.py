import numutils 
import tempfile,subprocess
from array import array 
import Bio.SeqIO, Bio.SeqUtils, Bio.Restriction
Bio.Restriction
import numpy 
from scipy import weave
from joblib import Memory  
from math import sqrt
import os,cPickle 
from numutils import arrayInArray


class Genome():
    def __init__(self,genomeFolder, genomeName = None, gapfile = None ):
        """Loads a genome and calculates certain genomic parameters"""        
        #setting up genome name  and folder
        self.genomeFolder = genomeFolder                 
        if genomeName == None: 
            folderName = os.path.split(genomeFolder)[1]
            self.type = folderName
            print "genome name inferred from genome folder name : %s" % self.type
        else: 
            self.type = genomeName  

        #defining memorized functions 
        mymem = Memory(cachedir = self.genomeFolder)
        self.loadChromosomeLength = mymem.cache(self.loadChromosomeLength)
        self.getBinnedGCContent = mymem.cache(self.getBinnedGCContent)
        self.getRsitesRfrags = mymem.cache(self.getRsitesRfrags)
        
        #detecting number of chromosomes        
        names = os.listdir(self.genomeFolder)
        names = [i for i in names if ((i[-2:] == "fa") and (i[:3] == "chr"))]
        if len(names) == 0: raise("No Genome files found")
        chrs = []
        for i in names: 
            try:                     
                chrs.append(int(i[3:-3]))
            except ValueError:
                pass
        self.chromosomeCount = max(chrs)+1        
        self.chromosomes = self.loadChromosomeLength()   #loading cached chromosome length
        self.chromosomeLength = self.chromosomes   #backwards compatibility - chromsomeLength should be preferred
        self.maxChromLem = max(self.chromosomes)  
        self.fragIDmult = self.maxChromLem + 1000   #to be used when calculating fragment IDs for HiC  
        self._parseGapfile(gapfile)  #parsing gap file                     
        
    def _parseGapfile(self,gapfile):
        "internal: parses .gap file to determine centromere positions"
        if gapfile == None:            
            gapfile = os.path.join(self.genomeFolder,"%s.gap" % self.type )            
        try: 
            gapfile = open(gapfile).readlines()
        except IOError: 
            print """Gap file not found! \n 
            Please provide a link to a gapfile or put a file genome_name.gap in a genome directory.
            If there is no gapfile, create a gapfile with "centromere" records for each chromosome.  
            """
            exit() 
        gaps = [i.split() for i in gapfile]
        centromeres = [i for i in gaps if i[7] == 'centromere']
        def chromnum(x):  #calculates number of chromosomes, X = 0 
            a = x[3:]
            if (a == "X") or (a == "x"):
                return 0
            else:
                try:
                    return int(a)
                except:return -1
        chromnums = [chromnum(i[1]) for i in centromeres]
        cm = max(chromnums)
        if cm+1 != self.chromosomeCount: 
            raise ValueError("Chromosome count mismatch between genome and gapfile")            
        for i in xrange(len(chromnums)): 
            if chromnums[i] == 0: chromnums[i] = cm + 1

        self.centromereStarts = numpy.zeros(cm + 1,int)
        self.centromereEnds = numpy.zeros(cm+1,int)
        for i in xrange(len(chromnums)):
            if chromnums[i] != -1:
                chrom = chromnums[i]
                self.centromereStarts[chrom-1] = int(centromeres[i][2])
                self.centromereEnds[chrom-1] = int(centromeres[i][3])
        self.centromeres = (self.centromereStarts + self.centromereEnds) / 2
        lowarms = numpy.array(self.centromereStarts)
        higharms = numpy.array(self.chromosomes) - numpy.array(self.centromereEnds)
        self.maximumChromosomeArm = max(lowarms.max() , higharms.max() )
        self.maximumChromosome = max(self.chromosomes)          
            
    def loadChromosomeLength(self):
        #Memorized function
        self.loadSequence()
        return numpy.array([len(self.genome["chr%d" % i] ) for i in xrange(1,self.chromosomeCount+1)])     
    
    def createMapping(self,resolution,chromosomeExtensionLength = 0 ):
        """Calculates chromosome start/end positions for whole-genome datasets"""
        self.resolution = resolution
        self.chromosomeExtensionLength = chromosomeExtensionLength
        self.chromosomeSizes = numpy.array([i/self.resolution + 1 + self.chromosomeExtensionLength for i in self.chromosomes])
        self.realChromosomeSizes = numpy.array([i/self.resolution + 1 for i in self.chromosomes])
        self.chromosomeStarts = numpy.r_[0,numpy.cumsum(self.chromosomeSizes)[:-1]]
        self.centromerePositions = self.chromosomeStarts + numpy.array([i/self.resolution for i in self.centromeres],int)         
        
        self.chromosomeEnds = numpy.cumsum(self.chromosomeSizes)
        self.realChromosomeEnds = self.chromosomeStarts + self.realChromosomeSizes        
        self.N = self.chromosomeEnds[-1]
        self.chromosomeCount = len(self.chromosomeEnds)
        self.chromosomeIndex = numpy.zeros(self.N,int)
        self.positionIndex = numpy.zeros(self.N,int)        
        for i in xrange(self.chromosomeCount):
            self.chromosomeIndex[self.chromosomeStarts[i]:self.chromosomeEnds[i]] = i
            self.positionIndex[self.chromosomeStarts[i]:self.chromosomeEnds[i]] = numpy.arange(-self.chromosomeStarts[i]+self.chromosomeEnds[i]) * self.resolution
                    
    def loadSequence(self):
        "loads genomic sequence if it wasn't loaded."
        if hasattr(self,"genome"): return 
        def loadGenome(folder):        
            names = ["chr%d" % i for  i in xrange(1,self.chromosomeCount)] + ["chrX"]
            genome = {}
            for i in names:
                genome[i] = Bio.SeqIO.read(open(os.path.join(self.genomeFolder,i+".fa")), 'fasta')
            genome["chr%d" % self.chromosomeCount] = genome.pop("chrX")
            print genome.keys()                
            return genome        
        self.genome = loadGenome("data/%s" % self.type)                    
        
    def getSequence(self,chromosome,start,end):
        if not hasattr(self,"genome"): self.loadSequence()        
        return self.genome["chr%d" % chromosome][start:end]
    
    def getGC(self,chromosome,start,end):
        "returns GC content of a region"
        seq = self.getSequence(chromosome,start,end)
        return Bio.SeqUtils.GC(seq.seq)
    
    def getBinnedGCContent(self,resolution):
        #Memorized function
        BinnedGC = [[] for _ in xrange(self.chromosomeCount)]
        for  chromNum in xrange(self.chromosomeCount):
            for j in xrange(self.chromosomes[chromNum]/resolution + 1):
                BinnedGC[chromNum].append(self.getGC(chromNum+1,j*resolution,(j+1)*resolution))
                print "Chrom:",chromNum,"bin:",j
        return  BinnedGC
    def getRsitesRfrags(self,enzymeName):
        """returns: tuple(rsiteMap,rfragMap) 
        Finds restriction sites and mids of rfrags for a given enzyme
        Note that there is one extra rsite at beginning and end of chromosome
        Note that there are more rsites than rfrags (by 1)"""
        
        #Memorized function
        self.loadSequence()
        enzymeSearchFunc = eval('Bio.Restriction.%s.search' % enzymeName)
        rsiteMap = {}
        rfragMap = {}        
        for i in self.genome.keys():
            rsites = numpy.r_[-1 ,(enzymeSearchFunc(self.genome[i].seq)),len(self.genome[i].seq)-1] + 1   #+1 is a convention 
            rfrags = (rsites[:-1] + rsites[1:]) / 2
            rsiteMap[i] = rsites
            rfragMap[i] = rfrags          
        return rsiteMap,rfragMap 
    
    def _calculateRsiteIDs(self,enzymeName):
        "Calculates rsite/rfrag positions and IDs for a given enzyme name and memorizes them"
        rsiteMap, rfragMap = self.getRsitesRfrags(enzymeName)
        #Now truncating one "fake" rsite at the end of each chr. so that number of rsites matches number of rfrags         
        for i in rsiteMap.keys(): rsiteMap[i] = rsiteMap[i][:-1]
        self.rsiteMap = rsiteMap 
        self.rfragMap  = rfragMap         
        rsiteIDs = [self.rsiteMap["chr%d" % chrom] + (chrom-1) * self.fragIDmult for chrom in xrange(1,self.chromosomeCount+1)]        
        rsiteChroms = [numpy.ones(len(self.rsiteMap["chr%d" % chrom]),int) *  chrom for chrom in xrange(1,self.chromosomeCount+1)]
        rfragIDs = [self.rfragMap["chr%d" % chrom] + (chrom-1) * self.fragIDmult for chrom in xrange(1,self.chromosomeCount+1)]
        self.rsiteIDs = numpy.concatenate(rsiteIDs)
        self.rsiteChroms = numpy.concatenate(rsiteChroms)
        self.rfragIDs = numpy.concatenate(rfragIDs)
        assert (len(self.rsiteIDs) == len(self.rfragIDs))
        
        
    
    def getFragmentDistance(self,fragments1,fragments2,enzymeName):
        "returns distance between fragments in... fragments. (neighbors = 1, etc. )"
        if not hasattr(self,"rsiteIDs"): self._calculateRsiteIDs(enzymeName)
        frag1ind = numpy.searchsorted(self.rsiteIDs,fragments1)-1
        frag2ind = numpy.searchsorted(self.rsiteIDs,fragments2)-1
        assert (fragments1[::100] == self.rfragIDs[frag1ind[::100]]).all()   #fragments are consistent
        distance = numpy.abs(frag1ind - frag2ind)
        del frag1ind,frag2ind
        ch1 = fragments1 / self.fragIDmult
        ch2 = fragments2 / self.fragIDmult
        distance[ch1 != ch2] = 1000000
        return distance
    
    def getPairsLessThanDistance(self,fragments1,fragments2,cutoffDistance,enzymeName):
        "returns all possible pairs (fragment1,fragment2) with fragment distance less-or-equal than cutoff"
        if not hasattr(self,"rsiteIDs"): self._calculateRsiteIDs(enzymeName)
        f1ID = numpy.searchsorted(self.rsiteIDs,fragments1)-1 
        f2ID = numpy.searchsorted(self.rsiteIDs,fragments2)-1
        assert (fragments1[::40] - self.rfragIDs[f1ID[::40]]).sum() == 0 
        assert (fragments2[::40] - self.rfragIDs[f2ID[::40]]).sum() == 0

        fragment2Candidates = numpy.concatenate([f1ID + i for i in (range(-cutoffDistance,0) + range(1,cutoffDistance+1))])        
        fragment1Candidates = numpy.concatenate([f1ID for i in (range(-cutoffDistance,0) + range(1,cutoffDistance+1))])                 
        mask = arrayInArray(fragment2Candidates,f2ID) 
         
        
        fragment2Real = fragment2Candidates[mask]
        fragment1Real = fragment1Candidates[mask]
        
        return  (self.rfragIDs[fragment1Real],self.rfragIDs[fragment2Real])
        
        
        
         
        
         
    



def liftOver(x):
    "chromosomes go from 0 to 22, 22 being 'X' chromosome"
     
    chain_file = "genomes/hg18ToHg19.over.chain"
    liftOverFile = "genomes/liftOver"    
    tmp_input = tempfile.NamedTemporaryFile()
    tmp_output = tempfile.NamedTemporaryFile()
    tmp_error = tempfile.NamedTemporaryFile()
    a = ["chr%d" % (i+1) for i in xrange(22) ] + ["chrX"]
    listOfChromosomeIntervals = [(a[int(i[0])],int(i[1]),int(i[2]),i[3]) for i in x]

    for row in listOfChromosomeIntervals:
        tmp_input.write(" ".join(str(x) for x in row))
        tmp_input.write("\n")
    tmp_input.write("\n")    
    tmp_input.flush()#it's magic ;)
    
    cmd = [liftOverFile, tmp_input.name, chain_file, 
           tmp_output.name, tmp_error.name]
    
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
    p.communicate()    
    tmp_output.seek(0, os.SEEK_SET)
    
    res = []

    def chromosome(i):
        try:
            return int(i[3:]) - 1
        except:
            if i[3] == "X":
                return 22
            else:
                print "strange chromosome found!!!: ",i
                return None         

    for row in tmp_output:
        row = row.strip().split("\t")
        row[0] = chromosome(row[0])        
        row[1] = int(row[1])
        row[2] = int(row[2])
        if row[0] != None: res.append(tuple(row))
    tmp_error.seek(0, os.SEEK_SET)
    d = tmp_error.read()  
    if len(d) > 10: print d[:1000],"these were first 1000 symbols of error output\n"
    tmp_input.close()
    tmp_output.close()
    tmp_error.close()
    return res

#print liftOver([(i,2000000,2000001,"bla") for i in xrange(23)])     

def parseLamin(filename):
    laminFile  = filename
    laminFile #Eclipse warning remover 
    chroms = numpy.zeros(4000000,int)
    positions = numpy.zeros(4000000,int)
    values = numpy.zeros(4000000,float)
    end = numpy.zeros(1,int)
    code = r"""
    #line 16 "binary_search.py"
    using namespace std;
    int chrom;    
    
    char line[260];
    const char * filename = laminFile.c_str();    
    FILE *myfile; 
    myfile = fopen(filename,"r");
    int breakflag = 0;
    fgets(line, 260, myfile);
    fgets(line, 260, myfile);
    
    int point = 0; 
    while (fgets(line, 260, myfile) != NULL)    
    {
      if (line[0] == 'v')
          {
          for (int j = 0;j<strlen(line);j++)
          {
          if (line[j] == 'X') 
              {              
              cout << "terminated with string " << line << endl;
              breakflag = 1; 
              break;
              }
          }
          if (breakflag == 1) break;                    
          sscanf(line,"variableStep chrom=chr%d",&chrom);
          cout << chrom << endl;          
          continue; 
          }
       
      double t;
      int pos;        
      sscanf(line,"%d\t%lf",&pos,&t);
      chroms[point] = chrom;
      positions[point] = pos;
      values[point] = t;
      point++;
    }
    end[0] = point;
    """
    support = """
    #include <math.h>
    #include <iostream>
    #include <fstream>      
    """
    
    weave.inline(code, ['laminFile',"chroms","positions","values","end" ], extra_compile_args=['-march=native -malign-double'],support_code =support )
    end = end[0]
    chroms = chroms[:end]
    positions = positions[:end]
    values = values[:end]
    chroms = numpy.array(chroms, float)
    positions = numpy.array(positions,float)
    print end 
    a = numpy.array([chroms,positions,values], dtype = float)
    cPickle.dump(a,open("lamin",'wb'),-1)

def liftOverLamin():
    a = cPickle.load(open("lamin")).T
    
    b = numpy.zeros((len(a),4),float)
    b[:,0] = a[:,0]-1
    b[:,1] = a[:,1]
    b[:,2] = a[:,1]+1
    b[:,3] = a[:,2]
    c = liftOver(b)
    d = numpy.array(c,'float')
    #cPickle.dump(d,open("timing/"))
    d[:,0] += 1
    cPickle.dump(d[:,numpy.r_[0,1,3]].T,open("lamin_HG19",'wb'),-1)
    
#liftOverLamin() 
    
#parseLamin("hg18.laminB1.txt")


def parseDNAseFile(DNAseFile):
    "code used for parsing fixedStep wig files"
    data = numpy.zeros(250000 * 25,float) 
    code = r"""
    #line 14 "binary_search.py"
    using namespace std;
    int chrom;    
    int pos;
    char line[60];
    const char * filename = DNAseFile.c_str();    
    FILE *myfile; 
    myfile = fopen(filename,"r");
    int breakflag = 0;          
    while (fgets(line, 60, myfile) != NULL)    
    {           
      if (line[0] == 'f')
          {
          for (int j = 0;j<strlen(line);j++)
          {
          if (line[j] == 'X') 
              {              
              cout << "terminated with string " << line << endl;
              breakflag = 1; 
              break;
              }
          }
          if (breakflag == 1) break;                    
          sscanf(line,"fixedStep chrom=chr%d start=%d step=1",&chrom,&pos);
          cout << chrom << endl;
          cout << pos << endl;
          continue; 
          }
       
      double t; 
      sscanf(line,"%lf",&t);                     
      data[250000 * (chrom - 1) + pos / 1000] += t;
      pos++;       
    }
    """
    support = """
    #include <math.h>
    #include <iostream>
    #include <fstream>      
    """
    weave.inline(code, ['DNAseFile',"data" ], extra_compile_args=['-march=native -malign-double'],support_code =support )
    cPickle.dump(data,open(file + ".dat",'wb'),-1)


def observedOverExpected(matrix):
    "Calculates observedOverExpected of any contact map"
    data = numpy.array(matrix, order = "C")
    N = data.shape[0]
    bins = numutils.logbins(1,N,1.2)
    bins = [(0,1)] + [(bins[i],bins[i+1]) for i in xrange(len(bins)-1)]
    bins = numpy.array(bins,order = "C")
    M = len(bins)
    M #Eclipse warning remover
    code = r"""
    #line 50 "binary_search.py"
    using namespace std;
    for (int bin = 0; bin < M; bin++)
    {
        int start = bins[2 * bin];
        int end = bins[2 * bin + 1];
        
        double ss = 0 ;
        int count   = 0 ;  
        for (int offset = start; offset < end; offset ++)
        {
            for (int j = 0; j < N - offset; j++)
            {
                ss += data[(offset + j) * N + j];
                count += 1;                            
            }
        }
        double meanss = ss / count;
        printf("%lf\n",meanss); 
        for (int offset = start; offset < end; offset ++)
        {
            for (int j = 0; j < N - offset; j++)
            {
                data[(offset + j) * N + j] /= meanss;                                             
                if (offset > 0) {data[(offset + j)  + j*N] /= meanss;}
            }
        }
    }
    
    """
    support = """
    #include <math.h>  
    """
    weave.inline(code, ['data', 'bins' , 'N' ,'M'], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
    return data
    
#mat_img(observedOverExpected(rescaled_map(Cload("/home/magus/evo/trajectories/globule_creation/32000_RW/crumpled1.dat"),500)))

def save(filename, data,doint=False):
    f = open(filename, 'w')
    f.write(str(len(data[0])) + "\n")
    for i in xrange(len(data[0])):
        if (doint == False):
            for j in xrange(len(data)):
                f.write(str(data[j][i]) + " ")
            f.write("\n")
        else:
            for j in xrange(len(data)):
                f.write(str(int(data[j][i])) + " ")
            f.write("\n")            
    f.close()
