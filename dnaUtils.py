import plotting 
import numutils 
import tempfile,subprocess
from copy import deepcopy  

from array import array 
import Bio.SeqIO, Bio.SeqUtils
import numpy 
from scipy import weave
from joblib import Memory  

from math import sqrt
import os,cPickle 


#os.chdir("/home/magus/HiC2011/")            

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
        self._parseGapfile(gapfile)  #parsing gap file             
        
        
    def _parseGapfile(self,gapfile):
        "internal: parses .gap file to determine centromere positions"
        if gapfile == None:            
            gapfile = os.path.join(self.genomeFolder,"%s.gap" % self.type )            
        try: 
            gapfile = open(gapfile).readlines()
        except IOError: 
            print "Gap file not found! \n Please provide a link to a gapfile or put a file genome_name.gap in a genome directory"
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
            raise("Chromosome count mismatch between genome and gapfile")
            
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
        self.maximumChromosomeArm =max(lowarms.max() , higharms.max() )
        self.maximumChromosome = max(self.chromosomes)  
        
            
    def loadChromosomeLength(self):
        self.loadSequence()
        return numpy.array([len(self.genome["chr%d" % i] ) for i in xrange(1,self.chromosomeCount+1)])
     

    
    def createMapping(self,resolution,chromosomeExtensionLength = 0 ):
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
        BinnedGC = [[] for _ in xrange(self.chromosomeCount)]
        for  chromNum in xrange(self.chromosomeCount):
            for j in xrange(self.chromosomes[chromNum]/resolution + 1):
                BinnedGC[chromNum].append(self.getGC(chromNum+1,j*resolution,(j+1)*resolution))
                print "Chrom:",chromNum,"bin:",j
        return  BinnedGC
                


def rad2(data):    
    def give_radius_scaling(data):
        N = len(data[0])
        target = int(N**(2/3.))
    
        coms = numpy.cumsum(data,1)
        coms2 = numpy.cumsum(data**2,1)
        def radius_gyration(len2):
            coms2d = (-coms2[:,:-len2]+coms2[:,len2:])/len2
            comsd = ((coms[:,:-len2]-coms[:,len2:])/len2)**2
            diffs = coms2d - comsd
            sums = numpy.sqrt(numpy.sum(diffs,0))
            return numpy.mean(sums)
        return radius_gyration(target)
     
    return give_radius_scaling(data)


     
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

    
#parseDNAseFile("DNAse/wgEncodeDukeDNaseSeqSignalHepg2V2.wig")
#parseDNAseFile("DNAse/wgEncodeDukeDNaseSeqSignalGm12878V2.wig")


def load(filename,center=True):
    "loads xyz polymer file using pure python... slow!"
    f = open(filename, 'r')
    lines = f.readlines()
    N = int(lines[0])            
    datax = array('f', [0.] * N )
    datay = array('f', [0.] * N )
    dataz = array('f', [0.] * N )
    for i in xrange(1, N + 1):
        line = lines[i].split()
        datax[i - 1] = float(line[0])
        datay[i - 1] = float(line[1])
        dataz[i - 1] = float(line[2])
    if center==False:
        return numpy.array([numpy.array(datax), numpy.array(datay), numpy.array(dataz)])
    datax,datay,dataz = numpy.array(datax),numpy.array(datay),numpy.array(dataz)
    diffs = (datax[0:N-1]-datax[1:N])**2+(datay[0:N-1]-datay[1:N])**2+(dataz[0:N-1]-dataz[1:N])**2
    diffs = numpy.abs(1-numpy.sqrt(diffs))
    if numpy.max(diffs) > 0.6:
        print "error, %lf at %s" % (numpy.max(diffs),file)
    datax,datay,dataz = numpy.array(datax),numpy.array(datay),numpy.array(dataz)
    datax -= datax.mean()
    datay -= datay.mean()
    dataz -= dataz.mean()
    return numpy.array([datax,datay,dataz])


def Cload(filename,center = True):
    "fast polymer loader using weave.inline"
    f = open(filename, 'r')
    line = f.readline()
    N = int(line)
    ret = numpy.zeros((3,N),order = "C",dtype = numpy.double)
    code = """
    #line 85 "binary_search.py"
    using namespace std;
    FILE *in;
    const char * myfile = filename.c_str(); 
    in = fopen(myfile,"r");
    int dummy;
    dummy = fscanf(in,"%d",&dummy);    
    for (int i = 0; i < N; i++)
    {
    dummy = fscanf(in,"%lf %lf %lf ",&ret[i],&ret[i + N],&ret[i + 2 * N]);
    if (dummy < 3){printf("Error in the file!!!");throw -1;}
    } 
    """
    support = """
    #include <math.h>  
    """
    weave.inline(code, ['filename', 'N' , 'ret' ], extra_compile_args=['-march=native -malign-double'],support_code =support )
    return ret
    

def intload(filename,center = False):
    "dedicated weave.inline for loading int polymers"
    f = open(filename, 'r')
    line = f.readline()
    f.close()
    N = int(line)
    ret = numpy.zeros((3,N),order = "C",dtype = numpy.int)
    code = """
    #line 85 "binary_search.py"
    using namespace std;
    FILE *in;
    const char * myfile = filename.c_str(); 
    in = fopen(myfile,"r");
    int dummy;
    dummy = fscanf(in,"%d",&dummy);    
    for (int i = 0; i < N; i++)
    {
    dummy = fscanf(in,"%ld %ld %ld ",&ret[i],&ret[i + N],&ret[i + 2 * N]);
    if (dummy < 3){printf("Error in the file!!!");throw -1;}        
    } 
    fclose(in);
    """
    support = """
    #include <math.h>  
    """
    weave.inline(code, ['filename', 'N' , 'ret' ], extra_compile_args=['-march=native -malign-double'],support_code =support )
    return ret


def giveIntContacts(data):    
    "give all contacts of a polymer on a cubic lattice"
    data = numpy.transpose(data)
    data -= numpy.min(data)
    M = numpy.max(data)+1
    if M>1500: raise ValueError("Polymer is to big, can't create bounding box!")
    N = len(data)
    tocheck = numpy.zeros(M*M*M,dtype = numpy.int32) - 1
    tocheck[data[:,0] + data[:,1] * M + data[:,2] * M * M] = numpy.arange(N,dtype = numpy.int32)
    tocheck.shape = (M,M,M)
    contacts1 = numpy.concatenate([tocheck[1:,:,:].ravel(),tocheck[:,1:,:].ravel(),tocheck[:,:,1:].ravel()]) 
    contacts2 = numpy.concatenate([tocheck[:-1,:,:].ravel(),tocheck[:,:-1,:].ravel(),tocheck[:,:,:-1].ravel()])        
    mask = (contacts1 != -1) * (contacts2 != -1)
    contacts1 = contacts1[mask]
    contacts2 = contacts2[mask]
    contacts3 = numpy.minimum(contacts1,contacts2)
    contacts4 = numpy.maximum(contacts1,contacts2)    
    return numpy.concatenate([contacts3[:,None],contacts4[:,None]],1)
    
 
def give_contacts(data,cutoff=1.7,maxContacts = 40 ):
    "returns contacts of a current structure with a current cutoff"    
    points = numpy.zeros((maxContacts*len(data[0]),2),int,order = "C")    
    data = numpy.asarray(data,float,order = "C")
    if len(data.shape) != 2: raise ValueError("Wrong dimensions of data")
    if len(data[0]) == 3: data = data.T
    if len(data) != 3: raise ValueError("Wrong size of data: %s,%s" % data.shape) 
    N = len(data[0])
    N #Eclipse warning remover 
    code = """
    #line 50 "binary_search.py"
        using namespace std; 
        int counter  = 0;
          double  d;
        for (int i=0;i<N;i++)
            {
            for (int j=i+1;j<N;j++)
                {
                d = dist(data,N,i,j); 
                if (d<= cutoff) 
                    {
                    int t = j - i;                     
                    points[2*counter] = i;
                    points[2*counter + 1 ] = j; 
                    counter ++ ;
                    //cout << i << " "   << j  << " "  << d <<  " " << counter << endl;   
                
                    }
                else if (d>4)
                    {
                    j +=  d-4; 
                    }
            }
            }
       """
    support = """
        #include <math.h>  
           double  dist(double* data, int N, int i,int j)
        {
        return sqrt(pow((data[0*N + i]-data[0*N + j]),2.) + pow((data[N+i] - data[N+j]),2.) + pow((data[2*N+i] - data[2*N+j]),2.)); 
        }
       """
    weave.inline(code, ['data', 'N' , 'cutoff' , 'points'], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
    k = numpy.max(numpy.nonzero(points[:,1]))
    return points[:k,:]


 
def give_contacts_arb(data,cutoff=1.4,maxContacts = 100):
    "returns contacts of a current structure with a current cutoff"
    points = numpy.zeros((maxContacts*len(data[0]),2),int,order = "C")    
    data = numpy.asarray(data,float,order = "C")
    if len(data.shape) != 2: raise ValueError("Wrong dimensions of data")
    if len(data[0]) == 3: data = data.T
    if len(data) != 3: raise ValueError("Wrong size of data: %s,%s" % data.shape) 
    N = len(data[0])
    N #Eclipse warning remover  
    code = """
    #line 50 "binary_search.py"
        using namespace std; 
    
        
        int counter  = 0;
         
        double  d;
        for (int i=0;i<N;i++)
            {
            for (int j=i+1;j<N;j++)
                {
                 
                if (dist(data,N,i,j) <= CUTOFF2) 
                    {
                    int t = j - i;                     
                    points[2*counter] = i;
                    points[2*counter + 1 ] = j; 
                    counter ++ ;
                    //cout << i << " "   << j  << " "  << d <<  " " << counter << endl;                   
                    }
            }
            }
       """
    support = """
        #include <math.h>
        #define CUTOFF2  %lf
        double sq(double x) {return x*x;} 
           double  dist(double* data, int N, int i,int j)
        {
         return sq(data[0*N + i]-data[0*N + j]) + sq(data[N+i] - data[N+j]) + sq(data[2*N+i] - data[2*N+j]);
        }
       """ % (cutoff*cutoff) 
    #
    weave.inline(code, ['data', 'N' , 'cutoff' , 'points'], extra_compile_args=['-march=native -malign-double'],support_code =support )
    k = numpy.max(numpy.nonzero(points[:,1]))
    return points[:k,:]


def give_distance_map(data,size = 1000):
    "returns  a distance map of a polymer with a given size"    
    toret = numpy.zeros((size,size),float,order = "C")
    data = numpy.asarray(data,order = "C")
    N = len(data[0])
    N #Eclipse warning remover 
    code = """
    #line 50 "binary_search.py"
        using namespace std; 
        int counter  = 0;
          double  d;
        for (int i=0;i<size;i++)
            {
            for (int j=0;j<size;j++)
                {
                d = dist(data,N,N*i/size,N*j/size);
                toret[size*i+j] = d;                 
            }
            }
       """ 
    support = """
        #include <math.h>  
           double  dist(double* data, int N, int i,int j)
        {
        return sqrt(pow((data[0*N + i]-data[0*N + j]),2.) + pow((data[N+i] - data[N+j]),2.) + pow((data[2*N+i] - data[2*N+j]),2.)); 
        }
       """
    weave.inline(code, ['data', 'N','size', 'toret'], extra_compile_args=['-march=native -malign-double'],support_code =support )
    return toret


def findSimplifiedPolymer(data):
    """a weave.inline wrapper for polymer simplification code
    Calculates a simplified shape of a polymer ring"""
    
    if len(data) != 3: data = numpy.transpose(data)
    if len(data) != 3: raise ValueError("Wrong dimensions of data")        
    datax = numpy.array(data[0],float, order = "C")
    datay = numpy.array(data[1], float,order = "C")
    dataz = numpy.array(data[2], float,order = "C")
    N = len(datax)
    ret = numpy.array([1])    
    datax,datay,dataz,N ##eclipse warning removal 
    code = r"""
    #line 290 "binary_search.py"
    int M = 0;
    int sum = 0;
    int t=0,s=0,k=0;
    int turn=0;
    bool breakflag;
    double maxdist=0;
    int a;
    position=vector<point>(N);
    newposition=vector<point>(N);

    for (i=0;i<N;i++)
    {
    position[i].x = datax[i] +  0.00001*(rand()%1000);
    position[i].y = datay[i] +0.00001*(rand()%1000);
    position[i].z  = dataz[i] +  0.00001*(rand()%1000);    
    }
    todelete = vector <int> (N);
    for (i=0;i<N;i++) todelete[i] == -2;    
    while (true)
        {
        maxdist = 0; 
        for (i=0;i<N-1;i++)
        {
        if (dist(i,i+1) > maxdist) {maxdist = dist(i,i+1);}        
        }
        //printf ("maxdist = %lf\n",maxdist);
        turn++;
        M=0;
        for (i=0;i<N;i++) todelete[i] = -2;
        for (int j=1;j<N-1;j++)  //going over all elements trying to delete
            {

            breakflag = false; //by default we delete thing
            for (k=0;k<N-1;k++)  //going over all triangles to check
                {
                double dd = dist(j,k);
                if (dd  < 3 * maxdist)
                {
        
                if (k < j-2 || k > j+1)
                    {
                    sum = intersect(position[j-1],position[j],position[j+1],position[k],position[k+1]);
                    if (sum!=0)
                        {
                        //printf("intersection at %d,%d\n",j,k);
                        breakflag = true; //keeping thing
                        break;
                        }
                    }
                }
                else k+= abs((int)((float)dd/(float)maxdist )- 3);
                }
            if (breakflag ==false)
            {
            todelete[M++] = j;
            position[j] = (position[j-1] + position[j+1])* 0.5;
            //printf("%d will be deleted at %d\n",j,k);
            j++;
            //break;
            }
            }
        t = 0;//pointer for todelete
        s = 0;//pointer for newposition
        if (M==0)
            {
            break;
            }
        for (int j=0;j<N;j++)
            {            
            if (todelete[t] == j)
                {
                t++;
                continue;
                }
            else
                {
                newposition[s++] = position[j];
                }            
            }
        N = s;
        M = 0;
        t = 0;
        position = newposition;            
        }
    ret[0] = N;
    """
    support = r"""
#line 400 "binary_search.py"    
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <ctime>
#include <omp.h>
#include <stdio.h>
using namespace std;
struct point{
    double x,y,z;
    point operator + (const point &p) const {
        return (point) {x+p.x, y+p.y, z+p.z};
    }
    point operator - (const point &p) const {
        return (point) {x-p.x, y-p.y, z-p.z};
    }
/* cross product */
    point operator * (const point &p) const {
        return (point) {y*p.z - z*p.y,
                        z*p.x - x*p.z,
                        x*p.y - y*p.x};
    }
    point operator * (const double &d) const {
        return (point) {d*x, d*y, d*z};
    }

    point operator / (const double &d) const {
        return (point) {x/d, y/d, z/d};
    }
};

vector <point> position;
vector <point> newposition;
vector <int> todelete;
int N;
int i; 
double dist(int i,int j);
double dotProduct(point a,point b);
int intersect(point t1,point t2,point t3,point r1,point r2);

inline double sqr(double x){
    return x*x;
}
inline double dist(int i,int j){
    return sqrt(dotProduct((position[i]-position[j]),(position[i]-position[j])));
}

inline double dist(point a,point b){
    return sqr(a.x-b.x)+sqr(a.y-b.y)+sqr(a.z-b.z);
}

inline double dotProduct(point a,point b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

int intersect(point t1,point t2,point t3,point r1,point r2)
{
point A,B,C,D,n;
int r;
double det,t,u,v,c1,d1,d2,d3;
B = t2 - t1;
C = t3 - t1;
D = r2 - t1;
A = r2 - r1;

d1 = (B.y*C.z-C.y*B.z);
d2 = (B.x*C.z-B.z*C.x);
d3 = (B.x*C.y-C.x*B.y);
det = A.x*d1-A.y*d2+A.z*d3;
if (det == 0) return 0;
if (det >0){
t = D.x*d1-D.y*d2+D.z*d3;
if (t<0 || t>det) return 0;
u = A.x*(D.y*C.z-C.y*D.z)-A.y*(D.x*C.z-D.z*C.x)+A.z*(D.x*C.y-C.x*D.y);
if (u<0 || u>det) return 0;
v = A.x*(B.y*D.z-D.y*B.z)-A.y*(B.x*D.z-B.z*D.x)+A.z*(B.x*D.y-D.x*B.y);
if (v<0 || v>det || (u+v)>det) return 0;
//printf("\n%lf,%lf,%lf, ",t/det,u/det,v/det);
n = B*C;
c1 = dotProduct(r1-t1,n);
if (c1>0) return 1;
else return -1;
}
else{
t = D.x*d1-D.y*d2+D.z*d3;
if (t>0 || t<det) return 0;
u = A.x*(D.y*C.z-C.y*D.z)-A.y*(D.x*C.z-D.z*C.x)+A.z*(D.x*C.y-C.x*D.y);
if (u>0 || u<det) return 0;
v = A.x*(B.y*D.z-D.y*B.z)-A.y*(B.x*D.z-B.z*D.x)+A.z*(B.x*D.y-D.x*B.y);
if (v>0 || v<det || (u+v)<det) return 0;
//printf("\n%lf,%lf,%lf, ",t/det,u/det,v/det);
n = B*C;
c1 = dotProduct(r1-t1,n);
if (c1>0) return 1;
else return -1;
}
}
//DNA conformation
"""    
    weave.inline(code, ['datax', 'datay','dataz', 'N','ret'], extra_compile_args=['-malign-double'],support_code =support )
    return ret[0]


def rescale_points(points, res = 100):
    "converts array of contacts to the reduced resolution contact map"   
    a = numpy.histogram2d(points[:,0],points[:,1],res)[0]
    a = a + numpy.transpose(a)
    return a 

 
def rescaled_map(data,res,cutoff = 1.4 ):
    "calculates a rescaled contact map of a structure"    
    t = give_contacts(data,cutoff)
    return rescale_points(t,res) 
    
def pure_map(data,N=3,cutoff=1.4):
    "calculates an all-by-all contact map of a structure"
    t = give_contacts(data,cutoff)
    N = len(data[0])
    contacts = numpy.zeros((N,N),int)
    for i in t:
        contacts[i[0],i[1]] = 1;
        contacts[i[1],i[0]] = 1;
    return contacts  
        

def observedOverExpected(matrix):
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


def cool_trunk(data):
    "somehow trunkates the globule so that the ends are on the surface"
    CUTOFF = 0.7
    CUTOFF2 = 3    
    datax, datay, dataz = data[0], data[1], data[2]
    N = len(data[0])
    found1, found2 = 0, 0
    def dist(i):
        return sqrt(datax[i] ** 2 + datay[i] ** 2 + dataz[i] ** 2)
    def sqdist(i, j):
        return sqrt((datax[i] - datax[j]) ** 2 + (datay[i] - datay[j]) ** 2 + (dataz[i] - dataz[j]) ** 2)
    def sqdist2(i, j, scale):
        return sqrt((scale * datax[i] - datax[j]) ** 2 + (scale * datay[i] - datay[j]) ** 2 + (scale * dataz[i] - dataz[j]) ** 2)
    breakflag = 0
    escapeflag = 0
    for i in xrange(N):
        exitflag = 0
        pace = 0.25 / dist(i)
        for j in numpy.arange(1, 10, pace):
            #print j
            escapeflag = 1
            for k in xrange(N):
                if k == i: continue
                escapeflag, breakflag
                
                if 0.001 < sqdist2(i, k, j) < CUTOFF:
                    breakflag = 1
                    print "breaking at", k
                    break
                if 0.001 < sqdist2(i, k, j) < CUTOFF2:                    
                    escapeflag = 0
            
            if breakflag == 1:
                breakflag = 0
                print i, dist(i), j
                break
            if escapeflag == 1:
                print i, dist(i)
                found1 = i
                exitflag
                exitflag = 1
                break
        if exitflag == 1: break
    for i in xrange(N - 1, 0, - 1):
        exitflag = 0
        pace = 0.25 / dist(i)
        for j in numpy.arange(1, 10, pace):
            #print j
            escapeflag = 1
            for k in xrange(N - 1, 0, - 1):
                if k == i: continue
                escapeflag, breakflag
                
                if 0.001 < sqdist2(i, k, j) < CUTOFF:
                    breakflag = 1
                    print "breaking at", k
                    break
                if 0.001 < sqdist2(i, k, j) < CUTOFF2:                    
                    escapeflag = 0
            
            if breakflag == 1:
                breakflag = 0
                print i, dist(i), j
                break
            if escapeflag == 1:
                print i, dist(i)
                found2 = i
                exitflag
                exitflag = 1
                break
        if exitflag == 1: break
        
    f1 = found1
    f2 = found2
    datax[f1 - 1], datay[f1 - 1], dataz[f1 - 1] = datax[f1] * 2, datay[f1] * 2, dataz[f1] * 2    
    datax[f2 + 1], datay[f2 + 1], dataz[f2 + 1] = datax[f2] * 2, datay[f2] * 2, dataz[f2] * 2
    a = (datax[f1 - 1], datay[f1 - 1], dataz[f1 - 1])
    b = (datax[f2 + 1], datay[f2 + 1], dataz[f2 + 1])
    c = (a[0] + b[0], a[1] + b[1], a[2] + b[2])
    absc = sqrt(c[0] ** 2 + c[1] ** 2 + c[2] ** 2)
    c = (c[0] * 40 / absc, c[1] * 40 / absc, c[2] * 40 / absc)
    datax[f1 - 2], datay[f1 - 2], dataz[f1 - 2] = c
    datax[f2 + 2], datay[f2 + 2], dataz[f2 + 2] = c        
    return [datax[f1 - 2:f2 + 3], datay[f1 - 2:f2 + 3], dataz[f1 - 2:f2 + 3]]
            
