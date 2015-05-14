from mirnylib.h5dict import h5dict
from mirnylib.numutils import completeIC
import matplotlib.pyplot as plt 
import numpy as np 

filename = "MyFolder/byChromosomeHiCDataset" 
#this only works with Hi-C dataset saved by chromosome

dataset = h5dict(filename,'r')  #open in the "read" mode

#a bit of a weird way to find chromosome number
keys = dataset.keys()
cisKeys = [i for i in keys if len(set(i.split()))==1] #extract keys of the type "a a"
numChromosomes = len(cisKeys)


for chromosome in range(numChromosomes):

    chromosomeHeatmap = dataset["{0} {0}".format(chromosome)] #extracting cis heatmap
    

    #This line executes proper Iterative Correction, which accounts for regions with low coverage
    #It only works for cis (symmetric) heatmaps. 
    correctedHeatmap,bias = completeIC(chromosomeHeatmap,returnBias=True)
    
    #if you want to see log-heatmap, uncomment below
    #plt.imshow(np.log(correctedHeatmap))
    #plt.show()
    
    """
    Your code here
    """
    
    #Below is example of how to save data to txt
    #np.savetxt("corrected_chr{0}.txt.gz".format(chromosome),correctedHeatmap)  #example of saving chromosomes


