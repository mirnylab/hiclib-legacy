from mirnylib import h5dict, genome

"""
Here you need to edit the path to the genome files (see mirnylib.genome file or online help for details)
Other components will automatically 

Note that names of fasta files should match exactly to the names of sequences in a fasta file
And that each chromosome has exactly the same sequence
"""
allGenomes = {}

def getGenome(name): 
    if name in allGenomes:
        return allGenomes[name]
    if name == "hg19":
        genome_db = genome.Genome("../data/hg19")
    elif name == "hg18":
        genome_db = genome.Genome("../data/hg18")
    elif name == "mm9": 
        genome_db = genome.Genome("../data/mm9")    
    elif name == "mm10": 
        genome_db = genome.Genome("../data/mm10")                    

    #You can also use genomes with only numbered and X chromosomes
    #genome_db = genome.Genome("../data/hg19", readChrms=["#","X"])
    elif name == "dm3": 
        #Drosophila (example of specifying exact chromosomal order)
        genome_db = genome. Genome("../data/dm3", 
                           readChrms=["2L","2R","3L","3R","4","X","2LHet","2RHet","3LHet","3RHet","XHet","YHet","U","Uextra","M"]
                           ,forceOrder=True)
    elif name == "cb10":
        genome_db = genome.Genome("../data/cb10", readChrms=["I","II","III","IV","V","X","M"], forceOrder=True)        
    else: 
        raise ValueError("Genome {0} not defined. Edit defineGenome.py and define it".format(name))
    allGenomes[name] = genome_db
    return genome_db

