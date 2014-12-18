from hiclib.fragmentHiC import HiCdataset
from mirnylib.h5dict import h5dict
import os
import sys
import numpy as np
from mirnylib.systemutils import setExceptionHook

if os.path.exists("test-1M.hm"):
    os.remove("test-1M.hm")
if os.path.exists("test-1M-byChr.hm"):
    os.remove("test-1M-byChr.hm")

workingGenome = "hg19"
genomeFolder = "../../fasta/hg19"
if not os.path.exists(genomeFolder):

    genomeFolder = sys.argv[1]
    if not os.path.exists(genomeFolder):
        raise StandardError("Please provide hg18 Genome folder in the code or as a first argument")


def source(ID):
    return os.path.join("%s-%s.hdf5" % (ID, workingGenome))  # determines path to the parsed file by ID


def refine_paper(filename, create=True):
    """filename[0] is a list of filenames of incoming files
    filename[1] is a folder for outgoing file"""
    if create == True:
        for onename in filename[0]:
            #Parsing individual files
            if not os.path.exists(onename):
                raise StandardError("path not found: %s" % onename)
            TR = HiCdataset("bla", genome=genomeFolder, enzymeName="HindIII",maximumMoleculeLength=500, inMemory=True)
            print "\nTesting loading new data without rsite information    "
            TR.parseInputData(dictLike=onename,
                              enzymeToFillRsites="HindIII")
            #assert len(TR.DS) == 856143

            #assert len(TR.ufragments) == 634572
            TR.save(onename + "_parsed.frag")

        #Merging files alltogether, applying filters
        TR = HiCdataset(filename[1] + "_merged.frag",enzymeName = "HindIII",
                        genome=genomeFolder, mode="w")
        TR.merge([i + "_parsed.frag" for i in filename[0]])

        TR = HiCdataset("refined", genome=genomeFolder,enzymeName = "HindIII",
                        mode="w", inMemory=True)

        print "\nTesting chunking during all tests"
        TR.chunksize = 30000
        #because we do many operations, we disable autoFlush here
        TR.load(filename[1] + "_merged.frag")

        print "\nTesting Rsite filter"
        TR.filterRsiteStart(offset=5)

        #assert len(TR.DS) == 832110

        print "\nTesting duplicate filter"
        TR.filterDuplicates()

        #assert len(TR.DS) == 830275

        print "\nTesting small/large and extreme fragment filter"
        TR.filterLarge()

        #assert len(TR.DS) == 825442
        TR.filterExtreme(cutH=0.005, cutL=0)
        TR.writeFilteringStats()

        #assert len(TR.DS) == 803845


    #-------------------------------------------
    TR.printMetadata(saveTo="metadata")
    import cPickle

    stop = False
    mdata = cPickle.load(open("sampleMetadata"))
    for i in sorted(mdata.keys()):
        if TR.metadata[i] != mdata[i]:
            print "Key {0} is not consistent: should be {1}, is {2}".format(i, mdata[i], TR.metadata[i])
            stop = True
    if stop == True:
        raise ValueError("Inconsistent metadata: see above")


    print "Testing allxall and by-chromosome heatmap counting diagonal twice"

    TR.printStats()
    print "----> saving allxall heatmap"
    TR.saveHeatmap(filename[1] + "-1M.hm", 1000000,
                   countDiagonalReads="twice")
    a = h5dict(filename[1] + "-1M.hm")
    st, end = TR.genome.chrmStartsBinCont[1], TR.genome.chrmEndsBinCont[1]
    st2, end2 = TR.genome.chrmStartsBinCont[2], TR.genome.chrmEndsBinCont[2]
    chrom1 = a["heatmap"][st:end, st:end]
    chrom12 = a["heatmap"][st:end, st2:end2]
    setExceptionHook()
    print "----> saving by chromosome heatmap"
    TR.saveByChromosomeHeatmap(
        filename[1] + "-1M.hm", resolution=1000000, includeTrans=True,
        countDiagonalReads="twice")

    b = h5dict(filename[1] + "-1M.hm")["1 1"]
    bb = h5dict(filename[1] + "-1M.hm")["1 2"]
    assert (b - chrom1).sum() == 0
    print "Cis heatmap consistent"
    assert (bb - chrom12).sum() == 0
    print 'Trans heatmap consistent'
    print  a["heatmap"][::10, ::10].sum()
    #assert  a["heatmap"][::10, ::10].sum() == 21800
    print "Heatmap sum correct\n"

    #---------------------------------
    print "Testing allxall and by-chromosome heatmap counting diagonal once"

    TR.saveHeatmap(filename[1] + "-1M.hm", 1000000,
                   countDiagonalReads="once")
    Ta = h5dict(filename[1] + "-1M.hm")
    st, end = TR.genome.chrmStartsBinCont[1], TR.genome.chrmEndsBinCont[1]
    st2, end2 = TR.genome.chrmStartsBinCont[2], TR.genome.chrmEndsBinCont[2]
    chrom1 = Ta["heatmap"][st:end, st:end]
    chrom12 = Ta["heatmap"][st:end, st2:end2]
    setExceptionHook()
    print "----> saving by chromosome heatmap"
    TR.saveByChromosomeHeatmap(
        filename[1] + "-1M-byChr.hm", resolution=1000000, includeTrans=True,
        countDiagonalReads="once")

    Tb = h5dict(filename[1] + "-1M-byChr.hm")["1 1"]
    Tbb = h5dict(filename[1] + "-1M-byChr.hm")["1 2"]
    assert ((Tb - chrom1) == 0).all()
    assert ((Tbb - chrom12) == 0).all()
    assert ((Tb + np.diag(np.diag(Tb))) == b).all()
    print "Diagonal counting methods are consistent\n"

    #------------------------------
    print "Testing updateGenome method"
    from mirnylib.genome import Genome
    removeChromIDs = np.array([0, 1, 1, 1, 1] + [0] * 17 + [1] + [0])
    #print ((removeChromIDs[TR.chrms1] == 1) + (removeChromIDs[TR.chrms2] == 1) ).sum()
    t = ((removeChromIDs[TR.chrms1] == 1) * (removeChromIDs[TR.chrms2] == 1)).sum() + ((removeChromIDs[TR.chrms1] == 1) * (TR.chrms2 == -1)).sum()
    newGenome = Genome(genomePath=genomeFolder, readChrms=["2",
                                                           "3", "4", "5", "X"])
    TR.updateGenome(newGenome)
    assert  TR.N == t

    a = h5dict(filename[1] + "-1M.hm")["heatmap"]


map(refine_paper,
    [
    [(source("test"),
      ), "test", "HindIII"]])

#os.remove("test_breaks.frag")
os.remove("test-hg19.hdf5_parsed.frag")
os.remove("test_merged.frag")
#os.remove("test-1M.hm")
#os.remove("test_refined.frag")
print "Test finished successfully!"
