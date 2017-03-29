from singleShared import combinedDf, singleDf, df, averageLoopsWithControl, averageDomains, doSaddle, eigenDict, doScaling, loopDict, genomDict, getCoverage
import pickle
from mirnylib.h5dict import h5dict
import numpy as np
import os
from multiprocessing import Pool
from mirnylib.systemutils import  setExceptionHook
setExceptionHook()

loopDomainResolution = 10000
loopWindowsSize = 8
tadSize = 30

loopDomainResolution = 5000
loopWindowsSize = 12
tadSize = 60




def doOne(ind):
    print(ind)
    cur = df.loc[ind]
    gen = cur["genome"]
    fnameref = os.path.join(gen, cur["filenames"] + "_refined.frag")
    fname10 = os.path.join(gen, cur["filenames"] + ".{0}.cool".format(loopDomainResolution))
    fname1k = os.path.join(gen, cur["filenames"] + ".200000.cool")

    coverages = getCoverage(fname10)
    pickle.dump(coverages, open("coverages/{0}".format(ind),'wb'))

    loops = averageLoopsWithControl(loopDict[gen][0], fname10, pad=loopWindowsSize)
    pickle.dump(loops, open("loops/{0}".format(ind),'wb'))

    tads = averageDomains(fname10, M=tadSize)
    pickle.dump(tads, open("tads/{0}".format(ind),'wb'))

    #scal = doScaling(fnameref)
    #pickle.dump(scal, open("scalingsNew/{0}".format(ind),'wb'))
    comp = doSaddle(fname1k, eigenDict[gen], gen)
    pickle.dump(comp, open("saddles/{0}".format(ind),'wb'))

doOne(df.index.values[0])

mypool = Pool(24)

mypool.map(doOne, df.index.values[1:])
#scal = doScaling("mm9/F123_ES_2014_Ren-all-HindIII_refined.frag")
#pickle.dump(scal, open("scalingsNew/F123_ES",'wb'))











