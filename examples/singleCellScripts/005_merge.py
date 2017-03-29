from hiclib.fragmentHiC import HiCdataset
import pickle
import pandas as pd
import os
from multiprocessing import Pool
df = pd.read_csv("statistics.csv", index_col=0)

def doOne(inData):
    coolResolutions = [10000000, 5000000, 2000000, 1000000, 500000, 200000, 100000, 40000, 20000, 10000, 5000, 2000,
                       1000]
    i, j = inData
    if i == "":
        return
    # if i not in ["total", "pronuc", "K562"]:
    #    continue
    print(i)
    gens = j["genome"]
    if len(gens) == 0:
        print("Genome not found")
        return
    genome = gens.values[0]
    out_file = "{1}/{0}_combined".format(i, genome)
    # if os.path.exists(out_file + "-10k_HighRes.byChr"):
    #    continue
    mygen = "/home/magus/HiC2011/data/{0}".format(genome)
    filenames = ["{1}/{0}_refined.frag".format(s, genome) for s in j["filenames"].values]

    # assert False not in list(map(os.path.exists, filenames))
    filenames = [i for i in filenames if os.path.exists(i)]
    if len(filenames) == 0:
        print("No filenames found!")
        return
    TR = HiCdataset("bla", mygen, "DpnII", inMemory=True, tmpFolder="/tmp", dictToStoreIDs="dict")
    TR.merge(filenames)
    TR.setSimpleHighResHeatmap()
    TR.writeFilteringStats()
    TR.printMetadata(saveTo="statistics/{1}/{0}_combined.frag".format(i, genome))
    pickle.dump(TR.metadata, open("statistics/{1}/{0}_combined.pkl".format(i, genome), 'wb'))
    TR.save("{1}/{0}_combined_refined.frag".format(i, genome))
    for res in coolResolutions:
        TR.saveCooler(out_file + ".{0}.cool".format(res), res)


mydf = df[df["type"] != "combined"]

mypool = Pool(16)
for field in ["type", "type2", "type3"]:

    splited = pd.groupby(mydf, mydf[field])
    print(set(mydf[field]))
    tomap = [(i,j) for i,j in splited]
    list(map(doOne, tomap))


        


