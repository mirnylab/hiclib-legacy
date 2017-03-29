import itertools
from qweutils.matrix import hierarchicalReorder
from hiclib.binnedData import binnedData
import os
import random
import pickle
from mirnylib.h5dict import h5dict
from mirnylib.numutils import coarsegrain, observedOverExpected
import numpy as np
import glob
import hiclib.hicManager as manager
import joblib
import matplotlib.pyplot as plt
import mirnylib
import pandas as pd
from mirnylib.genome import Genome
from mirnylib.plotting import nicePlot
from hiclib.hicShared import getResolution
from singleShared import doSaddle, doEigenvector, doSaddleError
mirnylib.plotting.setLongColorCycle()
import functools
mem = joblib.Memory(".")

df = pd.read_csv("statistics.csv", index_col=0)

df = df[np.array(["combined" in str(i) for i in df["filenames"].values])]

def mysum(x):
    return functools.reduce(lambda x,y:x+y, x)



HiCHuman = glob.glob("/home/magus/HiC2011/Erez2014/hg19/*inSitu*all*MboI*1000k.hm") + ["/home/magus/HiC2011/BingRen2012Human/hg19/HES-all-HindIII-1000k.hm" , "GC"]

HiCMouse = ["BingRen2013Mouse/mm9/F123_ES_2014_Ren-all-HindIII-mm9","Erez2014/mm9/CH12-LX_inSitu-all-MboI-mm9",
 "hadjurCohesin2012/mm9/NSC_Grande-mm9", "mouseZhang2012/mm9/Pro_B_ATM-all-HindIII-mm9",
 "BingRen2012Mouse/mm9/MouseES-all-HindIII-mm9", "hadjurCohesin2012/mm9/AST_Grande-mm9",
"BingRen2012Mouse/mm9/MouseES-NcoI-R1-NcoI-mm9", "BingRen2012Mouse/mm9/MouseESCort-all-HindIII-mm9",
            "Sperm2015/mm9/sperm-all-HindIII-mm9","Sperm2015/mm9/fibroblastControl-all-HindIII-mm9",
            #"DrosophilaSingleCell2015/alternativeFiltering/mm9/oocyte_combined-mm9"
            ]
HiCMouse = ["/home/magus/HiC2011/" + i[:-4] + "-1000k.hm" for i in HiCMouse] + ["GC"]
hicBig = {"hg19":HiCHuman, "mm9":HiCMouse}

for i in HiCHuman + HiCMouse:
    if not os.path.exists(i):
        print(i)

def saddleStrength(saddle, permutted):
    def strength(i):
        # comment return below to make kind of weird but realistic 50/50 compartment split 
        # rather than default 20% compartments
        # ------------!!!!!!!!!!!!!!!!!!!!!!-------------------
        return np.log(i[0,0] * i[-1,-1] / i[0,-1]**2)
        t1 = i[:3,:3].sum() + i[:2, :2].sum()
        t2 = i[2:, 2:].sum() + i[3:, 3:].sum()
        t3 = i[:3, 2:].sum() + i[:2, 3:].sum()
        return np.log(t1 * t2 / t3**2)
    s1 = strength(saddle)
    sp = [strength(i) for i in permutted]
    return s1, np.std(sp)



for gen in ["mm9"]:
    if gen == "mm9":
        cur = ["mm9/pronuc-female_combined-1000k.hm", "mm9/pronuc-male_combined-1000k.hm",
               "mm9/pronuc-w-o-inh-female_combined-1000k.hm", "mm9/pronuc-w-o-inh-male_combined-1000k.hm",
               #"mm9/pronuc-w-o-inh-both_combined-1000k.hm", "mm9/pronuc-w-o-inh-undef sex_combined-1000k.hm",
               "mm9/NSN-Hoechst_combined-1000k.hm", "mm9/SN-Hoechst_combined-1000k.hm",
               "../../BingRen2013Mouse/mm9/F123_ES_2014_Ren-all-HindIII-1000k.hm",
               "mm9/sperm-all-HindIII-1000k.hm"
               ]
        #Uncomment below to make extra plots ------------------!!!!!!!!!!!!!!!!!!!--------------
        cur = ["mm9/SN_combined-1000k.hm", "mm9/NSN_combined-1000k.hm",
               "mm9/SN-Hoechst_combined-1000k.hm", "mm9/NSN-Hoechst_combined-1000k.hm",]

        #cur = hicBig[gen]

    eig  = doEigenvector(hicBig[gen][0], gen)

    eig = "GC"
    names = []
    values = []
    errors = []
    plt.figure(figsize=(7,5))
    for i in cur:

        plt.title("Compartment strength in pronuclei" )
        if "bingren" in i.lower():
            sad, permutted = doSaddleError(i, eig, gen, correct=True)
        else:
            sad, permutted = doSaddleError(i, eig, gen)
        stren, err = saddleStrength(sad, permutted)
        values.append(stren)
        errors.append(err)
        names.append(i[4:].split("combined")[0])
    writer = pd.ExcelWriter('output1.xlsx')
    pd.DataFrame([values, errors], index = ["value","standard error"], columns=names).to_excel(writer,"3d")
    writer.close()

    args = np.argsort(values)
    names = [names[i] for i in args]
    values = [values[i] for i in args]
    errors = [errors[i] for i in args]
    index = np.arange(len(names))
    plt.bar(index, values, 0.8,
            yerr = errors,
            error_kw= {'ecolor': '0.3'})
    print(names)
    print(values)
    print(errors)
    plt.xticks(index, names, rotation="vertical")
    plt.show()
    #nicePlot(fs=10)








for gen in ["hg19","mm9"]:
    if gen == "mm9":
        #cur = ["mm9/oocyte_combined-1000k.hm", "mm9/pronucleus_combined-1000k.hm", "../../Sperm2015/mm9/sperm-all-HindIII-1000k.hm"]
        cur = ["mm9/oocyte_combined-1000k.hm",
            "mm9/SN_combined-1000k.hm", "mm9/NSN_combined-1000k.hm",
               "mm9/SN-Hoechst_combined-1000k.hm", "mm9/NSN-Hoechst_combined-1000k.hm",]
        print("yay")
    if gen == "hg19":
        cur = ["hg19/K562_combined-1000k.hm"]
        #cur = hicBig[gen]
    eigs  = [doEigenvector(i, gen) for i in hicBig[gen]]
    for i in cur:
        plt.figure(figsize=(6,4))
        plt.title("Compartment strength in " + i)
        names = []
        values = []
        errors = []
        for j,eig in enumerate(eigs):
            if "GC" in hicBig[gen][j]:
                namej = "GC"
            else:
                namej = hicBig[gen][j].split(gen)[1][1:20]
            sad, permutted = doSaddleError(i, eig, gen)
            stren, err = saddleStrength(sad, permutted)
            values.append(stren)
            errors.append(err)
            names.append(namej)
        args = np.argsort(values)
        names = [names[i] for i in args]
        values = [values[i] for i in args]
        errors = [errors[i] for i in args]

        writer = pd.ExcelWriter('figED3.xlsx')
        pd.DataFrame([values, errors], index=["value", "standard error"], columns=names).to_excel(writer, "ED_3a")
        writer.close()

        index = np.arange(len(names))

        #pd.DataFrame([values, errors], index=["value", "standard error"], columns=names).to_excel(writer, "supFig3a")

        plt.bar(index, values, 0.8,
                yerr = errors,
                error_kw= {'ecolor': '0.3'})
        plt.xticks(index, names, rotation="vertical")
        nicePlot(fs=10)





#for gen in ["hg19", "mm9"]:
for gen in []:
    print(df)
    mydf = df[df["genome"] == gen]
    fnames = [os.path.join(gen, i + "-1000k.hm") for i in mydf["filenames"].values] + hicBig[gen]
    eigs  = [doEigenvector(i, gen) for i in hicBig[gen]]

    for i, fname in enumerate(fnames):
        for j, eig in enumerate(eigs):
            namei = os.path.split(fname)[-1][:20]
            namej = hicBig[gen][j].split(gen)[1][1:11]
            #plt.subplot(len(fnames), len(eigs), i * len(eigs) + j + 1)
            sad = doSaddle(fname, eig, gen)
            sad = mysum(sad)
            plt.subplot(len(eigs), len(fnames), j * len(fnames) + i + 1)
            plt.imshow(np.log(sad), interpolation = "none", vmin = -0.5, vmax = 0.5)
            #plt.colorbar()
            #plt.title(fname[:5])
            plt.xlabel(namei, fontsize=5)
            plt.ylabel(namej, fontsize=5)
            #nicePlot(show=False, fs=5)

            print(i,j)


    nicePlot(fs=5)




for gen in ["hg19", "mm9"]:
    print(df)
    mydf = df[df["genome"] == gen]
    fnames = hicBig[gen] + sorted([os.path.join(gen, i + "-1000k.hm") for i in mydf["filenames"].values])
    fnames = [i for i in fnames if "GC" not in i]
    eigs  = [doEigenvector(i, gen) for i in hicBig[gen]]
    allsads = []
    for i, fname in enumerate(fnames):
        sads = []
        for j, eig in enumerate(eigs):
            #plt.subplot(len(fnames), len(eigs), i * len(eigs) + j + 1)
            sad = doSaddle(fname, eig, gen)
            sad = mysum(sad)
            sads.append(np.log((sad[0,0] * sad[-1,-1]) / (sad[0,-1] ** 2)))
            #nicePlot(show=False, fs=5)
        allsads.append(sads)
    #allsads, args1, args2 = hierarchicalReorder(allsads, useCorr=False)
    args1 = list(range(len(fnames)))
    args2 = list(range(len(hicBig[gen])))
    namesi = [os.path.split(fnames[j])[-1][:20] for j in args1]
    hicBig[gen][-1] = "{0}/GC".format(gen)
    namesj = [hicBig[gen][j].split(gen)[1][1:20] for j in args2]
    plt.imshow(np.transpose(allsads), interpolation="none", vmin = 0, vmax = 2)
    plt.colorbar()
    plt.xticks(range(len(namesi)), namesi, rotation="vertical")
    plt.yticks(range(len(namesj)), namesj)

    nicePlot(fs=7)

