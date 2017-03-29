import matplotlib
import mirnylib.plotting
import numpy as np
from scipy.stats import gaussian_kde
from mirnylib.plotting import nicePlot
import matplotlib.pyplot as plt
from singleShared import typeDict, df, secondTypeDict, singleDf
import pickle
import functools
from mirnylib.numutils import observedOverExpected, rank, ultracorrect, coarsegrain
from scipy.stats import pearsonr
import matplotlib.cm as cm

######## DISPLAY via SSH
#import os
#os.environ['DISPLAY'] = 'localhost:10.0'
#######################

td = typeDict
backPeers = {}
backMain = {}
for i,j in typeDict.items():
    for k in j:
        backPeers[k] = j
        backMain[k] = i

def norm(scal):
    x,y = scal
    mask = np.isfinite(y)
    y = y[mask]
    x = x[mask]
    norm = y[3]
    y = y / norm
    return x,y


def findTads(tads, target):
    means = []
    keys = sorted(tads.keys())
    for i in keys:
        onetads = np.concatenate(tads[i])
        mean = (onetads[:,1] - onetads[:,0]).mean()
        means.append(mean)
    print(means)
    ind = np.argmin(np.abs(np.array(means) - target))
    onetads = np.concatenate(tads[keys[ind]])
    mean = (onetads[:,1] - onetads[:,0]).mean()
    print(mean, target)
    return tads[keys[ind]]


def smartcg(inAr, by):
    cg1 = coarsegrain(inAr, by)
    newAr = np.zeros_like(inAr)
    for i in range(by):
        for j in range(by):
            newAr[i::by, j::by] = cg1
    return newAr / by**2


def nicePlot(ax="gca", fs=8, show=True):
    """
    replaces obsolete "niceShow" command, packs it with new features
    """
    if ax == "gca":
        ax = plt.gca()
    matplotlib.rcParams.update({'font.size': fs})

    legend = ax.legend(loc=0, prop={"size": fs })
    if legend is not None:
        legend.draw_frame(False)
    mirnylib.plotting.removeAxes(shift=0, ax=ax)

    try:
        plt.tight_layout(pad=0.3)
    except:
        pass
    if show:
        plt.show()



def load(pref, ind):
    return pickle.load(open("{0}/{1}".format(pref, ind), 'rb'))

def mysum(x):
    return functools.reduce(lambda x,y:x+y, x)

def statLoops(ind, cg):
    loops, controls = load("loops", ind)
    l1 = 5. * mysum(loops[::2])
    c1 = 1. * mysum(controls[::2])
    l2 = 5. * mysum(loops[1::2])
    c2 = 1. * mysum(controls[1::2])
    if (np.sum(c1 == 0) > 0) or (np.sum(c2 == 0) > 0):
        cg = 2
        print("Coarsegrained", ind)
    if cg != 1:
        l1 = smartcg(l1, cg)
        l2 = smartcg(l2, cg)
        c1 = smartcg(c1, cg)
        c2 = smartcg(c2, cg)
    r1 = l1 / c1
    r2 = l2 / c2
    r3 = (l1 + l2) / (c1 + c2)
    return r1, r2, r3

def statDomains(ind):
    tads = load("tads", ind)
    t1 = mysum(tads[::2])
    t2 = mysum(tads[1::2])
    t3 = observedOverExpected(t1 + t2)
    t1 = observedOverExpected(t1)
    t2 = observedOverExpected(t2)

    return t1, t2, t3

def statSaddles(ind):
    sads = load("saddles", ind)
    s1 = mysum(sads[::2])
    s2 = mysum(sads[1::2])
    s3 = s1 + s2
    def process(sad):
        sad = sad / np.mean(sad)
        sad = ultracorrect(sad)
        sad = np.log(sad)
        return sad

    return [process(i) for i in [s1, s2, s3]]

def statScalings(ind):
    s1, s2, s3 = load("scalingsNew", ind)
    return [norm(i) for i in [s1, s2, s3]]


colors = list(cm.rainbow(np.linspace(0, 1, 5)))


df = df.sort("cis", ascending = False)


def mainFigTadCompEtc(cells):


    matplotlib.rcParams.update({'font.size': 5})
    plt.figure(figsize=(7,5))
    b = 3
    a = len(cells)
    titles = [ "Even chromosomes", "Odd chromosomes"]
    heatmaps = {}
    for cellid, cell in enumerate(cells):
        plt.figure(1)
        sad = statSaddles(cell)


        plt.subplot(a,b, 3*cellid+1)
        #plt.title(cell)

        plt.imshow(sad[2], vmin = -0.5, vmax = 0.5, interpolation="none", cmap = "coolwarm")
        heatmaps["Saddle_{0}".format(cell)] = sad[2]
        plt.xticks([0,4],["active", "inactive"])
        plt.yticks([0,4],["active", "inactive"], rotation="vertical")
        plt.colorbar(label="Log enrichment")




        plt.subplot(a,b, 3*cellid+3)

        tad = statDomains(cell)[2]
        ar = np.arange(len(tad))
        mask = 1 / (1 + np.abs(ar[:,None] - ar[None,:])) ** 0.25
        tad = tad * mask

        heatmaps["TAD_{0}".format(cell)] = tad
        plt.imshow(tad, cmap = "fall", interpolation="none", vmin = 0.3, vmax=0.9)
        plt.colorbar(label="Effective contact probability")



        plt.subplot(a,b,3*cellid+2 )
        plt.title(cell)


        loop = statLoops(cell, 1)[2]
        heatmaps["Loop_{0}".format(cell)] = loop
        loop = np.log(loop)

        plt.imshow(loop, vmin = -0.8, vmax = 0.8, cmap="coolwarm" ,interpolation="none")
        plt.colorbar(label="Log enrichment")
        plt.xticks([0,8,15], ["-80kb","Upstream loop base", "+70kb"])
        plt.yticks([0,8,15], ["-80kb","Downstream loop base", "+70kb"], rotation="vertical")

    plt.show()
    exit()

mydf = df[(df["type"] == "SN")]

#run this ---------------------!!!!!!!!!!--------------
print(mydf.index.values[:3])
#mainFigTadCompEtc(["oocyte"] + list(mydf.index.values[:3]))

#run this ---------------------!!!!!!!!!!--------------
#mainFigTadCompEtc(["pronucleus-female", "pronucleus-male", "SN-Hoechst", "NSN-Hoechst"])

#exit()



def doOne(i):
    f1 = plt.figure(figsize=(3,10))
    a = 4
    b = 1
    plt.subplot(a,b,1)
    plt.title("saddle plot")
    sad = statSaddles(i)[2]
    plt.imshow(sad, vmin = -0.5, vmax = 0.5, cmap="coolwarm", interpolation="none")
    plt.xticks([0,4],["active", "inactive"])
    plt.yticks([0,4],["active", "inactive"], rotation="vertical")
    plt.colorbar(label="Log enrichment")

    plt.subplot(a,b,2)
    plt.title("Average TAD")

    tad = statDomains(i)[2]

    ar = np.arange(len(tad))
    mask = 1 / (1 + np.abs(ar[:,None] - ar[None,:])) ** 0.25
    tad = tad * mask


    plt.imshow(tad, cmap = "fall", interpolation="none", vmin = 0.3, vmax=0.9)
    plt.colorbar(label="Effective contact probability")

    plt.subplot(a,b,3)
    plt.title("Average loop")
    cg = 1

    loop = statLoops(i, cg)[2]
    loop = np.log(loop)
    #print(loop)
    plt.imshow(loop, vmin = -0.8, vmax = 0.8, interpolation="none", cmap = "coolwarm")
    plt.colorbar(label="Log enrichment")
    plt.xticks([0,8,15], ["-80kb","Upstream loop base", "+70kb"])
    plt.yticks([0,8,15], ["-80kb","Downstream loop base", "+70kb"], rotation="vertical")

    plt.subplot(a,b,4)
    if i in backPeers:
        peers = backPeers[i]
        first = True
        for peer in peers:
            scal = statScalings(peer)[2]
            if first:
                plt.plot(*scal, color = "gray", linewidth=0.5, label="All {0}".format(backMain[i]))
                first = False
            else:
                plt.plot(*scal, color = "gray", linewidth=0.5)
    if i in backMain:
        main = backMain[i]
        scal = statScalings(main)[2]
        plt.plot(*scal, label="average {0}".format(backMain[i]), color = "black")
    scal = statScalings(i)[2]

    plt.plot(*scal, label=str(i), color = "blue")
    try:
        plt.xscale("log")
        plt.yscale("log")
    except:
        pass
    plt.xlabel("Distance (bp)")
    plt.ylabel("Relative contact probability")
    nicePlot(show=False, fs=8)
    if i in backMain:
        extraname = backMain[i]
    else:
        extraname=""
    try:
        plt.savefig("total/{0}-{1}.pdf".format(i, extraname))
    except:
        plt.xscale("linear")
        plt.yscale("linear")
        plt.savefig("total/{0}-{1}.pdf".format(i, extraname))
    del f1


from multiprocessing import Pool
print(len(df))
mypool = Pool(20)
list(mypool.map(doOne, df.index.values))
[doOne(i) for i in df.index.values]
