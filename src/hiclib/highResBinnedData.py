# (c) 2012 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Maksim Imakaev (imakaev@mit.edu)

"""
This is a basic toolbox to perform high-resolution analysis of a Hi-C data.
By high resolution we meen that there are more than 10000-20000 bins per genome
- otherwise, binnedData is easily capable of the analysis.

It can perform basic filtering (poor coverage regions), and iterative correction at any resolution.

The main advantage of this class is that it supports both in memory and HDD storage,
and for HDD storage it supports both sparse and dense matrix logic.
In fact, sparse logic with HDF5-based storage in memory is already good (default settings).

.. note::
    This class loads data, saved by saveByChromosomeHeatmap method of fragmentHiC.
    It is not designed to load heatmaps saved by "saveHeatmap" class, because those
    can fit in memory, and should be analyzed by a more powerful binnedData class.





Class structure
---------------

Class defaultMatrix implements a wrapper around a 2D numpy.array matrix. This
class defines methods that will be applied to a Hi-C map between two
chromosomes. This class should serve as a template and backup for other Hi-C map
storage classes which will be inhereted from this class.

Class h5dictMatrix implements a subclass of a defaultMatrix, where the matrix is
in fact stored on the HDD. It leaves all other methods intact, as  all other
methods use only self.setData() and self.getData(). Therefore overriding getter
and setter is enough to move matrix storage from memory to an h5dict.
Use defaultMatrix when the speed is critical, if you're working with cis data only.

Class h5dictSparseMatrix implements a subclass of defaultMatrix, where all the
methods were overwrridden to operate in a sparse logic. It stores a Hi-C map as
three arrays: X coordinate, Y coordinate and Value. Therefore, for sparse
matrices all overriden methods will be much faster than default, as they will
bypass calling self.getData and creating a huge matrix in memory. It is
suggested to use h5dictMatrix for cis maps and h5dictSparseMatrix for trans
(between-chromosomal) maps.

Class HiResHiC allows to load a by-chromosome Hi-C map, with either cis only,
or cis and trans matrices. It then allows to do an iterative correction of the
map, and allows for some other filtering. Unlike binnedData, it does not support
multiple datasets at once, or complicated filtering, or PCA. It will support by-
chromosome domain finder from cis reads at some point as we converge on the
ideal algorithm to do this.

On my machine with 32GB RAM it was successfully used to perform IC of a human
Hi-C at 10kb resolution. It took it 20 minutes to load the data, 4 mins to
remove poor bins, and couple hours to perform IC, about 10 minutes per pass. I
also note that the data was in fact stored in memory for doing that, and it
never used more than 16GB of RAM... in fact, creating this dataset used more,
but this will be optimized later.
"""


from mirnylib.genome import Genome
import numpy as np
import warnings
from mirnylib.h5dict import h5dict
from mirnylib import numutils
from mirnylib.systemutils import setExceptionHook
from scipy.stats.stats import spearmanr
from hiclib import hicShared

setExceptionHook()


class defaultMatrix(object):
    """
    This is a template object which stores matrix in memory.
    Alternatively, matrix can be stored in an h5dict, either normally
    or in a sparse mode.

    All the methods should be first implemented here using getData() and setData()
    Then they shold be translated to sparse subclasses of this class.
    """
    def __init__(self, data=None, dictToSave=None, key=""):
        """
        Initializes the object that stores the Hi-C matrix.

        Parameters
        ----------

        data : 2D matrix in any format
            Hi-C matrix between two chromosomes
        dictToSave : dict, h5dict or any other dict-like structure
            Dict to store actual data. Should be h5dict for high resolution analysis.
            Is not needed for defaultMatrix, will be used by it's subclasses only.
        key : str or anything
            A key, unique for each pair of chromosomes, used to identify dataset
            in the dictToSave. Is provided by HiResHiC.
        """
        self._h5dict = dictToSave
        self._key = repr(key)

        if data is not None:
            self.setData(data)

    def getData(self):
        "Returns a matrix in a dense format (NxM array)"
        return self.data.copy()

    def setData(self, data):
        """Accepts the matrix to store. Here, just puts it in RAM.
        For further subclasses this will need to convert the
        matrix to the sparse form.
        """
        data = np.asarray(data, dtype=np.float64)
        assert len(data.shape) == 2
        self.data = data

    def getSumX(self):
        "returns sum of all values along the first (0) axis, i.e. sum of all rows"
        return np.sum(self.getData(), axis=0)

    def getSumY(self):
        "returns sum of all values along the second (1) axis, i.e. sum of all columns"
        return np.sum(self.getData(), axis=1)

    def getSums(self):
        "returns getSumX and getSumY without calling getData() twice"
        data = self.getData()
        sumX = np.sum(data, axis=0)
        sumY = np.sum(data, axis=1)
        return (sumX, sumY)

    def clearRows(self, rows):
        """Sets to 0 certain rows and columns

        Parameters
        ----------
        rows : tuple of two arrays (rows, columns)
            Two arrays which hold indices of rows and colums to be removed
        """
        rowsX, rowsY = rows
        data = self.getData()
        data[rowsX, :] = 0
        data[:, rowsY] = 0
        self.setData(data)

    def divideByVectorX(self, vectorX):
        """
        Divides each row by correspoinding value from vectorX
        """
        vectorX[vectorX == 0] = 1
        data = self.getData()
        data /= vectorX[None, :]
        self.setData(data)

    def divideByVectorY(self, vectorY):
        """
        Divides each column by correspoinding value from vectorY
        """

        vectorY[vectorY == 0] = 1
        data = self.getData()
        data /= vectorY[:, None]
        self.setData(data)

    def divideByVectors(self, vectors):
        """
        Divides each row and column by correspoinding
         value from vectors[0] and vectors[1]

        Does it without calling getData twice!
        """

        vecX = vectors[0]
        vecY = vectors[1]
        vecX[vecX == 0] = 1
        vecY[vecY == 0] = 1
        data = self.getData()
        assert data.shape[1] == len(vecX)
        assert data.shape[0] == len(vecY)
        data /= vecX[None, :]
        data /= vecY[:, None]
        self.setData(data)

    @property
    def shape(self):
        """Returns shape of the data.
        Should be overridden in subclasses
        not to load the data every time! """
        return self.getData().shape


class h5dictMatrix(defaultMatrix):
    """
    Changes the storage from memory to h5dict, keeping the matrix
    in the dense (regular) format.
    Just overrides getData, setData and shape to use h5dict.
    """
    def getData(self):
        data = self._h5dict[self._key]
        return data

    def setData(self, data):
        data = np.asarray(data, dtype=np.float64)
        self.savedShape = data.shape
        self._h5dict[self._key] = data

    @property
    def shape(self):
        return self.savedShape


class h5dictSparseMatrix(defaultMatrix):
    """
    Changes the storage from memory to h5dict,
    and changes matrix to sparse.
    All methods from defaultMatrix are overridden here!
    """

    def __init__(self, data=None, dictToSave=None, key=""):
        self._h5dict = dictToSave
        self._key = repr(key)

        self._keyx = self._key + "x"
        self._keyy = self._key + "y"
        self._keyv = self._key + "v"
        if data is not None:
            self.setData(data)

    @property
    def _X(self):
        return self._h5dict[self._keyx]

    @_X.setter
    def _X(self, data):
        self._h5dict[self._keyx] = data

    @property
    def _Y(self):
        return self._h5dict[self._keyy]

    @_Y.setter
    def _Y(self, data):
        self._h5dict[self._keyy] = data

    @property
    def _V(self):
        return self._h5dict[self._keyv]

    @_V.setter
    def _V(self, data):
        self._h5dict[self._keyv] = data

    @property
    def shape(self):
        return self.savedShape

    def getData(self):
        data = np.zeros(self.savedShape)
        data[self._X, self._Y] = self._V
        return data

    def setData(self, data):
        x, y = np.nonzero(data)
        self.savedShape = data.shape
        values = data[x, y]
        values = np.asarray(values, np.float64)
        self._X = x
        self._Y = y
        self._V = values

    def getSumX(self):
        return np.bincount(self._Y, weights=self._V, minlength=self.shape[1])

    def getSumY(self):
        return np.bincount(self._X, weights=self._V, minlength=self.shape[0])

    def getSums(self):
        X, Y, V = self._X, self._Y, self._V
        s1 = np.bincount(Y, weights=V, minlength=self.shape[1])
        s2 = np.bincount(X, weights=V, minlength=self.shape[0])
        return (s1, s2)

    def divideByVectorX(self, vecX):
        vecX[vecX == 0] = 1
        self._V = self._V / vecX[self._Y]

    def divideByVectorY(self, vecY):
        vecY[vecY == 0] = 1
        self._V = self._V / self.vecX[self._X]

    def divideByVectors(self, vecs):
        vecX = vecs[0]
        vecY = vecs[1]
        V = self._V
        V /= vecX[self._Y]
        V /= vecY[self._X]
        self._V = V

    def clearRows(self, rows):
        rowsX, rowsY = rows
        indexX = np.ones(self.shape[0], bool)
        indexY = np.ones(self.shape[1], bool)
        indexX[rowsX] = False
        indexY[rowsY] = False
        X = self._X
        Y = self._Y
        mask = indexX[X]
        mask *= indexY[Y]
        self._X = X[mask]
        self._Y = Y[mask]
        self._V = self._V[mask]


"""

a = np.zeros((20, 100), float)
a.flat[np.random.randint(0, 2000, 100)] = np.random.random(100)
dd = h5dict("bla", in_memory=True)
b = h5dictSparseMatrix(a, dd)
b2 = defaultMatrix(a)


sums = b.getSums()
sums2 = b2.getSums()
print [len(i) for i in sums],
print [len(i) for i in sums2]
print b.getData().sum(),
print b2.getData().sum()
b2.divideByVectors(sums)
b.divideByVectors(sums)
print ((b.getSumX() - b2.getSumX()) ** 2).sum()
print b2.getData().sum()
exit()

"""


class HiResHiC(object):
    """
    Class to store Hi-C data at high resolution,
    perform basic filtering and iterative correction.

    This class provides an efficient storage for the Hi-C data, allowing to retrive a
    heatmap between any chromosomes quickly.
    """

    def __init__(self, genome, resolution, storageFile="inMemory", mode="w"):
        """
        Initializes the high-resolution Hi-C data storage.

        Parameters
        ----------

        genome : folder or Genome object
            matching Genome object or folder to load it form
        resolution : int
            Resolution (number of bp per bin)
        storageFile : str (optional)
            File to store the h5dict.
            File will be created.
            By default stores in memory
        mode : "w", "w-" "r+" or "a", optional
            Access mode to h5dict (see h5dict manual)
        """

        inMemory = (storageFile == "inMemory")

        self._h5dict = h5dict(storageFile, mode=mode, in_memory=inMemory)

        if type(genome) == str:
            genome = Genome(genome, readChrms=["#", "X"])
        assert isinstance(genome, Genome)
        self.genome = genome

        self.resolution = resolution
        self.genome.setResolution(resolution)

        if self.genome.numBins < 7000:
            print("Total number of bins in the genome is just %d" % self.genome.numBins)
            warnings.warn("For low-resolution analysis use binnedData, as it provides"
                          "more analysis tools")

        M = self.genome.chrmCount
        self.cisKeys = [(i, i) for i in range(M)]
        self.transKeys = [(i, j) for i in range(M) for j in range(M) if j > i]
        self.allKeys = self.cisKeys + self.transKeys

        self.data = {}
        self._initChromosomes()

    def _initChromosomes(self):
        "internal: loads mappings from the genome class based on resolution"
        self.chromosomeStarts = self.genome.chrmStartsBinCont
        self.centromerePositions = self.genome.cntrMidsBinCont
        self.chromosomeEnds = self.genome.chrmEndsBinCont
        self.trackLength = self.genome.numBins

        self.chromosomeCount = self.genome.chrmCount
        self.chromosomeIndex = self.genome.chrmIdxBinCont
        self.positionIndex = self.genome.posBinCont
        self.armIndex = self.chromosomeIndex * 2 + \
            np.array(self.positionIndex > self.genome.cntrMids
                     [self.chromosomeIndex], int)

    def _checkConsistency(self):
        "Checks if shapes of all datasets are correct"
        self._hasData()
        for i in self.data:
            shape = self.data[i].shape
            chr1, chr2 = i
            expected = (self.genome.chrmLensBin[chr1],
                        self.genome.chrmLensBin[chr2])
            if (expected[0] != shape[0]) or (expected[1] != shape[1]):
                raise ValueError("Wrong dimensions of a chromosomes {2}: expected {0},"
                                 "got {1}".format(expected, shape, i))

    def _hasData(self):
        mykeys = list(self.data.keys())
        for i in range(self.genome.chrmCount):
            if ((i, i)) not in mykeys:
                print("Not all cis keys found in data!")
                print("Please load data!")
                raise ValueError("Please load data first")

    def _marginalError(self, marginals=None, percentile=99.9):
        """Checks after each pass of IC if marginals are close enough to 1.
        The error is calculated as the specified percentile of deviation
        from the mean marginal.
        """
        if marginals is None:
            if not hasattr(self, "marginals"):
                return 99999
            marginals = self.marginals
        marginals = np.concatenate(marginals)
        marginals = marginals[marginals != 0]
        error = np.percentile(np.abs(marginals - marginals.mean()), percentile)
        return error / marginals.mean()

    def getCombinedMatrix(self, force=False):
        """Combines all matrices into one giant matrix.

        Should be used for testing only, as it will run out of RAM at high resolution.
        """
        if self.genome.numBins > 15000:
            warnings.warn("Resolution is too high. Be careful with RAM!")
        if (self.genome.numBins > 30000) and (force is False):
            raise RuntimeError("Matrix will take more than 8GB RAM. Use force to override")
        t = np.zeros((self.genome.numBins, self.genome.numBins), float)
        for chr1, chr2 in self.data:
            data = self.data[(chr1, chr2)].getData()
            beg1, end1 = self.genome.chrmStartsBinCont[
                chr1], self.genome.chrmEndsBinCont[chr1]
            beg2, end2 = self.genome.chrmStartsBinCont[
                chr2], self.genome.chrmEndsBinCont[chr2]
            t[beg1:end1, beg2:end2] = data
            t[beg2:end2, beg1:end1] = data.T
        return t

    def setRowsToZero(self, rows):
        """Sets select rows to zero

        Parameters
        ----------
        rows : list of arrays/lists
            List of per-chromosome row numbers to be removed.
            Counting starts from 0 for each chromosome.
        """
        rows = [np.array(i) for i in rows]
        for chr1, chr2 in self.data:
            self.data[(chr1, chr2)].clearRows((rows[chr1], rows[chr2]))

    def loadData(self, dictLike,
                 keyFunction=lambda x: ("%d %d" % x),
                 mode="All",
                 cisProtocol=h5dictMatrix,
                 transProtocol=h5dictSparseMatrix):
        """
        Parameters
        ----------

        dictLike : dictionary-like structure or str
            either by-chromosome h5dict, generated by fragmentHiC, or
            any dictionary-like object with by-chromosome heatmaps
        keyFunction : function tuple->string
            Function to convert chromosome pairs (chr1, chr2) to a
            key in a dictLike. Default: "chr1 chr2"
        mode : "all" or "cis"
            Use or not use trans chromosome maps
        cisProtocol : class similar to defaultMatrix
            see below
        transProtocol : class similar to defaultMatirx
            see below


        cisProtocol and transProtocol should implement all functions, currently
        defined in the defaultMatrix protocol. If inhereted from defaultMatrix,
        it should implement proper get and set functions It cannot store the
        matrix itself in memory, and should forget it after any function call.
        """

        if mode.lower() not in ["cis", "all"]:
            raise ValueError("Mode can be only 'cis' or 'all'")

        if type(dictLike) == str:
            try:
                dictLike = h5dict(dictLike, 'r')
            except:
                raise ValueError("Cannot open h5dict at filename %s" %
                                 dictLike)

        for myKey in self.cisKeys:
            try:
                data = dictLike[keyFunction(myKey)]
            except KeyError:
                raise KeyError("Key {0} not found in h5dict".format(
                    keyFunction(myKey)))

            self.data[myKey] = cisProtocol(
                data, dictToSave=self._h5dict, key=myKey)

        if mode.lower() == "all":
            for myKey in self.transKeys:
                try:
                    data = dictLike[keyFunction(myKey)]
                except KeyError:
                    raise KeyError("Key {0} not found in h5dict".format(
                        keyFunction(myKey)))
                self.data[myKey] = transProtocol(
                    data, dictToSave=self._h5dict, key=myKey)

        self._checkConsistency()

    def getMarginals(self, normalizeForIC=False):
        """
        Returns a sum over each row/column, and saves them to self.marginals

        normalizeForIc=True will normalize them to mean 1.
        """
        self._hasData()
        marginals = [np.zeros(i, float) for i in self.genome.chrmLensBin]
        for chr1, chr2 in self.data:
            m2, m1 = self.data[(chr1, chr2)].getSums()
            marginals[chr1] += m1
            if chr1 != chr2:
                marginals[chr2] += m2

        self.marginals = marginals
        return marginals

    def divideByMarginals(self, marginals=None):
        """Divides matrix by a vector.
        If the vector is not provided, will divide it by a
        marginals calculated previously.
        """
        self._hasData()
        if marginals is None:
            marginals = self.marginals

        for chr1, chr2 in self.data:
            m2 = marginals[chr1]
            m1 = marginals[chr2]
            self.data[(chr1, chr2)].divideByVectors((np.sqrt(m1), np.sqrt(m2)))

    def iterativeCorrection(self, tolerance=1e-2):
        """
        Performs iterative correction in place.
        Tolerance is the maximum allowed relative deviation of the marginals.

        Tries to set biases to self.biases after IC.
        """
        self._hasData()
        curPass = 0
        marginals = np.ones(self.genome.numBins, float)
        while self._marginalError() > tolerance:
            m = self.getMarginals(normalizeForIC=True)
            marginals *= np.concatenate(m)
            self.divideByMarginals()
            print("Pass = %d, Error = %lf" % (curPass, self._marginalError()))
            curPass += 1
        self.biases = marginals
        return self.biases

    def removeDiagonal(self, m=1):
        """
        Removes diagonal from the data.
        m=0 is main diagonal, m=1 is main +-1, etc.
        """
        self._hasData()
        for i in self.cisKeys:
            data = self.data[i].getData()
            numutils.removeDiagonals(data, m)
            self.data[i].setData(data)

    def removePoorRegions(self, percent=0.5):
        """
        Removes "percent" percent of regions with low coverage.
        """
        self._hasData()
        marg = self.getMarginals(normalizeForIC=False)
        allMarg = np.concatenate(marg)
        cut = np.percentile(allMarg[allMarg > 0], percent)
        nonZeroMarg = allMarg[allMarg > 0]
        numRemoved = (nonZeroMarg < cut).sum()
        toRemove = [np.nonzero(i < cut)[0] for i in marg]
        self.setRowsToZero(toRemove)
        print("Removed %d poor bins" % numRemoved)

    def export(self, filename):
        mydict = h5dict(filename)
        for i in self.allKeys:
            data = self.data[i].getData()
            mydict["%d %d" % i] = data
        mydict["resolution"] = self.resolution

    def removeMadMax(self, madMax=3.0, adjustByChromMedian=True):
        margs = self.getMarginals()
        masks = hicShared.getMasksMadMax(margs, madMax, adjustByChromMedian)
        self.setRowsToZero([~mask for mask in masks])
        print("Removed {} poor bins".format(sum([(~mask).sum() for mask in masks])))

    def getCisEig(self, 
            robust=True,
            reIC = False, 
            numEigs = 3, 
            byArm=True, 
            verbose=True, 
            sortByGCCorr=False):

        cisEigVecs = {}
        cisEigVals = {}

        for chr1, _ in self.cisKeys:
            if verbose:
                print('Processing chromosome #{}'.format(chr1))

            hm = self.data[(chr1, chr1)]
            GC = self.genome.GCBin[chr1]

            if byArm:
                cent = self.genome.cntrMids[chr1] // self.resolution

                eigvecs_left, eigvals_left = hicShared.cisEigProperNorm(
                    hm.getData()[:cent, :cent], robust=robust, reIC=reIC, numEigs=numEigs, 
                    GC=GC[:cent], sortByGCCorr=sortByGCCorr)
                eigvecs_right, eigvals_right = hicShared.cisEigProperNorm(
                    hm.getData()[cent:, cent:], robust=robust, reIC=reIC, numEigs=numEigs,
                    GC=GC[cent:],sortByGCCorr=sortByGCCorr)

                cisEigVecs[chr1] = np.hstack([eigvecs_left, eigvecs_right])
                cisEigVals[chr1] = np.vstack([eigvals_left, eigvals_right])
            else:
                eigvecs, eigvals = hicShared.cisEigProperNorm(
                    hm.getData(), reIC=reIC, numEigs=numEigs, 
                    GC=GC, sortByGCCorr=sortByGCCorr)
                cisEigVecs[chr1] = eigvecs
                cisEigVals[chr1] = eigvals
        self.cisEigVecs = cisEigVecs
        self.cisEigVals = cisEigVals
        
        return cisEigVecs, cisEigVals

    def getTransEig():
        pass


