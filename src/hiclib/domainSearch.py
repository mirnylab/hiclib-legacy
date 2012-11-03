#(c) 2012 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Maksim Imakaev (imakaev@mit.edu)


import numpy
import numexpr
import mirnylib.numutils


class gradientDescentDomains:
    """
    A strange class to find domains by finning exponential model.
    Domains are fitted to the formula
    O(bserved) = Poisson(B1 * B2 * exp(D1 * D2)),
    where B1 and B2 are biases, D1 and D2 are domain vectors.


    """
    def __init__(self, data):
        self.O = data
        self.N = len(data)
        self.infCount = 0
        self.wrapPower = 1 / 2.

    def wrap(self, x):
        """
        "Wrapping" function that "compresses" the input vector so that
        when we take exp(B1*B2) we don't really overshoot
        """
        return numpy.sign(x) * ((numpy.abs(x) + 1.) ** (self.wrapPower) - 1)

    def unwrap(self, x):
        """
        "Unwrapping" function that reverts wrapping
        """
        return numpy.sign(x) * ((numpy.abs(x) + 1.) **
                                (1. / self.wrapPower) - 1)

    def diffWrap(self, x):
        """
        Wrapping modifies the derivative
        """
        return self.wrapPower * ((numpy.abs(x) + 1) ** (self.wrapPower - 1))

    def gradientPoissonModel(self, vec):
        N = self.N
        wrapDerivative = self.diffWrap(vec)
        vec = self.wrap(vec)

        B = vec[:N]
        D = vec[N:]
        O = self.O
        B1 = B[:, None]
        B2 = B[None, :]
        D1 = D[:, None]
        D2 = D[None, :]

        MOE = numexpr.evaluate("exp(B1 + B2 + D1 * D2)")
        B1, B2, D1, D2, O, MOE  # Eclipse warning remover
        gb = 2 * numexpr.evaluate("sum(O - MOE,axis = 1)")
            #Gradient for biases
        gc = 2 * (numexpr.evaluate("sum(O * D2 -  MOE * D2,axis = 1)")
                  )  # Gradient for domains

        t = numpy.r_[gb, gc]
        t *= wrapDerivative  # Gradient changed becaue of wrapping
        return -t

    def functionPoissonModel(self, vec):

        vec = self.wrap(vec)

        self.vec = vec
        N = self.N
        B = vec[:N]
        D = vec[N:]
        O = self.O
        B1 = B[:, None]
        B2 = B[None, :]
        D1 = D[:, None]
        D2 = D[None, :]
        M0 = numexpr.evaluate("(B1 + B2 + D1 * D2)")
        O, B1, B2, D1, D2, M0  # Eclipse warning removal
        #M0 = B[:,None] + B[None,:] + D[:,None] * D[None,:]
        M1 = numexpr.evaluate("O * M0 - exp(M0)")
        s = mirnylib.numutils.openmpSum(M1)

        print -s
        if numpy.isfinite(s) == False:
            self.infCount += 1
            if self.infCount == 10:
                raise ArithmeticError
        return -s

    def doSearch(self):
        import scipy.optimize
#        vec = numpy.random.random(2 * self.N)+0.5
#        a =  self.fun(vec)
#        vec[-1] = vec[-1] - 0.0001
#        b =  self.fun(vec)
#
#        c =  self.gradient(vec)[-1]
#        print (a-b) , c * 0.0001
#        raise
        energy = self.energyFunction
        gradient = self.gradientFunction
        scipy.optimize.fmin_cg(energy,
                               numpy.r_[(numpy.random.random(self.N) + 0.5),
                                        numpy.random.random(
                                            self.N) * 0.2 - 0.1],
                               gradient, gtol=0.5, maxiter=400)

        #g = scipy.optimize.fmin_cg(self.funExp,
        #numpy.random.random(2 * self.N)+0.5,self.gradientExp, gtol=1)
        #g = scipy.optimize.fmin_cg(self.funGeoff,
        #numpy.random.random(2 * self.N)+0.5,self.gradientGeoff,
        #gtol=5,maxiter = 55 )
        #g = scipy.optimize.fmin_cg(self.funGeoffOld,
        #numpy.random.random(2 * self.N)+0.5,
        #self.gradientGeoffOld, gtol=5,maxiter = 500 )
        self.B = self.vec[:self.N]
        self.D = self.vec[self.N:]
        return self.unwrap(self.vec)

    def test(self):

        init = numpy.r_[(numpy.random.random(
            self.N) + 0.5), numpy.random.random(self.N) * 0.2 - 0.1]
        funValue = self.functionPoissonModel(init)
        gradValue = self.gradientPoissonModel(init)

        init[0] += 0.001
        newFunValue = self.functionPoissonModel(init)

        print newFunValue - funValue, "difference in functions"
        print gradValue[0] * 0.001, "According to the gradient"
        d1 = newFunValue - funValue
        d2 = gradValue[0] * 0.001
        assert (abs(d1 - d2) / abs(d1) < 0.1)

        funValue = newFunValue

        init[-1] += 0.001
        newFunValue = self.functionPoissonModel(init)

        print newFunValue - funValue, "difference in functions"
        print gradValue[-1] * 0.001, "According to the gradient"
        d1 = newFunValue - funValue
        d2 = gradValue[-1] * 0.001
        assert (abs(d1 - d2) / abs(d1) < 0.1)


def exponentialDomains(heatmap, returnBiases=False):
    """Calculates domains using gradient descent algorithm

    Parameters
    ----------
    heatmap : NxN array
        Trans heatmap with faked cis reads
    returnBiases : bool, optional
        Retun a tuple (domains, biases)

    Returns
    -------
    A vector of domains, or a tuple (domains, biases).

    """
    vec = 0

    for _ in xrange(3):
        dom = gradientDescentDomains(heatmap)
        dom.test()
        dom.energyFunction = dom.functionPoissonModel
        dom.gradientFunction = dom.gradientPoissonModel
        try:
            vec = dom.doSearch()
            break
        except ArithmeticError:
            pass
    if returnBiases == True:
        return [vec[dom.N:], vec[:dom.N]]
    return vec[dom.N:]
