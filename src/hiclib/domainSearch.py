import numpy,numexpr 


    


class gradientDescentDomains:
    "rewrite using numexpr module for a much faster operation!" 
    def __init__(self,data):
        self.O = data
        self.N = len(data)
        self.infCount = 0  

    def gradientPoissonModel(self,vec):
        
        N = self.N 
        B = vec[:N]
        D = vec[N:]
        O = self.O
        B1 = B[:,None]
        B2 = B[None,:]
        D1 = D[:,None]
        D2 = D[None,:]
        
        MOE = numexpr.evaluate("exp(B1 + B2 + D1 * D2)")
        B1,B2,D1,D2,O,MOE  #Eclipse warning remover         
        gb = 2 * numexpr.evaluate("sum(O - MOE,axis = 1)")                        
        #print gb1.asarray()
        #print gb
        gc =  2 * (numexpr.evaluate("sum(O * D2 -  MOE * D2,axis = 1)"))        
        #gc =  gc1.asarray().reshape(N) 
        #gb = gb1.asarray().reshape(N)
        t = numpy.r_[gb,gc]
        #cm.shutdown()          
        #print t.shape         
        return -t


        
    def functionPoissonModel(self,vec):
        self.vec = vec 
        N = self.N 
        B = vec[:N]
        D = vec[N:]
        O = self.O
        B1 = B[:,None]
        B2 = B[None,:]
        D1 = D[:,None]
        D2 = D[None,:]
        M0 = numexpr.evaluate("(B1 + B2 + D1 * D2)")    
        O,B1,B2,D1,D2,M0 #Eclipse warning removal     
        #M0 = B[:,None] + B[None,:] + D[:,None] * D[None,:]        
        M1 = numexpr.evaluate("O * M0 - exp(M0)")        
        s =  M1.sum()
        
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
        scipy.optimize.fmin_cg(energy,   numpy.r_[(numpy.random.random( self.N) + 0.5),numpy.random.random( self.N) * 0.2 -0.1],gradient, gtol=10, maxiter = 400)
         
        #g = scipy.optimize.fmin_cg(self.funExp, numpy.random.random(2 * self.N)+0.5,self.gradientExp, gtol=1)
        #g = scipy.optimize.fmin_cg(self.funGeoff, numpy.random.random(2 * self.N)+0.5,self.gradientGeoff, gtol=5,maxiter = 55 )
        #g = scipy.optimize.fmin_cg(self.funGeoffOld, numpy.random.random(2 * self.N)+0.5,self.gradientGeoffOld, gtol=5,maxiter = 500 )
        self.B = self.vec[:self.N]
        self.D = self.vec[self.N:]
        return self.vec
    
    def test(self):
        self.funGeoffOld(numpy.r_[(numpy.random.random( self.N) + 0.5),numpy.random.random( self.N) * 0.2 -0.1])
        self.gradientGeoffOld(numpy.r_[(numpy.random.random( self.N) + 0.5),numpy.random.random( self.N) * 0.2 -0.1])

def exponentialDomains(heatmap,returnBiases = False ):
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
        dom.energyFunction = dom.functionPoissonModel
        dom.gradientFunction = dom.gradientPoissonModel
        try: 
            vec = dom.doSearch()
            break 
        except ArithmeticError:
            pass
    if returnBiases == True: 
        return [vec[dom.N:],vec[:dom.N]] 
    return vec[dom.N:]   
    
    