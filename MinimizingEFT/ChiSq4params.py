import numpy as np
import math
import matplotlib.pylab as plt
import matplotlib.image as mpimg

class ChiSq4params(object):

    def __init__(self, data, t, H0, a0, t0, Om, Ophi, dt):
        self.dat = data
        self.t = t
        self.avals, self.Hvals, self.Hdot, self.Omega, self.Lambda= [], [], [], [], []
        self.avals.append(a0)
        self.anext = a0
        self.H = H0
        self.Hvals.append(H0)
        self.chivals, self.params1, self.params2,  self.params3, self.params4 = [], [], [], [], []
        self.Hdot.append(0.)       
        self.dt = dt
        self.H0, self.a0, self.t0, self.Om, self.Ophi = H0, a0, t0, Om, Ophi

    def __call__(self, params):
        D = []
        ypred = self.function(params)[0]  # predicted value from fitted function. 
        for i in range(len(ypred)):
            arg = (math.pow((ypred[i] - self.dat[i]), 2))/self.dat[i]
            D.append(arg) # error
        chi = sum(D)
        self.chivals.append(chi)
        self.params1.append(params[0])
        self.params2.append(params[1])
        self.params3.append(params[2])
        self.params4.append(params[3])
        return chi

    def Hubble(self, ai, H, i, params):
        alpha, OM0, beta, L0 = params[0], params[1],  params[2], params[3]
        OM = OM0*(math.pow(self.t[i]/self.t0, alpha))
        if OM ==0: OM = -0.00001
        OMprime = OM0*(alpha)*(math.pow(self.t[i]/self.t0, alpha - 1))/self.t0
        L = L0*(math.pow(self.t[i]/self.t0, beta))
        term1 = math.pow(self.H0, 2)*self.Om*(math.pow(ai,-3))/(OM)
        term2 = L/(3*OM)
        term3 = H*OMprime/OM
        H2 = term1 - term2 - term3
        H = math.sqrt(abs(H2))
        return H, OM, L

    def function(self, params):
        avals, Hvals, Hdot, Omega, Lambda= [], [], [], [], []
        anext, H = self.a0, self.H0
        for i in range(len(self.t)):
            avals.append(anext)
            Hvals.append(H)
            val = i
            temp = self.Hubble(anext, H, val, params)
            H, OM, L = temp[0], temp[1], temp[2]
            Omega.append(OM)
            Lambda.append(L)
            atemp = anext*(1 + H*self.dt)
            Hdot.append((H - Hvals[i-1])/self.dt)
            anext = anext*(1 + H*self.dt)
        return avals, Hvals, Hdot, Omega, Lambda

    def check(self):
        val = min(n for n in self.chivals if n!=min(self.chivals))
        i = self.chivals.index(val)
        altparams = [self.params1[i], self.params2[i], self.params3[i], self.params4[i]]
        print altparams, self.chivals[i]
        return altparams
        
    
