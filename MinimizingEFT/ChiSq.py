import numpy as np
import math
import matplotlib.pylab as plt
import matplotlib.image as mpimg

class chiSq(object):

    def __init__(self, data, t, H0, a0, t0, Om, Ophi, dt, Lambda0, beta):
        self.dat = data
        self.t = t
        self.Lambda0, self.beta = Lambda0, beta
        self.avals, self.Hvals, self.Hdot, self.Omega, self.Lambda= [], [], [], [], []
        self.avals.append(a0)
        self.anext = a0
        self.H = H0
        self.Gamma = 0.
        self.Hvals.append(H0)
        self.chivals, self.params1, self.params2 = [], [], []
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
        return chi

    def HOmega(self, ai, H, i, params):
        alpha, OM0 = params[0], params[1]
        Lambda0, beta = self.Lambda0, self.beta
        OM = OM0*(math.pow(self.t[i]/self.t0, alpha))
        OMprime = OM0*(alpha)*(math.pow(self.t[i]/self.t0, alpha - 1))/self.t0
        L = Lambda0*(math.pow(self.t[i]/self.t0, beta))
        term1 = math.pow(self.H0, 2)*self.Om*(math.pow(ai,-3))/(OM)
        term2 = L/(3*OM)
        term3 = H*OMprime/OM
        term4 = self.Gamma/(3*OM)
        H2 = term1 - term2 - term3 + term4
        H = math.sqrt(abs(H2))
        return H, OM, L

    def HLambda(self, ai, H, i, params):
        beta, L0 = params[0], params[1]
        OM0, alpha = 1., 0.
        OM = OM0*(math.pow(self.t[i]/self.t0, alpha))
        OMprime = OM0*(alpha)*(math.pow(self.t[i]/self.t0, alpha - 1))/self.t0
        L = L0*(math.pow(self.t[i]/self.t0, beta))
        term1 = math.pow(self.H0, 2)*self.Om*(math.pow(ai,-3))/(OM)
        term2 = L/(3*OM)
        term3 = H*OMprime/OM
        term4 = self.Gamma/(3*OM)
        H2 = term1 - term2 - term3 + term4
        H = math.sqrt(abs(H2))
        return H, OM, L

    def HGamma(self, ai, H, i, params):
        delta, G0 = params[0], params[1]
        OM0, alpha, L0, beta = 1., 0., self.Lambda0, self.beta
        OM = OM0*(math.pow(self.t[i]/self.t0, alpha))
        OMprime = OM0*(alpha)*(math.pow(self.t[i]/self.t0, alpha - 1))/self.t0
        L = L0*(math.pow(self.t[i]/self.t0, beta))
        Gamma = G0*(math.pow(self.t[i]/self.t0, delta) )
        term1 = math.pow(self.H0, 2)*self.Om*(math.pow(ai,-3))/(OM)
        term2 = L/(3*OM)
        term3 = H*OMprime/OM
        term4 = Gamma/(3*OM)
        H2 = term1 - term2 - term3 + term4
        H = math.sqrt(abs(H2))
        return H, OM, L

    def function(self, params):
        avals, Hvals, Hdot, Omega, Lambda= [], [], [], [], []
        anext, H = self.a0, self.H0
        for i in range(len(self.t)):
            avals.append(anext)
            Hvals.append(H)
            val = i
              #temp = self.HOmega(anext, H, val, params)
            temp = self.HLambda(anext, H, val, params)
            #temp = self.HOmega(anext, H, val, params)
            H, OM, L = temp[0], temp[1], temp[2]
            Omega.append(OM)
            Lambda.append(L)
            atemp = anext*(1 + H*self.dt)
            Hdot.append((H - Hvals[i-1])/self.dt)
            anext = anext*(1 + H*self.dt)
        return avals, Hvals, Hdot, Omega, Lambda

    def check(self):
        i = self.chivals.index(min(self.chivals))
        print self.chivals[i]
        scaledchi, scaledalpha, scaledomo = [], [], []
        for j in range(len(self.chivals)):
            if self.chivals[j] <= 0.5:
                scaledchi.append(self.chivals[j])
                scaledalpha.append(self.params1[j])
                scaledomo.append(self.params2[j])
        plt.scatter(scaledalpha, scaledomo, c=scaledchi)
        plt.colorbar(label = r'$\chi^2$')
        plt.xlabel(r'$\alpha$')
        plt.ylabel(r'$\Omega_0$')
        plt.xlim([-3,3])
        plt.ylim([-3,3])
        plt.plot(0,1, marker = '*', color = 'm', markersize = 15)
        plt.title(r'$\chi^2$ for parameters $\alpha$ and $\Omega_0$')
        plt.show()
