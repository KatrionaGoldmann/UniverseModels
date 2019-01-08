import numpy as np
import math
import matplotlib.pylab as plt
import matplotlib.image as mpimg

class chiSq(object):

    def __init__(self, data, t, H0, a0, t0, Om, Ophi, dt, Lambda0, beta):
        self.dat = data
        self.t = t
        self.Lambda0, self.beta = Lambda0, beta
        self.a, self.Hval, self.Hd, self.Omega, self.Lambda= [], [], [], [], []
        self.phi, self.pdot, self.overp, self.ratio = [], [], [], []
        self.anext = a0
        self.H = H0
        self.Gamma = 0.
        self.chivals, self.params1, self.params2 = [], [], []
        self.n, self.p = 0., 3./7.
        self.dt = dt
        self.H0, self.a0, self.t0, self.Om, self.Ophi = H0, a0, t0, Om, Ophi
        self.gam, self.gamdot, self.pidot, self.pi, self.R = [], [], [], [], []

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
        Lambda0, beta = self.Lambda0, 0.
        G0, delta = 1E-5, 0.
        OM = OM0*(math.pow(self.t[i]/self.t0, alpha))
        Gamma = G0*(math.pow(self.t[i]/self.t0, delta))
        OMprime = OM0*(alpha)*(math.pow(self.t[i]/self.t0, alpha - 1))/self.t0
        L = Lambda0*(math.pow(self.t[i]/self.t0, beta))
        term1 = math.pow(self.H0, 2)*self.Om*(math.pow(ai,-3))/(OM)
        term2 = L/(3*OM)
        term3 = Gamma/(3*OM)
        term4 = H*OMprime/OM
        H2 = term1 - term2 + term3 - term4
        H = math.sqrt(abs(H2))
        return H, OM, L, Gamma


    def function(self, params):
        avals, Hvals, Hdot, Omega, Lambda= [], [], [], [], []
        anext, H = self.a0, self.H0
        for i in range(len(self.t)):
            self.a.append(anext)
            self.Hval.append(H)
            val = i
            temp = self.HOmega(anext, H, val, params)
            #temp = self.HLambda(anext, H, val, params)
            #temp = self.HGamma(anext, H, val, params)
            H, OM, L, G = temp[0], temp[1], temp[2], temp[3]
            self.Omega.append(OM)
            self.Lambda.append(L)
            self.gam.append(G)
            self.gamdot.append((G - self.gam[-1])/self.dt)
            atemp = anext*(1 + H*self.dt)
            if i >1: self.Hd.append((H - self.Hval[-1])/(self.dt))
            else: self.Hd.append(0.)
            anext = anext*(1 + H*self.dt)
        return self.a, self.Hval, self.Hd, self.Omega, self.Lambda

    def check(self):
        i = self.chivals.index(min(self.chivals))
        print(self.chivals[i])
        scaledchi, scaledalpha, scaledomo = [], [], []
        for j in range(len(self.chivals)):
            if self.chivals[j] <= 0.5:
                scaledchi.append(self.chivals[j])
                scaledalpha.append(self.params1[j])
                scaledomo.append(self.params2[j])
        plt.scatter(scaledalpha, scaledomo, c=scaledchi)
        plt.colorbar(label = r'$\chi^2$')
        plt.xlabel(r'$\delta$')
        plt.ylabel(r'$\Gamma_0$')
        plt.xlim([-3,3])
        plt.ylim([-3,3])
        plt.title(r'$\chi^2$ for parameters $\delta$ and $\Gamma_0$')
        plt.show()

    def phirout(self):
        p = self.p
        n = self.n
        for i in range(len(self.a)):
            nbert = ((-1*math.pow(self.Hval[i],2) - 2*self.Hd[i])*p - 3*self.Hval[i]*n )*self.dt + n
            nbloom = (-1*(3*math.pow(self.Hval[i],2) + 2*self.Hd[i])*p - 4*self.Hval[i]*n )*self.dt + n
            n = nbloom
            p = n*self.dt + p
            self.pdot.append(n)
            self.phi.append(p)
        k, m2 = 1., 1.
        for i in range(len(self.phi)):
            dp = -m2*(2*math.pow(k/self.a[i], 2)*self.phi[i] + 6*self.Hval[i]*self.pdot[i] + 6*math.pow(self.Hval[i],2)*self.phi[i])
            self.overp.append(dp)
            self.ratio.append(dp/self.phi[i])
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        ax1.plot(self.t, self.Hval, 'b')
        ax1.set_title('\n $\phi$')
        ax1.set_xlabel('Time (Gyr)')
        ax1.set_ylabel('$\phi$')
        ax2.plot(self.t, self.Hd, 'm')
        ax2.set_title('\n Overdensity', size = 12)
        ax2.set_xlabel('Time (Gyr)')
        ax2.set_ylabel(r' $\delta \rho$')
        ax3.plot(self.t, self.phi, 'c')
        ax3.set_title('\n Scale Factor', size = 12)
        ax3.set_xlabel('Time (Gyr)')
        ax3.set_ylabel(r' a(t)')
        plt.suptitle('Einstein-de Sitter Model \n \n', size = 14)
        fig.set_tight_layout(True)
        #plt.savefig('EinsteindeSitterBert.png')
        plt.show()


    def pirout(self):
        pi = 0.
        x = 0.1
        k = 1.
        pival = 0.
        for i in range(len(self.Hval)):
            c, cdot =  self.gam[i]/2, self.gamdot[i]/2
            dR1 = -1*(24*self.Hval[i]+12*self.Hd[i]+2*math.pow(1/self.a[i], 2))*self.phi[i]
            dR2 = -30*self.Hval[i]*self.pdot[i]
            dR3 = -6*(self.pdot[i] - self.pdot[i-1])/self.dt
            dR = dR1 + dR2 + dR3
            Ric1 = -6*(self.Hd[i] + math.pow(self.Hval[i], 2))
            Ric2 = -2*math.pow(k/self.a[i], 2)*self.phi[i]
            Ric3 = 12*(self.Hd[i] + math.pow(self.Hval[i], 2))*self.phi[i]
            Ric4 = 6*((self.pdot[i]-self.pdot[i-1])/self.dt)
            Ric5 = 24*self.Hval[i]*self.pdot[i]
            Ricci = math.pow(self.a[i], -2)*(Ric1 + Ric2 + Ric3 + Ric4 + Ric5)
            self.R.append(Ricci)
            omprime = (self.Omega[i] - self.Omega[i-1])/self.dt
            x1 = (c + cdot + 6*c*self.Hval[i])*self.phi[i] + 3*self.pdot[i]
            x2 = -1*(c*math.pow(k/self.a[i], 2) - (omprime*(self.R[i] - self.R[i-1])/(4*self.dt)) - 3*c*self.Hd[i])*pival/c
            x3 = -1*((cdot + 3*c*self.Hval[i])/c)*x
            x4 = dR*omprime
            x = (x1 + x2 + x3 + x4)*self.dt + x
            self.pidot.append(x)
            pival = (x - self.pidot[-1])/self.dt
            if i> 1:
                self.pi.append(pival)
                pival = (x - self.pidot[-1])/self.dt
            else:
                self.pi.append(0.)
                pival = 0.
        plt.plot(self.t, self.pidot)
        plt.title(r'$\pi$ Field Overdensity')
        plt.ylabel(r'$\pi$')
        plt.xlabel('Time (Gyr)')
        plt.show()

















'''
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
        return H, OM, L, Gamma
'''
