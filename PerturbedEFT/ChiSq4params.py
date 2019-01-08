import numpy as np
import math
import matplotlib.pylab as plt
import matplotlib.image as mpimg

class ChiSq4params(object):

    def __init__(self, data, t, H0, a0, t0, Om, Ophi, dt, name):
        self.dat = data
        self.t = t
        self.a, self.Hval, self.Hd, self.chivals = [],[],[],[]
        self.Omega, self.Lambda, self.Gamma, self.Gamdot = [],[],[],[]
        self.anext, self.dt, self.H = a0, dt, H0
        self.p, self.n = 3./7., 0.
        self.params0, self.params2, self.params4 = [],[],[]
        self.params1, self.params3, self.params5 = [],[],[]
        self.phi, self.pdot, self.overp, self.ratio = [],[],[],[]
        self.R, self.pi, self.pidot = [],[],[]
        self.H0, self.a0, self.t0, self.Om, self.Ophi = H0, a0, t0, Om, Ophi
        self.name = name
        self.k, self.m = 1., 1.

    def __call__(self, params):
        D = []
        ypred = self.function(params)[0]  # LambdaCDM values
        for i in range(len(ypred[0:len(self.dat)])):
            arg = (math.pow((ypred[i] - self.dat[i]), 2))/self.dat[i]
            D.append(arg)
        chi = sum(D)
        self.chivals.append(chi)
        self.params0.append(params[0])
        self.params1.append(params[1])
        self.params2.append(params[2])
        self.params3.append(params[3])
        self.params4.append(params[4])
        self.params5.append(params[5])
        return chi

    def Hubble(self, ai, H, i, params):
        alpha, beta, delta = params[0], params[2],  params[4]
        OM0, L0, G0 = params[1], params[3], params[5]
        OM = OM0*(math.pow(self.t[i]/self.t0, alpha))
        if OM == 0.: OM = -0.00001
        OMprime = OM0*(alpha)*(math.pow(self.t[i]/self.t0, alpha - 1))/self.t0
        L = L0*(math.pow(self.t[i]/self.t0, beta))
        Gamma = G0*(math.pow(self.t[i]/self.t0, delta))
        term1 = math.pow(self.H0, 2)*self.Om*(math.pow(ai,-3))/(OM)
        term2 = -1*L/(3*OM) + Gamma/(3*OM) - H*OMprime/OM
        H2 = term1 + term2
        H = math.sqrt(abs(H2))
        return H, OM, L, Gamma

    def function(self, params):
        avals, Hvals, Hdot, Omega, Lambda= [], [], [], [], []
        anext, H = self.a0, self.H0
        for i in range(len(self.t)):
            self.a.append(anext)
            self.Hval.append(H)
            val = i
            temp = self.Hubble(anext, H, val, params)
            H, OM, L, G = temp[0], temp[1], temp[2], temp[3]
            self.Omega.append(OM)
            self.Lambda.append(L)
            self.Gamma.append(G)
            self.Gamdot.append((G - self.Gamma[-1])/self.dt)
            atemp = anext*(1 + H*self.dt)
            if i > 0: self.Hd.append((H - self.Hval[-1])/self.dt)
            else: self.Hd.append(0.)
            anext = anext*(1 + H*self.dt)
        return self.a, self.Hval, self.Hd

    def phirout(self):
        p = -1*self.p
        n = self.n
        for i in range(len(self.a)):
            nbert = ((-1*math.pow(self.Hval[i],2) - 2*self.Hd[i])*p - 3*self.Hval[i]*n )*self.dt + n
            nbloom = (-1*(3*math.pow(self.Hval[i],2) + 2*self.Hd[i])*p - 4*self.Hval[i]*n )*self.dt + n
            n = nbloom
            p = n*self.dt + p
            self.pdot.append(n)
            self.phi.append(p)
        k, m2 = self.k, self.m
        for i in range(len(self.phi)):
            dp = -m2*self.Omega[i]*(2*math.pow(k/self.a[i], 2)*self.phi[i] + 6*self.Hval[i]*self.pdot[i] + 6*math.pow(self.Hval[i],2)*self.phi[i])
            self.overp.append(dp)
            self.ratio.append(dp/  (math.pow(self.H0, 2)*math.pow(self.a[i], -3.)/(8*math.pi*6.67E-11))  )

    def plotting(self, params):
        self.function(params)
        self.phirout()

        fig, ((ax1, ax2, ax4)) = plt.subplots(1, 3)
        ax1.plot(self.t, self.a, 'b')
        ax1.set_title('\n \n Scale Factor \n', size = 13)
        ax1.set_xlabel('Time (Gyr)')
        ax1.set_ylabel('a(t)')
        ax2.plot(self.t, self.phi, 'm')
        ax2.set_title('\n \n Gravitational Potential \n', size = 13)
        ax2.set_xlabel('Time (Gyr)')
        ax2.set_ylabel(r'$\phi$')
        ax4.plot(self.t, self.ratio, 'y')
        ax4.set_title('\n \n Overdensity \n', size = 13)
        ax4.set_xlabel('Time (Gyr)')
        ax4.set_ylabel(r'$\delta$')
        title = self.name + ' model perturbations'
        plt.suptitle(title, size = 14)
        fig.set_tight_layout(True)
        filename = self.name + 'pertsFINALII.png'
        filename.replace("$", "")
        filename.replace("\ ", "")
        plt.savefig(filename.replace(" ", ""))
        plt.show()

    def check(self):
        val = min(n for n in self.chivals if n!=min(self.chivals))
        i = self.chivals.index(val)
        altparams = [self.params0[i], self.params1[i], self.params2[i], self.params3[i], self.params4[i], self.params5[i]]
        print('second lowest chi-squared: ', self.chivals[i], '\nwith ', altparams, '\n')
        return altparams

        '''
            def pirout(self):
                p = self.p
                n = self.n
                pi, x, k, pival = 0., 0.1, self.k, 0.
                for i in range(len(self.Hval)):
                    c, cdot =  self.Gamma[i]/2, self.Gamdot[i]/2
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
                    x2 = -1*(c*math.pow(k/self.a[i], 2)-(omprime*(self.R[i] - self.R[i-1])/(4*self.dt))-3*c*self.Hd[i])*pival/c
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
        '''
