import numpy as np
import math

class chiSq(object):

    def __init__(self, file):
        self.dat = np.genfromtxt(file) 
        # t = dat[:,0], H = dat[:,1] err = dat[:,2]
        
    def __call__(self, params):
        E = []
        ypred = self.function(params)  # predicted value from fitted function. 
        D = np.subtract(self.dat[:,1], ypred) # error between data and theory
        E = np.zeros(shape = (499,499))
        for i in range (0,499):
            E[i,i] = (2.)**2
        Einv = np.linalg.inv(E)
        Dtran = np.transpose(D)
        temp = np.dot(Einv, D)
        result = np.dot(Dtran, temp)
        return result

    # So we can easily change what function we are fitting. 
    def function(self, params):
        G = 6.67e-11
        t0 = 0.01
        time = self.dat[:,0]
        alpha, beta, L0 = params[0], params[1], 1.
        H = []
        rho = 0.3
        Hubble = 1. 
        for i in range(len(time)):
            tratio = time[i]/t0
            H.append(Hubble + 55)
            om = math.pow(tratio, alpha)
            omdot = alpha*math.pow(tratio, alpha - 1)
            Lambda = L0*math.pow(time[i], beta)
            Hubble = -1*omdot/(2*om) + math.sqrt( math.pow(omdot/om, 2)/4  -  (0.3 - Lambda)*G/(3*om) )
        return H

 

        


