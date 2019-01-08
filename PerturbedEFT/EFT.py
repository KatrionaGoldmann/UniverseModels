import numpy as np
import math

class EFT(object):

    # Initial conditions for the universe. 
    def __init__(self, omega_m, omega_phi, H0, a0, t0, iterations, dt):
        self.omega_m = omega_m
        self.omega_phi = omega_phi
        self.dt, self.H0 , self.a0, self.it = dt, H0, a0, iterations
        self.a, self.Hval, self.Hd = [],[],[]
        self.t = np.arange(t0, t0 + iterations*dt, dt) 
        
   # Concordence model: inferred solution for Lambda-CDM.
    def BackgroundTheory(self):
        anext, H = self.a0, self.H0
        for i in range(len(self.t)):
            self.a.append(anext)
            self.Hval.append(H)
            if i > 0: self.Hd.append((H - self.Hval[i-1])/(self.dt))
            else:  self.Hd.append(0.)
            brack = (self.omega_m)*(math.pow(anext,-3))  + self.omega_phi
            H2 = (math.pow(self.H0, 2))*brack
            H = math.sqrt(H2)
            anext = anext*(1 + H*self.dt)
        return self.a, self.Hval, self.Hd

