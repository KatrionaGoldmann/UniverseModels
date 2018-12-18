import numpy as np
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class EFT(object):

    # Initial conditions for the universe. 
    def __init__(self, omega_m, omega_phi, H0, a0, t0, iterations, dt):
        self.omega_m = omega_m
        self.omega_phi = omega_phi
        self.dt = dt
        self.H0 = H0
        self.a0 = a0
        self.it = iterations
        self.file = file
        self.t = np.arange(t0, t0 + iterations*dt, dt) 
        
   # Concordence model: inferred solution for Lambda-CDM.
    def BackgroundTheory(self):
        avals, Hvals, Hdot= [], [], []
        anext, H = self.a0, self.H0
        for i in range(len(self.t)):
            avals.append(anext)
            Hvals.append(H)
            brack = (self.omega_m)*(math.pow(anext,-3)) + self.omega_phi
            H2 = (math.pow(self.H0, 2))*brack
            H = math.sqrt(H2)
            anext = anext*(1 + H*self.dt)
        return avals, Hvals, Hdot


