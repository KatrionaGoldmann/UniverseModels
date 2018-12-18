from ChiSq import *
from ChiSq4params import *
from scipy.optimize import minimize
from EFT import *

c = 3 # Defines the parameters you are fitting: 2= Om, alpha 3= L0, beta or 4 = all four
Om, Ophi, m0 =  0.3, 0.7, 8*math.pi*6.67E-11
a0, H0, t0, dt, it= 0.15, 1/13.7, 1., 0.001, 2000
t = np.arange(t0, t0 + it*dt, dt)
OM0, Lambda0, beta = 1., -1*math.pow(H0, 2)*Ophi*3, 0.
Eft = EFT(Om, Ophi, H0, a0, t0, it, dt)
data = Eft.BackgroundTheory()[0]
chi1 = chiSq(data, t, H0, a0, t0, Om, Ophi, dt, Lambda0, beta)
chi2 = chiSq(data, t, H0, a0, t0, Om, Ophi, dt, Lambda0, beta)
chi4 = ChiSq4params(data, t, H0, a0, t0, Om, Ophi, dt)

if c ==2:
    theoryfun = chi1.function([0., OM0])
    madeupfun = chi1.function([1.5, 3.])
    value = minimize(chi1, [0., OM0], method = 'TNC', tol = 0.000001, bounds = ((-3, 3), (-3, 3)))
    chi1.check()
    # Finding chi over alpha, omega space.
    alphainits, om0inits = np.arange(-3, 3, 0.1), np.arange(-3, 3, 0.1)
    for i in range(len(alphainits)):
        print alphainits[i]
        for j in range(len(om0inits)):
            try: value2 = minimize(chi2, [alphainits[i], om0inits[j]], method = 'TNC')
            except OverflowError: print "couldn't manage ", i, j
    chi2.check()

if c ==3:
    theoryfun = chi1.function([0., Lambda0])
    value = minimize(chi1, [0., Lambda0], method = 'TNC', tol = 0.000001, bounds = ((-3, 3), (-3, 3)))
    chi1.check()
    # Finding chi over beta, Lambda space.
    betainits, L0inits = np.arange(-1, 1, 0.1), np.arange(-1, 1, 0.1)
    for i in range(len(betainits)):
        print betainits[i]
        for j in range(len(L0inits)):
            try:value2 = minimize(chi2, [betainits[i], L0inits[j]], method = 'TNC')
            except OverflowError: print "couldn't manage ", i, j
    chi2.check()

if c ==5:
    theoryfun = chi1.function([0., 0.])
    value = minimize(chi1, [0., 0.], method = 'TNC', tol = 0.000001, bounds = ((-3, 3), (-3, 3)))
    chi1.check()
    # Finding chi over beta, Lambda space.
    lambdainits, G0inits = np.arange(-3, 3, 0.1), np.arange(-3, 3, 0.1)
    for i in range(len(lambdainits)):
        print lambdainits[i]
        for j in range(len(G0inits)):
            try: value2 = minimize(chi2, [lambdainits[i], G0inits[j]], method = 'TNC', bounds=((-5, 5), (-5, 5)))
            except OverflowError as er: print "couldn't manage", i, j
    chi2.check()

if c==4 :
    theoryfun = chi4.function([0., OM0, 0., Lambda0])
    value = minimize(chi4, [0., OM0, 0., Lambda0], method = 'TNC', tol = 0.000001)
    altparams = chi4.check()
    altfun = chi4.function(altparams)

print 'Success?: ', value.success, value.message
if c==2: print 'The minised alpha = ', value.x[0],  ', \t Omega0 = ', value.x[1]
if c==3: print 'The minised beta = ', value.x[0],  ', \t Lambda0 = ', value.x[1]
if c==5: print 'The minised delta = ', value.x[0],  ', \t Gamma0 = ', value.x[1]
if c==4: print 'beta = ', value.x[2],  ', \t Lambda0= ', value.x[3]
print 'The minised chi^2 = ', value.fun
if c==2 or c==3 or c==5: new = chi1.function(value.x)
if c==4: new = chi4.function(value.x)

# Plotting the function
plt.plot(t, new[0], 'r', label = 'Minimised function')
plt.plot(t, madeupfun[0], 'g', label = 'alpha = 0, beta = 0, Lambda = 0.7')
plt.plot(t, data, 'b', label = 'Lambda-CDM')
if c ==4: plt.plot(t, altfun[0])
plt.legend(loc = 2)
plt.title('Expansion of the universe.')
plt.ylabel('a(t)')
plt.xlabel('time (Gyr)')
plt.grid(True)
plt.show()
