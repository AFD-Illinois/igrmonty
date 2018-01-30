from pylab import *
from subprocess import call
import os
import units
cgs = units.get_cgs()
ME = cgs['ME']
CL = cgs['CL']
KB = cgs['KBOL']
HPL = cgs['HPL']
QE = cgs['QE']

os.chdir('../')
call(['python', 'build.py', 'sphere'])
call(['./grmonty', '10000'])
call(['python', 'package.py', 'spectrum.dat'])

os.chdir('test/')
import pickle
data = pickle.load(open('../grmonty_out.p', 'rb'))
nu = data['nu']
dOmega = data['dOmega']
nuLnu = (data['nuLnu']*dOmega[:,None]/(4.*pi)).sum(axis=0)

import matplotlib.pyplot as plt

plt.figure(figsize=(14,10))
ax = plt.subplot(2,1,1)
ax.step(nu, nuLnu, where='mid')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim([1.e8, 1.e24]); ax.set_ylim([1.e23, 1.e41])
ax.set_xlabel('nu (Hz)'); ax.set_ylabel('nuLnu (erg s^-1)')

Ne     = 9.99456e+07
Thetae = 1.
R      = 5.90816458767e+14

dV = 4./3.*pi*R**3
def Jnu(nu):
  Te = Thetae*ME*CL**2/KB
  rel = (1. + 4.4e-10*Te)
  gff = 1.2
  x = HPL*nu/(KB*Te)

  jv = 2**5*pi*QE**6/(3*ME*CL**3)
  jv *= (2*pi/(3*KB*ME))**(1./2.)
  jv *= Te**(-1./2.)*Ne**2
  jv *= exp(-x)*rel*gff
  return jv

nuLnu_sol = nu*Jnu(nu)*dV
ax.plot(nu, nuLnu_sol, color='k')

ax = plt.subplot(2,1,2)
ax.plot(nu, fabs(nuLnu_sol - nuLnu)/nuLnu_sol, color='k')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim([1.e8, 1.e24]); ax.set_ylim([1.e-3, 1.e1])
ax.set_xlabel('nu (Hz)'); ax.set_ylabel('Fractional Error')

plt.savefig('sphere.png')

