import numpy as np
import matplotlib.pyplot as plt

ME   = 9.1093897e-28
CL   = 2.99792458e10
HPL  = 6.6260755e-27
LSUN = 3.827e33

fnam = 'spectrum.dat'

data = np.loadtxt(fnam)
nu = 10.**(data[:,0])*ME*CL**2/HPL

NVAR = 8
nbin = (len(data[0])-1)/NVAR
print nbin

nuLnu = np.zeros([nbin, len(nu)])
for n in xrange(nbin):
  nuLnu[n,:] = data[:,NVAR*n+1]*LSUN

ax = plt.subplot(1,1,1)
for n in xrange(nbin):
  ax.step(nu, nuLnu[n], where='mid')
#ax.step(nu, nuLnu.mean(axis=0), where='mid', color='k')
#ax.step(nu, nuLnu[-1], where='mid', color='k', linewidth=2)
ax.set_xscale('log'); ax.set_yscale('log')
#ax.axvline(1000.*ME*CL**2/HPL, color='k', linestyle='--')
nuLnu_max = nuLnu.max()
ax.set_ylim([1.e-10*nuLnu_max, 1.e1*nuLnu_max])
ax.set_xlim([1.e8, 1.e24])
plt.show()

