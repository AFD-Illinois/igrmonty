from pylab import *
import sys
import os

ME   = 9.1093897e-28
CL   = 2.99792458e10
HPL  = 6.6260755e-27
LSUN = 3.827e33

fnam = sys.argv[1]
NVAR = 8

data = np.loadtxt(fnam)
nu = 10.**(data[:,0])*ME*CL**2/HPL

nbin = (len(data[0])-1)/NVAR

nuLnu = np.zeros([nbin, len(nu)])
for n in xrange(nbin):
  nuLnu[n,:] = data[:,NVAR*n+1]*LSUN

print nbin

tout = open('spec.txt', 'wb')
tout.write('# nu nuLnu(th bins)\n')
Nd = len(nu)
for n in xrange(Nd):
  tout.write('%g ' % nu[n])
  for i in xrange(nbin):
    tout.write('%g ' % nuLnu[i][n])
  tout.write('\n')
tout.close()

