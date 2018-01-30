import numpy as np
import pickle
import sys

ME   = 9.1093897e-28
CL   = 2.99792458e10
HPL  = 6.6260755e-27
LSUN = 3.827e33

if len(sys.argv) != 2:
  print 'Usage: '
  print '  python package.py [filename]'
  sys.exit()

fnam = sys.argv[1]

data = np.loadtxt(fnam)
nu = 10.**(data[:,0])*ME*CL**2/HPL

nbin = (len(data[0])-1)/7

nuLnu = np.zeros([nbin, len(nu)])
dOmega = np.zeros(nbin)
for n in xrange(nbin):
  nuLnu[n,:] = data[:,8*n+1]*LSUN
  dOmega[n] = data[0,8*n+2]

out = {}
out['nu'] = nu
out['nuLnu'] = nuLnu
out['dOmega'] = dOmega

pickle.dump(out, open('grmonty_out.p', 'wb'))

