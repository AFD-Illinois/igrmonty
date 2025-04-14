"""
  
  plot.py

  Plots grmonty ".dat" spectrum files passed as arguments to script.
  See parameters below initial imports for
    (int) THBIN : which bin in elevation to plot. usually in {0,1,2,3,4,5}
    (dbl) nuMin : minimum frequency (cgs) to show 
    (dbl) nuMax : maximum frequency (cgs) to show
    (dbl) nuLnuMinThres : minimum value for nuLnu (cgs) to show when multiplied by max data value
    (dbl) nuLnuMaxThres : maximum value for nuLnu (cgs) to show when multiplied by max data value

$ python plot.py path/to/grmonty/data/*.dat  

"""

import sys
import numpy as np
import matplotlib.pyplot as plt

THBIN = -1 # nb: -1 == "last bin"; THBIN = 1 bin ~ 15 -> 30 degrees.
nuMin = 1.e7
nuMax = 1.e25
nuLnuMaxThres = 1.e1
nuLnuMinThres = 1.e-5

if __name__ == "__main__":

  fnames = [ x for x in sys.argv[1:] if x[-4:] == ".dat" ]

  for fname in fnames:

    print("plotting {0:s}".format(fname))

    dat = np.loadtxt(fname,skiprows=1).T
    lnu = dat[0]
    nuLnu = dat[1:]
    nuLnu = nuLnu[THBIN]

    plt.close('all')
    plt.step(lnu,nuLnu)
    plt.xlim(nuMin,nuMax)
    plt.ylim(nuLnu.max()/1.e5,nuLnu.max()*1.e1)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('nu [cgs]')
    plt.ylabel('nuLnu [cgs]')
    plt.savefig(fname.replace(".dat",".png"))
