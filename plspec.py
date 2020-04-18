"""
  
$ python plspec.py path/to/spectrum/files

makes a nice plot from hdf5 file.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py

ME   = 9.1093897e-28
CL   = 2.99792458e10
HPL  = 6.6260755e-27
LSUN = 3.827e33

MANY_SPEC = True

if __name__ == "__main__":

  fnamelist = sys.argv[1:]

  for fname in fnamelist:
    print("plotting spectrum for {0:s}".format(fname))

    with h5py.File(fname, "r") as fp:

      nu = None
      nuLnu = None

      # load data
      if "githash" in fp.attrs.keys():
        nu = np.power(10.,fp["output"]["lnu"]) * ME * CL * CL / HPL
        nuLnu = np.array(fp["output"]["nuLnu"]) * LSUN
        if MANY_SPEC:
          nuLnu = nuLnu[:,:,-1]
        else:
          nuLnu = nuLnu[:,-1]
      else:
        nu = np.power(10.,fp["ebin"]) * ME * CL * CL / HPL
        nuLnu = fp["nuLnu"][:,0] * LSUN

      # plot
      plt.close("all")
      ax = plt.subplot(1,1,1)
      print(nuLnu.shape)
      if MANY_SPEC:
        ax.step(nu, nuLnu.sum(axis=0), "k", label="total")
        ax.step(nu, nuLnu[0,:], label="(synch) base")
        ax.step(nu, nuLnu[1,:], label="(synch) once")
        ax.step(nu, nuLnu[2,:], label="(synch) twice")
        ax.step(nu, nuLnu[3,:], label="(synch) > twice")
        ax.step(nu, nuLnu[4,:], label="(brems) base")
        ax.step(nu, nuLnu[5,:], label="(brems) once")
        ax.step(nu, nuLnu[6,:], label="(brems) twice")
        ax.step(nu, nuLnu[7,:], label="(brems) > twice")
      else:
        ax.step(nu, nuLnu, "k", label="total")

      # formatting
      nuLnu_max = nuLnu.max()
      ax.set_xscale("log")
      ax.set_yscale("log")
      ax.set_xlim([1.e8, 1.e24])
      ax.set_ylim([1.e-10 * nuLnu_max, 1.e1 * nuLnu_max])

      ax.set_xlim([1.e9,1.e24])

      ax.set_xlabel(r"$\nu$ (Hz)", fontsize=16)
      ax.set_ylabel(r"$\nu L_\nu$ (erg s$^{-1}$)", fontsize=16)

      # saving
      plt.legend()
      plt.grid()
      plt.savefig(fname.replace(".h5",".png"))
