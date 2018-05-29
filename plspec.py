#!/usr/bin/python
#
# first argument is dump number.
# you must edit the file to change which variable is plotted!
#
# import
import numpy as np
import matplotlib.pyplot as plt
#
# open spectrum file
data = np.loadtxt("spectrum.dat")
tdata = np.transpose( data )
#
lw = tdata[0,:]	# log of photon energy in electron rest-mass units
#
nbins = 6
nspec = len(lw)
#
i = np.arange(0,nbins,1)
#
nvars = 8
nLn = np.log10(tdata[1+i*nvars,:] + 1.e-30)
tauabs = tdata[2+i*nvars,:]
tauscatt = tdata[3+i*nvars,:]
x1av = tdata[4+i*nvars,:]
x2av = tdata[5+i*nvars,:]
x3av = tdata[6+i*nvars,:]
nscatt = tdata[7+i*nvars,:]
#
# normalize 
me = 9.1e-28
c = 3.e10
h = 6.626e-27
lw = lw + np.log10(me*c*c/h)   # convert to Hz from electron rest-mass energy
Lsol = 3.83e33  
nLn = nLn + np.log10(Lsol)  # convert to erg/s from Lsol
#
plt.step(lw, nLn[0])
plt.step(lw, nLn[1])
plt.step(lw, nLn[2])
plt.step(lw, nLn[3])
plt.step(lw, nLn[4])
plt.step(lw, nLn[5])
#
# labels
plt.rc('text', usetex=True) 
plt.rc('font', size = 16)
plt.rc('font', family='times new roman')
plt.xlabel('$\\nu [{\\rm Hz}]$', weight='bold', fontsize=20)
plt.ylabel('$\\nu L_\\nu [{\\rm erg\\,\\,s}^{-1}]$', weight='bold', fontsize=20)
#
minlognu = 9
maxlognu = 22
#
eV = 1.602e-12
mum = 1.e-4 # 1 micron
ang = 1.e-8 # 1 angstrom

def plotenergy( freq, lab ):
	tmplnu = (np.log10(freq) - minlognu)/(maxlognu - minlognu)
	plt.figtext(tmplnu, 0.88, lab, rotation=-90, fontsize=15, va='top')
	return

plotenergy(30.e3*eV/h, '30 keV')
plotenergy(1.e3*eV/h, '1 keV')
plotenergy(c/(5000.*ang), '5000 $\\AA$')
plotenergy(c/(2.*mum), '2 $\\mu$m')
plotenergy(c/(10.*mum), '10 $\\mu$m')
plotenergy(c/0.1, '1 mm')
#
plt.xlim((minlognu, maxlognu))
plt.ylim((33.5-5, 33.5+2))
#
# plot on screen, where you can save to file after viewing plot
#plt.show()
# or, uncomment to save directly...
plt.savefig('tst.pdf')
