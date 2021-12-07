#!/usr/bin/python
"""
    write batch script for running multiple bhoss files
    
    e.g.) > ./write_grmonty_script.py dir=/ccs/home/astrodoo/PRJW/EHT_simulations/EHT_grmonty_files2/sane_a0 startn=1655 endn=2150 dscale=1.1e-14 nsph=10000

        dir: directory that contains data
        startn: start dump number
        endn: end dump number
        nsph: the number of superphoton
        nskip: (optional) if given, it picks up the data with given step interval
"""

import sys
import os, stat

noLog = True

if (len(sys.argv)<6):
    raise ValueError('dir / nsph / dscale / startn / endn / nskip')

for i,iarg in enumerate(sys.argv):
    if ("dir=" in iarg):
        datadir = iarg.split('=')[-1]
    elif ("nsph=" in iarg):
        nsph = int(iarg.split('=')[-1])
    elif ("dscale=" in iarg):
        dscale = float(iarg.split('=')[-1])
    elif ("startn=" in iarg):
        startn = int(iarg.split('=')[-1])
    elif ("endn=" in iarg):
        endn = int(iarg.split('=')[-1])
    elif ("nskip=" in iarg):
        nsample = int(iarg.split('=')[-1])

if (len(sys.argv)==6):
    nskip=5

# find the data 
allfiles = os.listdir(datadir)
allfiles = [file for file in allfiles if '.hdf5' in file]

# find the name of the grmonty files 
GRMHD_files = [file for file in allfiles if 'grmonty' in file]
dummy_GR = GRMHD_files[0].split('_')

GR_fname = dummy_GR[0]
for i in range(1,len(dummy_GR)-1):
    GR_fname += '_'+dummy_GR[i]

#print(GR_fname)

# write the script
fsh  = 'run_grmonty.sh'
ferr = 'run_grmonty.err'

f = open(fsh,'w')
f.write('#!/bin/sh\n')
f.write('export OMP_NUM_THREADS=40\n')

for i in range(startn,endn+1,nskip):
    GRMHD_file = '%s_%d.hdf5'%(GR_fname,i)
    GRMHD_file = os.path.join(datadir,GRMHD_file)

    if (noLog):
        f.write('./grmonty %d %s %e\n'%(nsph, GRMHD_file, dscale))
    else:
        ferr = '%s_%d.err'%(GR_fname,i)
        if (i == startn):
            f.write('./grmonty %d %s %e >& %s \n'%(nsph, GRMHD_file, dscale, ferr))
        else:
            f.write('./grmonty %d %s %e >> %s \n'%(nsph, GRMHD_file, dscale, ferr))

f.close()
os.chmod(fsh,stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
